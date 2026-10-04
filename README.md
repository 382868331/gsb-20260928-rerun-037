# 基因变异分析系统 — Python 基础流程（候选汇总与证据检查）

接收人工参考序列（FASTA）与已比对读段（SAM），汇总候选变异并记录每个
候选的支持与反对证据。**仅供研究流程使用，不产生任何医疗结论或临床报告。**

## 设计要点

- **统一方向与参考坐标**：所有读段一律按参考正向轴解释。`POS` 是相对
  参考 5'→3' 的起点，CIGAR 按参考坐标正序遍历，SEQ/QUAL 直接与参考
  逐位比较；FLAG 0x10 仅用于方向统计，绝不再对 SEQ 做反向互补。因此
  一条与参考完全一致（如 AACG）的反向读段产生 **零** 个变异，与正向
  读段完全相同。缺口（D/N）、软/硬截断（S/H/P）与配对重叠位置在正、
  反向读段中均落在同一参考坐标；截断碱基不产生证据。
- **配置与代码分离**：样本、参考版本（assembly）与全部分析参数集中在
  JSON 配置文件（见 `examples/config.example.json`），不写入临时脚本。
- **稳定变异身份**：等位基因先做去冗余与左移归一化（`variant_pipeline/normalize.py`），
  再生成规范键 `(chrom, pos, ref, alt)` 与 HGVS 风格 ID；后续注释功能沿用同一身份。
- **按片段而非读次数计票**：成对读段的重叠区域先在片段（fragment）级合并，
  同一分子只投一票；重叠 mates 等位冲突时该位置弃权并记为反对证据。
- **多维可信度，而非出现次数**：可信度综合 ALT 比例、碱基质量、方向平衡、
  覆盖深度、低质量证据比例和冲突证据（权重见
  `variant_pipeline/evidence.py:SCORE_WEIGHTS`）。
- **方向偏差检验**：ALT/REF × 正/反链 2×2 表做双侧 Fisher 精确检验。
- **低质量位置**：低于 `min_base_quality` 的 ALT 证据单独计数并降低可信度，
  不计入有效支持。
- **覆盖不足单独说明**：低于 `low_coverage_depth` 的区间在报告与 JSON 中
  以 `low_coverage_regions` 单独列出；候选状态可为 `LOW_COVERAGE`。

## 方向约定核对（独立标准记录）

`tests/test_orientation_contract.py` 中的 SAM 记录全部独立手写，不复用
`examples/generate_example_reads.py`，避免样例生成器沿用分析端的错误
约定。它们用已知参考逐项核对：反向完全匹配（AACG 零变异）、正/反向
覆盖一致、反向 SNV/插入/缺失的参考坐标、正反向读段两端软截断、以及
配对重叠的单片段计票与冲突弃权（保留 `MATE_OVERLAP_CONFLICT` 弃权
原因）。

## 版本与历史批次

分析规则版本记录在包常量 `variant_pipeline.ANALYSIS_VERSION`
（当前 `linked-evidence-units-v5`，包版本 `0.5.0`）。

`--outdir` 是**版本化结果库**的根目录：每次 `call` 都把结果发布为一个
新的不可变版本目录 `outdir/runs/vNNNN/`（含 `manifest.json`），绝不覆盖
或重算旧版本，因此同一输出位置反复分析新输入也是安全的。发布经临时
目录 + 原子重命名完成：发布失败不会留下半成品版本，已有版本不受影响。
完全相同的请求（输入文件、配置、规则版本的指纹一致）不会重复分析，
直接指向已发布版本；需要强制再分析时用 `--force`。

每个版本的 manifest 记录输入指纹（配置、reads、参考 FASTA 的 SHA-256）、
规则版本和全部输出文件的哈希，输入、解释规则与输出一一对应：

- `variant-pipeline list --outdir OUT` 列出版本；
- `variant-pipeline verify --outdir OUT [--version vNNNN]` 重算哈希，
  核对历史版本未被再次分析改写；
- `variant-pipeline export --outdir OUT --version vNNNN --dest DIR`
  在完整性核对通过后导出任一历史版本（重开进程后同样可用）；
- `variant-pipeline compare --outdir OUT --from vA --to vB`
  比较两个版本的候选差异。

## 版本间差异比较

`compare` 按稳定变异身份对齐两个版本的候选，区分**新增 / 消失 /
支持数量或状态变化**，并根据 manifest 归因：输入（reads/参考）变化、
方向解释规则（analysis_version/参数）变化，或两者同时变化。对两个版本
都存在的候选，报告把各自版本的支持/反对证据（含共享片段与软截断
留下的原始记录）并排列出，**绝不跨版本合并或拼接支持计数**。

## 原始片段来源与共享歧义

每个片段的证据身份为 `来源文件SHA-256前缀::QNAME`。因此两个批次即使使用相同
QNAME，也不会被误认为同一份片段；同一批片段中途保存后重复导入，会按该身份去重，
不增加任何独立支持，也不会覆盖此前保存在版本库外部的人工选择。

候选输出同时给出：
- `independent_alt_fragments`：只支持该候选的独立片段数；
- `shared_ambiguous_alt_fragments`：同一原始片段同时支持互斥候选的共享歧义数；
- `ambiguity_groups`：共享片段、候选集合及来源记录；
- `fragment_sources`：QNAME、来源 ID/名称和批次 ID，供复核连接关系。

`call --extra-reads 批次2.sam` 可在一次分析中加入新批次。比较版本时，
`ambiguity_tracking` 将旧/新批次中的共享片段分成已解决、仍保留和新引入三类；
比较结果还分别保留两个版本自己的临时采用决定，不把规则版本的决定相互覆盖。

临时采用不修改不可变分析版本：
- `variant-pipeline select --outdir OUT --version v0001 --variant ID` 保存选择；
- `variant-pipeline selected-view --outdir OUT --version v0001` 查看该选择下的证据视图；
- `variant-pipeline withdraw --outdir OUT --version v0001` 撤回选择。

选择记录位于 `selections/`，原 `runs/vNNNN` 输出哈希保持不变。导出版本时，
除原有输出和 manifest 外还包含 `shared_evidence_audit.json`，其中保留足以复核
共享关系的输入哈希、片段来源和歧义分组记录。

## 同一原分子的关联片段

关联关系是一份独立 JSON（可在配置中用 `associations_path` 指定，或调用
`call --associations relations.json`）。每条关系显式记录两端的证据 ID、
两端来源、分子身份、关系 ID/版本、方向和距离；证据 ID 含来源哈希前缀，
因此同名但来源或分子身份不同的片段不会自动配对。端点按（来源，片段标识）
共同定位：`a`/`b` 可以写来源限定 ID（`来源::QNAME`），也可以只写 QNAME，
此时在该端声明的来源内解析；不同来源的同名片段各自保留，交换导入顺序
不会改变关系状态、有效组合或计数，只有同一（来源，QNAME）的重复导入
才会被折叠并记为 `duplicate_import`：

```json
{
  "relations": [{
    "association_id": "unit-0001",
    "molecule_id": "original-molecule-1",
    "a": "<source-hash>::read-a",
    "b": "<source-hash>::read-b",
    "source_id_a": "<source-hash>",
    "source_id_b": "<source-hash>",
    "orientation": "opposite",
    "min_distance": 80,
    "max_distance": 500,
    "relation_id": "mate-rule-2026",
    "relation_version": "1"
  }]
}
```

输出的 `evidence_units` 同时保留：
- 每个端片段单独可支持的候选（`single_fragment_candidates`）；
- 经方向、距离和来源核对后仍成立的组合（`valid_combinations`、`linked_candidates`）；
- 被组合关系排除的候选（`excluded_candidates`）；
- 一端迟到（`pending_endpoint`）、片段匹配不足（`match_insufficient`）和
  关联依据矛盾（`relation_inconsistent`），不会静默丢弃任一端；
- **分子内部等位冲突（`molecule_conflict`）**：同一原分子的两端在同一位点
  （参考区间重叠）给出**不同等位**时，如一端 T>A、另一端 T>G，或重叠的插入与
  删除，双方候选**都保留**，单元标记为 `molecule_conflict` 而不是 `resolved`：
  不删除任一端、不按导入顺序选答案，该位点由该分子弃权一次（双方候选各记一次
  `conflicting_overlap_fragments`，均不计 ALT 支持），冲突依据保存在
  `conflicting_combinations`（两端等位、重叠坐标），不满足方向/距离的组合保存在
  `excluded_combinations`。不重叠、兼容的候选仍可正常成立；两端在同一位点给出
  相同等位（纯合）不算冲突。
- 重复导入、多个拆分片段和部分缺失均归入同一个可说明证据单元。

一个证据单元可以支持多个不同位置的候选，但对任一候选最多形成一份独立支持，
不会因为有多个拆分端而重复计票；一端迟到时整个关系仍只算一个分子，迟到端到达
后若与另一端冲突则按分子内部冲突处理。撤回关系时，原单片段歧义恢复（每个端
重新成为独立片段），而不是新增证据。

人工采用可指定证据单元和规则：`select --evidence-unit unit-0001
--relation-id mate-rule-2026 --rationale ...`。再次选择会把旧选择完整放入
`history`，不覆盖旧选择；在冲突单元上交换采用（先选一端等位、再改选另一端）
时，冲突依据（`conflicting_combinations`）、被排除组合与支持计数保持一致，
未被选中的冲突等位标记为 `intra_molecule_conflict`，不会被当作按导入顺序排除。
撤回后导出仍保留采用记录，审计可从支持数量追到原分子、端片段、关系 ID/版本、
冲突依据以及排除依据。

## 输出

每个候选包含支持证据（ALT 片段数、正/反向、ALT 比例、平均碱基质量、
低质量 ALT 片段数）与反对证据（REF 片段数、其他等位片段数、重叠冲突片段数、
方向偏差 p 值），以及深度、标记、评分分量和状态
（`PASS` / `SUSPECT` / `REJECTED` / `LOW_COVERAGE`）。

每个版本目录 `runs/vNNNN/` 内含：
- `candidates.json`：完整结构化结果（含稳定变异身份，供注释阶段 join）
- `candidates.tsv`：候选表
- `coverage.tsv`：逐位片段深度
- `report.txt`：中文可读报告，末尾单列覆盖不足区间
- `manifest.json`：输入指纹、规则版本与输出哈希（完整性核对依据）

## 使用

```bash
.venv/bin/python examples/generate_example_reads.py
.venv/bin/python -m variant_pipeline call \
  --config examples/config.example.json --outdir out
.venv/bin/python -m variant_pipeline list --outdir out
.venv/bin/python -m variant_pipeline compare --outdir out --from v0001 --to v0002
.venv/bin/python -m variant_pipeline select --outdir out --version v0001 --variant synth-1:g.15G>A
.venv/bin/python -m variant_pipeline selected-view --outdir out --version v0001
.venv/bin/python -m variant_pipeline withdraw --outdir out --version v0001
.venv/bin/python -m variant_pipeline verify --outdir out
.venv/bin/python -m variant_pipeline export --outdir out --version v0001 --dest out-v1
```

## 测试

```bash
.venv/bin/python -m pytest -q
```

测试用已知替换、短插入、短缺失、方向偏差、低质量证据、重叠读段合并/冲突
以及覆盖不足的合成数据进行核对；方向约定另由独立手写的标准 SAM 记录核对。

## 免责声明

Research use only. 输出是未经临床验证的候选汇总，用于研究流程开发，
不构成医疗诊断、临床报告或专业遗传解释的替代品。
