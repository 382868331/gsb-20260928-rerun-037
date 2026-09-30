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
（当前 `orientation-forward-v2`，包版本 `0.3.0`）。

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
