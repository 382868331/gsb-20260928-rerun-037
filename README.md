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

## 版本与历史批次（不可变）

分析规则版本记录在包常量 `variant_pipeline.ANALYSIS_VERSION`
（当前 `orientation-forward-v2`，包版本 `0.3.0`）。**每次分析都发布为一个
独立、不可变的批次**，即使 `--outdir` 与上次完全相同也不会覆盖旧结果：

- 批次目录：`<outdir>/batches/<UTC时间戳>_<内容ID前12位>/`，内含
  `candidates.json/.tsv`、`coverage.tsv`、`report.txt`、`manifest.json`，
  以及 `input/` 下当次使用的 `config.json`、`reference.fasta`、`reads.sam`
  原始快照。输入、解释规则版本（`analysis_version`）与输出在同一批次内
  一一对应。
- `manifest.json` 记录全部文件的 SHA-256、输入哈希、参数、规则版本与
  `content_id`（规则版本+三份输入的联合哈希）；`content_id` 相同的重复
  请求会生成**新批次**并通过 `duplicate_of` 指向前序批次，旧批次原封不动。
- 发布是原子操作：先写隐藏暂存目录，最后整体 rename；发布失败只清理
  暂存目录，`index.json` 也以原子替换更新，已发布版本不受影响。
- 根目录 `index.json` 按时间顺序登记全部批次，进程重开后仍可列出、
  比较、导出和校验任意历史版本。

### 跨版本差异比较

`diff` 只读取两个已冻结批次的 `candidates.json` 与 `manifest.json`，
**绝不重新分析旧输入**，因此旧内容不会被新代码改写。输出列出候选的
新增（appeared）、消失（disappeared）与支持/反对数量、可信度、状态的
逐字段变化，并依据两批次的输入哈希与 `analysis_version` 归因：

- 仅输入哈希不同 → `input_change`（输入变化，如新读段、参考或参数）；
- 仅规则版本不同 → `interpretation_rule_change`（方向解释规则变化）；
- 两者都变 → 同时列出两个归因；相同重复批次 →
  `unchanged_inputs_and_rules`。

两批次的支持记录（包括共享片段合并结果与软截断位置的原始证据）始终
各自随本批次 `input/reads.sam` 与候选字段独立保留，差异比较只做并列
对照，**不拼接两个版本的证据**。

### 历史导出与完整性核对

`export` 直接复制冻结批次文件（含输入快照），不重算；`verify` 只重算
存储文件的 SHA-256 并与 manifest 比对，用于核对旧内容有没有被再次分析
或被外部改写，校验本身也不会触发分析。

## 输出

每个候选包含支持证据（ALT 片段数、正/反向、ALT 比例、平均碱基质量、
低质量 ALT 片段数）与反对证据（REF 片段数、其他等位片段数、重叠冲突片段数、
方向偏差 p 值），以及深度、标记、评分分量和状态
（`PASS` / `SUSPECT` / `REJECTED` / `LOW_COVERAGE`）。

- `candidates.json`：完整结构化结果（含稳定变异身份，供注释阶段 join）
- `candidates.tsv`：候选表
- `coverage.tsv`：逐位片段深度
- `report.txt`：中文可读报告，末尾单列覆盖不足区间

## 使用

```bash
.venv/bin/python examples/generate_example_reads.py
.venv/bin/python -m variant_pipeline call \
  --config examples/config.example.json --outdir out
.venv/bin/python -m variant_pipeline list --outdir out
.venv/bin/python -m variant_pipeline diff \
  --outdir out --from <旧批次id或前缀> --to <新批次id或前缀>
.venv/bin/python -m variant_pipeline export \
  --outdir out --batch <批次id或latest> --dest /path/to/export
.venv/bin/python -m variant_pipeline verify --outdir out --batch latest
```

`list` 与 `diff` 支持 `--json`；批次 id 支持唯一前缀或 `latest`。
对同一输出根目录重复执行 `call` 是安全的：新分析进入新批次目录，旧批次
与 `index.json` 历史记录保持可用。

## 测试

```bash
.venv/bin/python -m pytest -q
```

测试用已知替换、短插入、短缺失、方向偏差、低质量证据、重叠读段合并/冲突
以及覆盖不足的合成数据进行核对；方向约定另由独立手写的标准 SAM 记录核对。

## 免责声明

Research use only. 输出是未经临床验证的候选汇总，用于研究流程开发，
不构成医疗诊断、临床报告或专业遗传解释的替代品。
