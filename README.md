# 基因变异分析系统 — Python 基础流程（候选汇总与证据检查）

接收人工参考序列（FASTA）与已比对读段（SAM），汇总候选变异并记录每个
候选的支持与反对证据。**仅供研究流程使用，不产生任何医疗结论或临床报告。**

## 设计要点

- **统一方向约定（SAM 标准）**：`SEQ`/`QUAL` 在 SAM 中始终按正向参考链存储，
  反向比对只由 FLAG `0x10` 标记，CIGAR 按参考坐标递增排列。解析时直接用
  存储的碱基与参考比较，**不再二次反向互补**——此前反向完全匹配读段被再次
  反向互补，会把与参考一致的读段误报为多个变异。缺口（`D/N`）、软截断
  （`S`，不产生证据）和配对重叠对正、反向读段都按同一套参考坐标解释。
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

## 输出

每次分析在载荷与报告中记录 `rules_version`（见
`variant_pipeline/__init__.py:RULES_VERSION`）。旧批次保留不动；新规则结果请写入
新的输出目录（如 `out-v2/`），不要覆盖旧目录。

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
```

## 测试

```bash
.venv/bin/python -m pytest -q
```

测试用已知替换、短插入、短缺失、方向偏差、低质量证据、重叠读段合并/冲突
以及覆盖不足的合成数据进行核对。

## 免责声明

Research use only. 输出是未经临床验证的候选汇总，用于研究流程开发，
不构成医疗诊断、临床报告或专业遗传解释的替代品。
