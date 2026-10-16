# Cardiometabolic-inflammatory-burden-and-CVD-mortality
Its incremaental predictive value beyond age and sex was uncertain,its specificity for CVD over cancer death was not established, and molecular and genetic evidence supported some components but not a causal effect of the composite
[README.md](https://github.com/user-attachments/files/33548065/README.md)
# CMIBI 与死亡结局重分析项目

更新日期：2026 年 10 月 16 日。

本项目用于核查和重分析 NHANES 心血管代谢炎症负担指数（CMIBI）与死亡的关联，并提供历史队列验证、缺失协变量多重插补、转录组分析及成分层面的孟德尔随机化（MR）结果。最新完整报告包含 **52 页、36 张表、25 幅图**，已加入新稿要求的补充图 S1。

本 README 说明交付文件与复现方法，不替代完整统计报告。文件中的结果来自实际计算；未根据显著性筛选插补结果或遗传工具组合。

## 1. 文件入口

本节路径以本机交付目录 `outputs/` 为基准；若只解压代码附件，完整 Word 报告与独立投稿图文件需另行下载。附件内部的结果可在 `results/` 找到。

| 文件或目录 | 内容 |
| --- | --- |
| `CMIBI_完整重分析报告.docx` | 完整统计报告、主表、图形、方法、局限与参考文献核查；建议首先阅读 |
| `CMIBI_代码数据与图表附件.zip` | 代码、结果、部分分析对象与来源证据的打包备份 |
| `CMIBI_复现附件/` | ZIP 的对应文件夹；目录结构见下文 |
| `final_analysis/` | 本机最新分析结果，包括 CSV、图形和核查记录 |
| `Figure_S1_CMIBI_CVD_spline.png` | 补充图 S1，600 dpi PNG |
| `Figure_S1_CMIBI_CVD_spline.pdf` | 补充图 S1，矢量 PDF |
| `Figure_S1_legend.txt` | 补充图 S1 的完整英文图注 |
| `重分析附件说明.md` | 附件阅读与历史更新说明 |

早期审核报告、旧英文替换段落和旧 ZIP 保留用于追溯，不应优先于最新完整报告。新稿的补充材料编号与报告图表编号不同，使用时请按图注对应。

## 2. 复现附件结构

```text
CMIBI_复现附件/
├── README.md
├── results/       # 汇总表、逐工具/逐基因结果、图形和核查输出
├── objects/       # 开发数据、插补、历史验证、转录组等分析对象
├── code/
│   ├── scripts/   # 最新新增和更新后的 R/Node 分析脚本
│   └── ...        # 主重分析程序、下载工具及早期代码备份
├── mr_inputs/     # 候选工具、LD clumping 结果、查询表和日志
└── sources/       # PubMed XML、GWAS 元数据、来源核对及上游方法代码
```

同名脚本若同时存在于 `code/` 和 `code/scripts/`，新增分析以 `code/scripts/` 中的最新版本为准；主程序为 `code/reanalyse_cmibi.R`。参考文献核查的 Python/Node 辅助程序保存在 `sources/`。

附件没有重复打包全部原始 NHANES、GEO、大型 GWAS 文件、完整 LD 参考或本地 R 包库。它是代码和分析结果备份，不是可以在任意电脑直接双击运行的安装包。

## 3. 数据与分析范围

- **开发队列：** NHANES 2005–2006、2007–2008、2009–2010、2015–2016、2017–2018；分数完整且符合死亡链接等条件者 10,868 人，M3 完整病例 9,794 人。
- **历史验证：** NHANES 1999–2000、2001–2002、2003–2004；合格样本 5,199 人，M3 完整病例 4,155 人；已核查与开发队列无 SEQN 重叠。
- **随访：** 使用体检起点的 `PERMTH_EXM/12`，公开死亡链接截至 2019 年末；分析全因、CVD 和癌症死亡。
- **缺失处理：** 两队列分别进行 30 份、每份 20 次迭代的协变量插补，不插补死亡、随访或 CMIBI 组成项。
- **转录组：** GSE98895 和 GSE9006，使用提交的表达矩阵，进行尺度检查、适当转换、limma 差异表达及 camera Hallmark 分析；不是从原始 CEL/IDAT 重新处理。
- **MR：** 实际 LD clumping 后的单变量 MR、方向性、多效性和留一法分析，以及 HDL/TG/LDL 联合多变量 MR；未完成 ApoB 条件模型。

历史队列使用开发队列固定的 CMIBI 构建参数和标准差。十年预测模型在 2005–2008 子样本开发后锁定系数，再应用于历史样本，不在验证集中重新拟合开发系数。

## 4. 本机运行环境

| 项目 | 本次使用 |
| --- | --- |
| 操作系统与命令行 | Windows，PowerShell |
| 项目工作目录 | `C:/Users/HUAWEI/Documents/Codex/2026-10-15/new-chat` |
| 原始数据目录 | `E:/NHANES` |
| R | `C:/Program Files/R/R-4.5.1/bin/Rscript.exe` |
| 本地 R 包库 | 项目下 `work/Rlib` |
| 主要包 | survey、survival、mice、limma、statmod、data.table、jsonlite、foreign、ggplot2、ggrepel、officer、flextable |
| MR 工具 | PLINK 1.9；1000 Genomes EUR LD 参考 |
| 辅助工具 | 支持 ESM 的 Node.js；Python 3；Word 渲染使用 LibreOffice 和 Poppler |

部分已记录包版本见 `results/software_versions.csv`。报告中的主环境为 R 4.5.1、survey 4.5、survival 3.8-3、limma 3.64.3。这些是本次使用的版本，不代表必须升级到其他版本。

Windows 下含中文路径的 R 脚本应先设置合适的 UTF-8 locale。本项目常用的脚本开头为：

```r
.libPaths(c("work/Rlib", .libPaths()))
Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
options(survey.lonely.psu = "adjust")
```

## 5. 复现顺序

以下路径以原项目工作目录为基准。在新位置复现时，先恢复原始数据与依赖，并检查各脚本内的绝对路径；将最新脚本放回 `work/scripts/`，下载辅助工具放回 `work/`，主程序放回 `outputs/`。不能直接在压缩包中运行。

1. **主队列核验和模型。** 运行 `outputs/reanalyse_cmibi.R`，生成分数参数、主模型与 `work/audit_models.rds`。按需要运行 `verify_sources.R`、`additional_checks.R`、`prediction_check.R` 和 `discrimination_check.R`。
2. **开发队列插补。** 运行 `work/scripts/imputation.R`，随后运行 `mi_diagnostics.R`。30×20 插补可能耗时十余分钟，具体取决于机器和依赖环境。
3. **历史验证。** 通过 `work/download_inputs.mjs` 补足缺少的官方输入，运行 `external_validation.R`；随后运行 `historical_imputation.R` 和 `historical_mi_diagnostics.R`。
4. **新增统计。** 运行 `new_requests.R`、`new_figures.R` 和 `verify_new.R`，生成随访 IQR、分类 Table 1、三个结局的四分位模型、正式死因差异比较，以及两队列正确分开的缺失情况表。
5. **转录组。** 准备 GEO 表达矩阵、平台注释和 Hallmark 集合，运行 `transcriptomics.R`；最后运行 `publication_figures.R` 完善火山图标签布局。
6. **MR。** 准备 GWAS、PLINK 和 EUR LD 文件，依次运行 `prepare_mr.mjs`、`clump_mr.R`、`extract_mr_outcomes.mjs`、`extract_eur.mjs`、`analyse_mr.R`、`multivariable_mr.R` 和 `mr_supplement_figures.R`。血压单位核查另运行 `check_bp_units.R`，所需官方 transformed 文件路径见脚本。
7. **其他图形和补充图 S1。** 运行 `report_figures.R`、`final_extra_figures.R`、`flow_diagram.R` 等对应绘图脚本；运行 `figure_S1.R` 生成 S1、逐点预测、结点和模型检验记录。
8. **参考文献核查。** 辅助程序顺序为 `verify_references.mjs`、`refine_references.mjs`、`format_references.py`、`verify_pubmed_and_units.py`；下载的 NCBI XML 保存在 `sources/`，核查表位于 `results/reference_audit.csv`。该核查针对最初稿件的 43 条参考文献，不等同于已逐条核对新稿重新编号后的所有引文。
9. **报告与质量检查。** 运行 `verify_final.R` 和 `verify_new.R`，再运行 `build_report.R`。报告构建会调用 `report_new_requests.R` 与 `report_S1.R`。构建完成后运行 `sanitize_report.py`，再渲染并检查所有修改页面。报告渲染工具的本机路径需单独恢复。

原项目中的 PowerShell 调用示例：

```powershell
Set-Location 'C:/Users/HUAWEI/Documents/Codex/2026-10-15/new-chat'
& 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe' 'work/scripts/figure_S1.R'
```

此例仅重绘 S1，前提是已有 `work/audit_models.rds` 和所需 R 包；它不是从零重建全部分析的命令。

## 6. 常用结果文件

下列文件名相对于附件的 `results/` 或本机 `final_analysis/`：

| 文件 | 用途 |
| --- | --- |
| `followup_IQR.csv` | 两队列及 M3 完整病例的随访中位数、IQR 与加权分位数 |
| `table1_categorical_quartiles.csv` | 四分位分类基线表：未加权 n、非缺失分母中的抽样加权百分比 |
| `quartile_definitions.csv` | 固定四分位界点及趋势检验所用组中位数 |
| `M3_quartiles_all_outcomes.csv` | 三个结局 Q2–Q4 vs Q1 的 HR、CI、P 和 P-trend |
| `CVD_cancer_formal_comparison.csv` | 联合死因模型中的 CVD/癌症 HR 比值与正式检验 |
| `multiple_imputation_cox.csv` | 开发队列 Rubin 合并结果 |
| `historical_multiple_imputation_cox.csv` | 历史队列 Rubin 合并结果 |
| `missingness.csv`、`historical_missingness.csv` | 两队列分别记录的缺失情况 |
| `MR_GWAS_effect_units.csv` | 各暴露效应单位、解释与来源 |
| `BP_unit_source_verification.csv` | 本地血压 β、SE 与官方 transformed 文件的核查汇总 |
| `mr_results_corrected.csv`、`mr_diagnostics.csv` | 单变量 MR 及诊断 |
| `multivariable_mr.csv` | HDL、TG、LDL 联合模型结果 |
| `reference_audit.csv`、`references_Vancouver.txt` | 原稿逐条文献核对与工作格式引文 |
| `Figure_S1_curve_data.csv`、`Figure_S1_knots.csv` | S1 曲线的逐点 HR/CI 及结点 |
| `Figure_S1_model_details.txt` | S1 的样本量、参考值、自由度及检验结果 |
| `verification_checks.txt`、`new_verification_checks.txt` | 程序性质量核查记录 |

开发队列按固定的未加权 CMIBI 四分位分组，每组 2,717 人。分类 Table 1 中的 n 未加权，百分比抽样加权且排除该变量缺失者；缺失行单列人数。PIR 描述性分组为 <1、1–<2、≥2，M3 仍使用连续 PIR。

## 7. 补充图 S1

S1 对应新稿的 **Natural spline of CMIBI and CVD death**，已加入完整报告第 52 页，报告内图号为 25。

- 直接使用既有调查加权 M3 自然立方样条模型，没有为绘图另拟一个模型。
- 完整病例 n=9,794，CVD 死亡 312 例。
- 四结点位于开发合格队列 CMIBI 的第 5、35、65、95 百分位。
- 以第 10 百分位（−1.32879 个开发队列 SD）为 HR=1 参考。
- 图示范围为第 1–99 百分位；阴影为逐点 95% CI，不是同时置信带。
- CI 使用 77 个设计自由度；原模型 Wald F 总体与非线性检验的分母自由度为 60。
- 总体 P=0.035327，非线性 P=0.30053；未检出非线性不证明严格线性。

投稿可使用 600 dpi PNG 或矢量 PDF，并使用 `Figure_S1_legend.txt` 中的完整图注。

## 8. 结果解释边界

1. CVD 与癌症的正式 HR 比值为 1.263（95% CI 0.955–1.669），P=0.100，不能因为一个结局显著、另一个不显著而宣称效应不同。
2. 历史队列插补后的癌症死亡 HR 为 1.286（1.020–1.620），而开发队列癌症结果不明确；须保留该差异，不能声称已排除癌症关联。
3. CVD 的四分位趋势检验 P=0.095；连续模型显著不代表四分位趋势显著。
4. 历史验证仍属同一 NHANES 调查体系；十年 AUC 增量区间跨 0，不能宣称已证明临床预测增益或超过 PREVENT。
5. SBP、DBP GWAS 使用逆正态转换残差，OR 不能按每 mmHg 解释。这两个暴露文件为多祖源汇总，欧洲 LD 参考是近似，且原血压 GWAS 调整了 BMI。
6. MR 针对部分成分及疾病结局，不证明复合 CMIBI 对 CVD 死亡的因果效应；存在异质性与多效性。脂质 MVMR 不等同于 ApoB 调整。
7. 转录组证据在人群间不一致，儿童糖尿病资料不能直接验证成人 CMIBI 机制。
8. 插补依赖 MAR 等假设，链轨迹和探索性收敛统计不能证明缺失机制正确。

## 9. 数据来源与可用性

原始输入来自 NHANES、GEO、GLGC、FinnGen、GWAS Catalog、GIGASTROKE 等项目。具体 accession、单位及方法来源以完整报告和附件元数据为准。

- [NHANES 官方站点](https://www.cdc.gov/nchs/nhanes/)
- [GSE98895](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98895)
- [GSE9006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9006)
- [GLGC 2013](https://csg.sph.umich.edu/willer/public/lipids2013/)
- [欧洲卒中 GCST90104539](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104539/)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)

欧洲卒中原始压缩文件的本次核对 MD5 为 `a967878d8efd6014444452137b07b7cf`。LD 参考使用 1000 Genomes EUR 的 `.bed/.bim/.fam` 文件。

原始数据未被本次分析覆盖。本地交付不等同于已公开发布；目前没有可声明的公开代码仓库 URL 或项目 DOI。上游源码与数据应遵守各自许可证和使用条款，不能将整个附件未经核查统一声明为某种开源许可证。

## 10. 投稿前需作者确认

作者顺序、正式单位、通讯邮箱、资助名称和编号、利益冲突及实际 CRediT 贡献仍需作者确认。新提供的 v2 稿件列出了候选期刊，但本项目尚未据此完成最终期刊格式适配或对外投稿。原稿 43 条文献的核查编号不应直接当作 v2 的引文编号。

README、完整报告及代码附件服务于文章重写和复现；作者应结合最终稿确认分析人群、变量定义、图表编号与引用顺序一致。
