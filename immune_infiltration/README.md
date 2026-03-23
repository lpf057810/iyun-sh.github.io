# 免疫浸润分析流程 (Immune Infiltration Analysis Pipeline)

## 概述

本流程通过集成多种算法进行肿瘤免疫微环境浸润水平评估，支持8种主流免疫浸润分析方法。分析完成后自动生成Word格式报告。

## 支持的方法

| 编号 | 方法 | 说明 |
|------|------|------|
| 1 | CIBERSORT | 基于线性支持向量机算法估算22种免疫细胞类型的相对比例 |
| 2 | EPIC | 利用混合线性模型估计多种免疫细胞群体的绝对丰度 |
| 3 | ESTIMATE | 计算 Stromal Score、Immune Score 和 ESTIMATE Score |
| 4 | IPS | 免疫表性评分 (Immunophenoscore) |
| 5 | MCPcounter | MCPcounter 免疫细胞计数 |
| 6 | ssGSEA | 基因集富集分析 (需要签名文件) |
| 7 | TIMER | TIMER 免疫浸润评估 (需要tissue参数) |
| 8 | xCell | xCell 细胞富集评分 |

## 快速开始

### 1. 配置文件

编辑 `config/config.ini`:

```ini
[Methods]
methods = 1-8

[Input]
expression = /path/to/expression.csv
group = /path/to/group.csv
signature = /path/to/signature.rds  # 仅ssGSEA需要

[Output]
base_dir = /path/to/result

[Cache]
enabled = true
force = false
# cache_dir 可留空：程序会自动使用 <base_dir>/qs2 存放 .qs2 缓存
```

### 2. 运行分析

```bash
cd /media/desk16/share/secure/immune_infiltration
./run_immune_infiltration.sh --config config/config.ini
```

### 3. 查看结果

- 结果保存在配置的 `base_dir` 目录（各方法子目录及 **`qs2/`** 下的 `.qs2` 缓存文件）
- 报告保存在 `report/` 目录
- 日志保存在 `logs/` 目录（含 `run_immune_infiltration.*.log` 与 **`sessionInfo.*.txt`**）

## 输入文件格式

### 表达矩阵 (CSV)

```
SYMBOL,GSM3591984,GSM3591985,...
Gene1,10.5,12.3,...
Gene2,8.2,9.1,...
```

### 分组文件 (CSV)

```
sample,group
GSM3591984,Control
GSM3591985,IBD
...
```

### 基因签名文件 (RDS)

**仅 ssGSEA 方法需要此文件**。用于定义免疫细胞类型的基因标志物。

RDS格式的列表对象:
```r
cellmarkers <- list(
  T_cell = c("CD3D", "CD3E", ...),
  B_cell = c("CD19", "CD20", ...),
  ...
)
saveRDS(cellmarkers, "signature.rds")
```

> **注意**：如果不使用 ssGSEA 方法，可以不配置此文件（或留空）。

## 高级用法

### 运行特定方法

```bash
# 单个方法
./run_immune_infiltration.sh --config config/config.ini --methods 1

# 多个方法（逗号分隔）
./run_immune_infiltration.sh --config config/config.ini --methods 1,3,5

# 连续范围
./run_immune_infiltration.sh --config config/config.ini --methods 1-4

# 跳着选择
./run_immune_infiltration.sh --config config/config.ini --methods 1,2,4,6,8
```

### 强制重新计算

```bash
./run_immune_infiltration.sh --config config/config.ini --force
```

### 查看帮助

```bash
./run_immune_infiltration.sh --help
```

## 输出说明

### 目录结构

```
immune_infiltration/
├── config/
│   └── config.ini                 # 配置文件
├── scripts/
│   ├── run_immune_infiltration.sh # 主入口脚本
│   ├── run_immune_infiltration.py # Python主程序
│   ├── generate_plots.py          # 图表生成
│   ├── generate_docx.py           # 报告生成(Python)
│   ├── generate_immune_report.R    # 报告生成(R)
│   ├── write_session_info.R        # 导出 R sessionInfo 到 logs/
│   ├── cibersort.R
│   ├── epic.R
│   ├── estimate.R
│   ├── ips.R
│   ├── mcpcounter.R
│   ├── ssgsea.R
│   ├── timer.R
│   └── xcell.R
├── test_data/
│   ├── GSE126124.dat.csv         # 测试数据（表达矩阵）
│   ├── GSE126124.group.csv       # 测试数据（分组）
│   └── immune_signature.rds       # 测试数据（签名文件）
├── result/                        # 分析结果（根据配置）
│   ├── qs2/                       # .qs2 中间缓存（自动创建，随 base_dir 变化）
│   ├── cibersort_output/
│   │   ├── 01.ciber_res.csv
│   │   ├── 02.stat.cibersort.csv
│   │   ├── 01.cibersort.stacked_bar.png/pdf
│   │   ├── 02.cibersort_boxplot.png
│   │   └── 03.DE.cibersort.png
│   ├── epic_output/
│   ├── estimate_output/
│   ├── ips_output/
│   ├── mcpcounter_output/
│   ├── ssgsea_output/
│   ├── timer_output/
│   └── xcell_output/
├── report/                        # 报告目录
│   └── immune_infiltration.日期时间.docx
└── logs/                         # 日志目录
    ├── run_immune_infiltration.日期时间.log   # Python 主流程日志
    └── sessionInfo.日期时间.txt             # R sessionInfo()（复现环境）
```

**本仓库说明**：`scripts/`（含 `scripts/cache/` 预置 `.qs2`）、`config/`（含 **`config.ini`** 与各示例 ini）、`docs/`、`shiny/`、`test_data/` 示例输入为完整内容；**`result/`、`report/`、`logs/` 在仓库内仅占位（`.gitkeep`）**，本地运行生成文件由 `.gitignore` 忽略。

### 各方法输出文件

- `{method}_output/01.{method}_res.csv` - 免疫浸润结果
- `{method}_output/02.stat.{method}.csv` - 统计结果（组间比较）
- `{method}_output/01.{method}_stacked_bar.png` - 堆叠图
- `{method}_output/02.{method}_boxplot.png` - 箱线图
- `{method}_output/03.DE.{method}.png` - 差异分析图

### 统计方法

- 组间比较使用 Mann-Whitney U 检验
- 多重检验校正使用 Benjamini-Hochberg 方法
- 显著性阈值: p < 0.05

## 依赖

### R 包

- optparse
- IOBR
- data.table
- tidyverse
- ggplot2
- rstatix
- GSVA (用于ssGSEA)
- officer (用于生成Word报告)
- flextable (用于生成Word表格)

### Python 包

- pandas
- numpy
- matplotlib
- seaborn
- scipy
- python-docx

## 故障排除

### 问题1: 数据质控失败

**症状**: 表达矩阵第一列名称应为 'SYMBOL'

**解决**: 确保表达矩阵第一列名为 SYMBOL，分组文件包含 sample 和 group 列

### 问题2: R脚本找不到函数库

**症状**: `cannot open the connection ... immune_functions.R`

**解决**: 确保Python调用R时工作目录正确，或使用绝对路径

### 问题3: ssGSEA失败

**症状**: `Error: long flag "arrays" is invalid`

**解决**: 确保配置文件指定了signature文件

### 问题4: 缓存问题

**解决**: 使用 `--force` 参数强制重新计算

## 测试数据

项目自带测试数据位于 `test_data/` 目录：

- `GSE126124.dat.csv` - 表达矩阵（78样本，18837基因）
- `GSE126124.group.csv` - 分组文件（21 Control + 57 IBD）

运行测试：
```bash
./run_immune_infiltration.sh --config config/config.ini
```

## 许可

本项目仅供研究使用。
