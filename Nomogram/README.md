# Nomogram 列线图分析工具

预测模型列线图（Nomogram）分析工具，支持 Logistic 回归模型构建、列线图绘制、校准曲线、DCA 和 ROC 分析，自动生成 SCI 标准的 Word 报告。

## 目录结构

```
Nomogram/
├── scripts/
│   ├── Nomogram.R                  # 主分析脚本
│   └── generate_nomogram_report.R  # 报告生成脚本
│   └── config/
│       └── config.ini              # 配置文件
├── test_data/                      # 测试数据
├── results/                        # 分析结果
├── logs/                          # 日志文件
└── report/                        # Word 报告
```

仓库中 `results/`、`report/`、`logs/` 三个目录会随代码一并提交（内含 `.gitkeep` 占位），本地运行产生的输出文件由 `.gitignore` 忽略，不进入版本库。

## 命令行使用

```bash
cd /media/desk16/share/secure/Nomogram

# 一键入口（推荐）
./run_nomogram.sh \
  -d <表达矩阵.csv> \
  -g <分组文件.csv> \
  -G <基因列表.csv> \
  -o results

# 或直接调用 R 脚本
Rscript scripts/Nomogram.R \
  -d <表达矩阵.csv> \
  -g <分组文件.csv> \
  -G <基因列表.csv> \
  -o results
```

### 参数说明

| 参数 | 简写 | 说明 | 必需 |
|------|------|------|------|
| `--data-file` | `-d` | 表达矩阵文件（CSV，第一列为基因名） | 是 |
| `--group-file` | `-g` | 分组文件（CSV，含 sample 和 group 列） | 是 |
| `--gene-file` | `-G` | 基因列表文件（CSV，含 gene 列） | 是 |
| `--comparison` | `-c` | 比较组，如 "Normal,Tumor" | 否（自动识别） |
| `--output-dir` | `-o` | 输出目录 | 否（默认 results/Nomogram） |
| `--prefix` | `-p` | 输出文件前缀 | 否 |
| `--calibrate-b` | - | 校准曲线 bootstrap 次数（默认1000） | 否 |
| `--quiet` | `-q` | 静默模式 | 否 |

### 使用示例

```bash
# 自动识别比较组（自动取第一组 vs 其余组）
Rscript scripts/Nomogram.R \
  -d test_data/GSE126124.dat.csv \
  -g test_data/GSE126124.group.csv \
  -G test_data/04.gene.csv

# 手动指定比较组
Rscript scripts/Nomogram.R \
  -d data.csv -g group.csv -G genes.csv \
  -c "Control,IBD" -o results -p MyStudy
```

## 输入文件格式

### 表达矩阵（CSV）

- **第一列**：基因名（SYMBOL）
- **其余列**：样本表达值

```
SYMBOL,Sample1,Sample2,Sample3
GeneA,10.5,12.3,8.7
GeneB,5.2,6.1,4.8
```

### 分组文件（CSV）

```
sample,group
Sample1,Control
Sample2,Control
Sample3,IBD
```

### 基因列表（CSV）

```
gene
GeneA
GeneB
GeneC
```

## 输出文件

| 文件 | 说明 |
|------|------|
| 01.nomogram.pdf/png | 列线图 |
| 02.calibrate.pdf/png | 校准曲线 |
| 03.dca.pdf/png | 决策曲线（DCA） |
| 04.roc.pdf/png | ROC 曲线 |
| 04.roc_resultss.csv | ROC 指标（AUC、敏感度、特异度、PPV、NPV、最佳阈值） |
| 05.model_coefficients.csv | 模型系数（包含 Coefficient、OR、95%CI、Expression Range、Points） |
| nomogram_YYYYMMDD_HHMMSS.docx | SCI 标准 Word 报告 |

## 报告内容

报告包含方法（Methods）和结果（Results）两大部分，结果部分所有文字根据实际数据动态生成：

**方法**：R软件版本及CRAN链接、rms/rmda/pROC包版本及参数

**结果**：
- 表1 模型系数（Coefficient、OR、95%CI、Expression Range、Points）及解读
- 图1 列线图及解读
- 图2 ROC曲线及解读（AUC、C-index、敏感度、特异度、PPV、NPV）
- 图3 校准曲线及解读（MAE、H-L检验p值）
- 图4 决策曲线（DCA）及解读

## 数据质控

脚本自动检查：

- 表达矩阵格式（基因为行 vs 基因为列，自动检测并转置）
- 样本量与基因数比例（建议样本/基因比 ≥ 3）
- 组间样本平衡
- 模型过拟合风险（C-index > 0.95 时警告）

## 配置文件

使用 `config/config.ini` 配置默认参数：

```ini
[input]
data_file = test_data/your_data.csv
group_file = test_data/your_group.csv
gene_file = test_data/your_genes.csv

[analysis]
calibrate_b = 1000
lrm_maxit = 1000

[directories]
result = results/Nomogram
```

配置好后直接运行：
```bash
Rscript scripts/Nomogram.R
```

## 常见问题

### Q: 表达矩阵格式错误？
A: 确保第一列为基因名（SYMBOL），其余列为样本。脚本会自动检测并提示。

### Q: 校准曲线无法收敛？
A: 样本量不足导致。尝试：1) 增加样本量；2) 减少基因数量；3) 降低 bootstrap 次数（`--calibrate-b 100`）

### Q: C-index 接近 1.0？
A: 可能过拟合。检查：样本/基因比是否 ≥ 3？两组样本是否平衡？

### Q: 报告字体？
A: 报告中中文使用 SimSun（宋体），英文和数字使用 Times New Roman。

## 环境要求

- R >= 4.0
- 依赖包：data.table, dplyr, ggplot2, rms, ResourceSelection, rmda, pROC, survival, optparse, officer, flextable, magrittr

## 安装依赖

```R
install.packages(c("data.table", "dplyr", "ggplot2", "rms", "ResourceSelection",
                   "rmda", "pROC", "survival", "optparse", "officer", "flextable", "magrittr"))
```

## 许可证

仅供研究使用。
