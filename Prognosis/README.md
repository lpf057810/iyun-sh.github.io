# 预后分析 (Prognosis)

**项目目录**: `/media/desk16/share/secure/Prognosis`

> **重要**: 每次进行预后分析前，请确认项目目录和输入文件路径。

## 功能

- **Cox回归分析** - 单因素Cox + PH假定检验
- **Lasso特征选择** - 自动筛选关键预后基因
- **风险模型构建** - 计算风险评分，分高低风险组
- **KM生存分析** - 生存曲线可视化
- **时间ROC分析** - 1/3/5年预测效能评估
- **热图绘制** - 预后基因表达可视化
- **Word报告** - 自动生成分析报告

## 使用流程

1. **确认项目目录**: `/media/desk16/share/secure/Prognosis`
2. **准备输入文件**:
   - 表达矩阵 (CSV): 第一列为基因名(SYMBOL)
   - 分组文件 (CSV): 包含 sample, group 列
   - 生存数据 (CSV): 包含 sample, OS.time, OS 列
   - 基因列表 (CSV): 包含 gene 列
3. **运行分析**: 使用 `./run_prognosis.sh`
4. **查看结果**: 在 output 目录查看图表和报告

## 输入文件格式

### 表达矩阵
```
SYMBOL,Sample1,Sample2,Sample3,...
GeneA,5.2,6.1,4.8,...
GeneB,3.4,3.9,4.2,...
```

### 分组文件
```
sample,group
Sample1,Tumor
Sample2,Normal
```

### 生存数据
```
sample,OS.time,OS
Sample1,365,1
Sample2,730,0
```
- OS: 生存状态 (1=死亡, 0=存活)
- OS.time: 生存时间 (天)

### 基因列表
```
gene
GeneA
GeneB
GeneC
```

## 使用方法

### 一键运行 (推荐)

```bash
cd /media/desk16/share/secure/Prognosis
./run_prognosis.sh \
  -w . \
  -e test_data/expression.csv \
  -g test_data/group.csv \
  -s test_data/survival.csv \
  -G test_data/genes.csv
```

### 分步运行

```bash
cd /media/desk16/share/secure/Prognosis

# 步骤1: Cox + Lasso
Rscript scripts/Cox.R \
  --expr test_data/expression.csv \
  --group test_data/group.csv \
  --survival test_data/survival.csv \
  --genes test_data/genes.csv \
  --output result/Cox

# 步骤2: 风险模型 + Word报告
Rscript scripts/RiskModel.R \
  --expr test_data/expression.csv \
  --group test_data/group.csv \
  --survival test_data/survival.csv \
  --lasso-coef result/Cox/03.Lasso_Coefficients.csv \
  --output result/RiskModel \
  --report
```

### 单独生成报告

```bash
cd /media/desk16/share/secure/Prognosis
Rscript scripts/generate_prognosis_report.R
```

## 输出文件

### result/Cox/

| 文件 | 说明 |
|------|------|
| `01.univariate_cox_all.csv` | 所有基因Cox结果 |
| `01.univariate_cox_ph.csv` | PH假定通过的基因 |
| `02.univariate_cox_0.05_forest.pdf/png` | 森林图 |
| `03.lasso.CV.pdf/png` | Lasso交叉验证图 |
| `03.lasso.Coef.pdf/png` | Lasso系数路径图 |
| `03.Lasso_Coefficients.csv` | Lasso筛选的基因系数 |
| `03.lasso_genes.csv` | 筛选的基因列表 |

### result/RiskModel/

| 文件 | 说明 |
|------|------|
| `01.rs_dat.csv` | 风险评分数据 |
| `01.KM.pdf/png` | KM生存曲线 |
| `02.timeROC.pdf/png` | 时间ROC曲线 |
| `03.SurvStat.pdf/png` | 生存状态分布图 |
| `04.RiskScore.pdf/png` | 风险评分分布图 |
| `05.HeatMap.pdf/png` | 预后基因热图 |

### report/

| 文件 | 说明 |
|------|------|
| `Prognosis_Report_*.docx` | Word分析报告 |

## 与列线图板块衔接

Prognosis分析完成后，可使用Nomogram板块进一步构建列线图：

```bash
cd /media/desk16/share/secure/Nomogram
./run_nomogram.sh -d expr.csv -g group.csv -G genes.csv
```

Nomogram板块会基于同样的基因构建列线图、校准曲线、DCA等可视化。

## 测试数据

测试数据位于 `test_data/` 目录:

```bash
cd /media/desk16/share/secure/Prognosis
./run_prognosis.sh -w . -e test_data/expression.csv -g test_data/group.csv -s test_data/survival.csv -G test_data/genes.csv
```

## 注意事项

- 分组文件中需要包含 `Tumor` 列才会筛选肿瘤样本
- Lasso系数文件用于第二阶段风险模型分析
- 测试数据基因数较少，结果可能不显著
- 报告一键生成在 RiskModel 步骤添加 `--report` 参数

## 项目结构

```
Prognosis/
├── scripts/
│   ├── Cox.R                    # Cox回归分析
│   ├── RiskModel.R              # 风险模型
│   ├── generate_prognosis_report.R  # 报告生成
│   ├── ikl_function.R          # 辅助函数
│   ├── Nomogram.R               # 列线图(独立板块)
│   └── run_prognosis.sh         # 入口脚本
├── test_data/                    # 测试数据
│   ├── expression.csv
│   ├── group.csv
│   ├── survival.csv
│   └── genes.csv
├── result/                      # 分析结果
│   ├── Cox/
│   └── RiskModel/
├── report/                      # 报告输出
└── README.md
```
