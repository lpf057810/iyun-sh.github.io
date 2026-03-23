# 免疫浸润分析流水线

批量执行多种免疫浸润分析方法的Python脚本，支持并行执行和灵活的配置管理。

## 功能特点

- 支持多种免疫浸润分析方法
- 并行执行，默认8核
- 智能配置验证和错误处理
- 缓存机制，提高重复运行效率
- 详细的日志记录

## 支持的方法

| 编号 | 方法名称 | 说明 | 特殊要求 |
|------|---------|------|---------|
| 1 | cibersort | CIBERSORT免疫细胞反卷积 | 支持permutations参数 |
| 2 | epic | EPIC免疫浸润分析 | - |
| 3 | estimate | ESTIMATE免疫评分 | - |
| 4 | ips | IPS免疫表型评分 | - |
| 5 | mcpcounter | MCPcounter免疫细胞计数 | - |
| 6 | ssgsea | ssGSEA基因集富集分析 | **需要signature文件** |
| 7 | timer | TIMER免疫浸润 | 需要tissue参数 |
| 8 | xcell | xCell细胞富集评分 | - |

## 安装要求

### R环境

**R 版本**: >= 4.0.0

**R 包及版本要求**：

| 包名 | 版本 | 来源 | 说明 |
|------|------|------|------|
| IOBR | >= 1.0.0 | GitHub | 免疫浸润分析核心包 |
| optparse | >= 1.7.0 | CRAN | 命令行参数解析 |
| data.table | >= 1.14.0 | CRAN | 高效数据处理 |
| tidyverse | >= 1.3.0 | CRAN | 数据处理套件 |
| dplyr | >= 1.0.0 | CRAN | 数据操作 |
| qs2 | >= 1.0.0 | CRAN | 快速序列化 |
| ggplot2 | >= 3.4.0 | CRAN | 数据可视化 |
| rstatix | >= 0.7.0 | CRAN | 统计分析 |
| RColorBrewer | >= 1.1.0 | CRAN | 配色方案 |
| patchwork | >= 1.1.0 | CRAN | 图形组合 |
| WGCNA | >= 1.70 | CRAN | 相关性分析 |
| GSVA | >= 1.46.0 | Bioconductor | ssGSEA方法需要 |
| digest | >= 0.6.0 | CRAN | 哈希计算 |

### Python环境

| 包名 | 版本 | 说明 |
|------|------|------|
| Python | >= 3.6 | 运行环境 |
| configparser | >= 5.0 | 配置文件解析 |

### 安装命令

**CRAN 包**：
```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
install.packages(c("optparse", "data.table", "tidyverse", "dplyr",
                   "qs2", "ggplot2", "rstatix", "RColorBrewer",
                   "patchwork", "WGCNA", "digest"))
```

**Bioconductor 包**：
```r
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("IOBR", "GSVA"))
```

**IOBR 包**（需从 GitHub 安装）：
```r
devtools::install_github("IOBR/IOBR")
```

## 快速开始

### 1. 准备输入文件

**必需文件**:
- `expr_tpm.csv`: 表达矩阵，第一列为基因名(SYMBOL)，列为样本
- `group.csv`: 分组信息，包含sample和group列

**可选文件**:
- `gene_list.csv`: 基因列表，包含gene列（用于相关性分析）
- `gene_signatures.rds`: ssGSEA基因签名文件（RDS格式）

### 2. 配置文件

复制示例配置文件：
```bash
cp config.example.ini config.ini
```

编辑`config.ini`文件，配置输入输出参数：
```ini
[Methods]
methods = 1-8  # 执行所有方法

[Input]
expression = expr_tpm.csv
group = group.csv
genes =  # 可选
signature = gene_signatures.rds  # ssgsea必需

[Execution]
parallel_enabled = true
max_workers = 8
```

### 3. 运行分析

**使用运行脚本（推荐）**:
```bash
# 运行分析
/media/share/secure/immune_infiltration/run_immune_infiltration.sh \
    --config config.ini \
    -o /media/share/output/immune_result

# 删除输出目录（只有创建者可删除）
/media/share/secure/immune_infiltration/run_immune_infiltration.sh \
    --config config.ini \
    -o /media/share/output/immune_result \
    --rm
```

**直接使用 Python 脚本**:
```bash
python run_immune_infiltration.py --config config.ini -o output_dir
```

**仅运行特定方法**:
```bash
# 方法1: cibersort, 方法7: timer
python run_immune_infiltration.py --config config.ini -o output_dir --methods 1,7

# 使用范围: 执行方法1-5
python run_immune_infiltration.py --config config.ini -o output_dir --methods 1-5
```

**自定义并行数**:
```bash
python run_immune_infiltration.py --config config.ini -o output_dir --workers 4
```

**强制重新计算**:
```bash
python run_immune_infiltration.py --config config.ini -o output_dir --force
```

**预览将要执行的命令**:
```bash
python run_immune_infiltration.py --config config.ini -o output_dir --dry-run
```

**详细输出模式**:
```bash
python run_immune_infiltration.py --config config.ini -o output_dir --verbose
```

## 配置文件说明

### Methods节
```ini
[Methods]
# 方法编号: 1=cibersort, 2=epic, 3=estimate, 4=ips,
#           5=mcpcounter, 6=ssgsea, 7=timer, 8=xcell
# 支持范围(1-8)和逗号分隔(1,3,5)两种格式
methods = 1-8
```

### Input节
```ini
[Input]
# 表达矩阵文件(CSV格式,第一列为基因名SYMBOL) [必需]
expression = expr_tpm.csv

# 分组信息文件(CSV格式,包含sample和group列) [必需]
group = group.csv

# 基因列表文件(可选,用于相关性分析)
# CSV格式,包含gene列
# 支持的方法: cibersort(1), epic(2), mcpcounter(5), ssgsea(6), timer(7), xcell(8)
genes = 

# ssGSEA基因签名文件(可选,但运行ssgsea(6)时必需)
# RDS格式,包含基因集
signature = gene_signatures.rds
```

### Output节
```ini
[Output]
# 输出基础目录
base_dir = immune_infiltration_results

# 是否为每个方法创建子目录
create_subdirs = true
```

### Execution节
```ini
[Execution]
# 是否启用并行执行
parallel_enabled = true

# 最大并行任务数(默认8核)
max_workers = 8
```

### Cache节
```ini
[Cache]
# 是否启用缓存(提高重复运行速度)
enabled = true

# 缓存目录
cache_dir = cache

# 是否强制重新计算(忽略缓存)
force = false
```

### 方法特定参数
```ini
[Method1.CIBERSORT]
# CIBERSORT排列次数(默认1000)
perm = 1000

[Method7.TIMER]
# TIMER组织类型(如: LUAD, LUSC, GBM等)
tissue = LUAD
```

## 输出文件

每个方法会在输出目录下生成以下文件：
- `01.{method}_res.csv`: 免疫浸润结果
- `01.{method}.stacked_bar.pdf/png`: 堆叠柱状图
- `02.stat.{method}.csv`: 统计分析结果
- `03.DE.{method}.boxplot.pdf/png`: 箱线图
- `04.Correlation.gene.csv/pdf`: 基因相关性分析（如果提供了genes文件）

## 错误处理

### 配置验证失败
程序会在启动时验证配置文件，如果发现问题会立即报错：
- 表达矩阵或分组文件不存在
- ssgsea方法缺少signature文件
- 基因列表文件配置了但不存在

### 执行失败
如果某个R脚本执行失败：
- 错误信息会记录到日志文件
- 程序会继续执行其他方法
- 最后会显示所有失败的方法

### 常见问题

**Q: ssgsea方法报错需要signature文件？**
A: ssgsea方法需要RDS格式的基因签名文件，请在配置文件的[Input]节中配置`signature = gene_signatures.rds`，或从方法列表中移除ssgsea(6)。

**Q: 如何获取ssgSEA的基因签名文件？**
A: 可以使用IOBR包中的数据：
```r
library(IOBR)
cellmarkers <- IOBR::iobr_get_sig(type = "GeneSet")
saveRDS(cellmarkers, "gene_signatures.rds")
```

**Q: 某个方法执行失败，但我想继续执行其他方法？**
A: 程序默认会继续执行其他方法，最后会显示所有失败的方法，不会因为单个方法失败而停止。

**Q: 如何查看详细的执行日志？**
A: 使用`--verbose`参数可以查看更详细的输出信息：
```bash
python run_immune_infiltration.py --config config.ini --verbose
```

**Q: 修改了配置后如何强制重新计算？**
A: 使用`--force`参数可以强制重新计算，忽略缓存：
```bash
python run_immune_infiltration.py --config config.ini --force
```

## 命令行参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--config` | 配置文件路径(必需) | `--config config.ini` |
| `-o, --output` | 输出目录(覆盖配置文件) | `-o /path/to/output` |
| `--methods` | 覆盖要执行的方法 | `--methods 1,7` 或 `--methods 1-5` |
| `--workers` | 覆盖并行任务数 | `--workers 4` |
| `--force` | 强制重新计算 | `--force` |
| `--dry-run` | 仅显示命令不执行 | `--dry-run` |
| `--verbose` | 详细输出模式 | `--verbose` |
| `--rm` | 删除输出目录(仅创建者) | `--rm` |
| `--help` | 显示帮助信息 | `--help` |

## 输出目录权限

执行完成后，输出目录会自动设置为 `1777` 权限（所有人可读写），并创建 `.owner` 文件记录创建者：

```bash
# 查看权限
ls -la /media/share/output/immune_result/
# drwxrwxrwt  3 iyunlyl iyunlyl 4096 ... ./
# -r--r--r--  1 iyunlyl iyunlyl    7 ... .owner

# 删除输出目录（只有创建者可以删除）
/media/share/secure/immune_infiltration/run_immune_infiltration.sh \
    --config config.ini \
    -o /media/share/output/immune_result \
    --rm
```

## 日志文件

执行日志会保存到`run_immune_infiltration.log`文件，包含：
- 配置验证信息
- 每个方法的执行状态
- 错误和警告信息

## 示例输出

```
2024-01-29 16:00:00 - INFO - 免疫浸润分析流水线启动
2024-01-29 16:00:00 - INFO - 配置文件: config.ini
2024-01-29 16:00:00 - INFO - 验证配置...
2024-01-29 16:00:00 - INFO - ✓ 表达矩阵: expr_tpm.csv
2024-01-29 16:00:00 - INFO - ✓ 分组信息: group.csv
2024-01-29 16:00:00 - INFO - ✓ 基因签名: gene_signatures.rds

⚠️  提示: 基因相关性分析
  以下方法支持基因相关性分析: cibersort, epic, mcpcounter, timer, xcell
  如需分析基因与免疫细胞的相关性，请在配置文件[Input]中指定genes文件

2024-01-29 16:00:00 - INFO - 将要执行的方法: cibersort, epic, estimate, ips, mcpcounter, ssgsea, timer, xcell
2024-01-29 16:00:00 - INFO - 验证R脚本...
2024-01-29 16:00:00 - INFO - 开始执行任务...
2024-01-29 16:00:00 - INFO - 并行执行任务 (max_workers=8)...
2024-01-29 16:05:00 - INFO - ✓ cibersort 完成
2024-01-29 16:06:00 - INFO - ✓ epic 完成
...
2024-01-29 16:15:00 - INFO - 执行完成

============================================================
执行总结
============================================================

成功: 8/8
失败: 0/8
```

## 许可证

MIT License

---

## 代码保护与维护

### 使用方法（普通用户）

```bash
# 执行免疫浸润分析
sudo -u iyunlyl /media/share/secure/immune_infiltration/run_immune_infiltration.sh \
  --config config.ini \
  -o /media/share/output/immune_result

# 删除输出（只有创建者可删除）
sudo -u iyunlyl /media/share/secure/immune_infiltration/run_immune_infiltration.sh \
  --config config.ini \
  -o /media/share/output/immune_result \
  --rm
```

### 维护者开发流程

1. **在开发目录修改源码**
   ```bash
   vim ~/standardized_workflow/immune_infiltration/run_immune_infiltration.sh
   vim ~/standardized_workflow/immune_infiltration/scripts/*.R
   ```

2. **同步到生产目录**
   ```bash
   rsync -av --exclude='*.R' ~/standardized_workflow/immune_infiltration/ \
     /media/share/secure/immune_infiltration/
   ```

3. **设置权限**
   ```bash
   chmod 700 /media/share/secure/immune_infiltration/scripts
   chmod 400 /media/share/secure/immune_infiltration/scripts/*.R
   ```

### 详细保护方案

参见：[../../standardized_workflow/README_代码保护方案.md](../../standardized_workflow/README_代码保护方案.md)
