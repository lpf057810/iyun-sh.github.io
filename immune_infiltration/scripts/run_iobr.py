#!/usr/bin/env python3
"""
免疫浸润分析 - 直接用 Python + IOBR 生成结果和图表
"""
import os
import subprocess
import sys

def run_r_analysis():
    """运行 R 分析"""
    expr_file = "/media/desk16/share/secure/immune_infiltration/GSE126124.dat.csv"
    group_file = "/media/desk16/share/secure/immune_infiltration/GSE126124.group.csv"
    output_dir = "/media/desk16/share/secure/immune_infiltration/immune_infiltration_results"

    r_script = """
library(IOBR)
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
expr <- read.csv("{expr_file}", row.names=1)
group <- read.csv("{group_file}")

# 准备分组信息
group$sample <- as.character(group$sample)
group$group <- as.character(group$group)

# 创建输出目录
dir.create("{output_dir}/cibersort_output", recursive=TRUE)
dir.create("{output_dir}/epic_output", recursive=TRUE)
dir.create("{output_dir}/estimate_output", recursive=TRUE)

# CIBERSORT
cat("Running CIBERSORT...\\n")
ciber <- deconvo_tme(eset=expr, method="cibersort", perm=100)
names(ciber) <- gsub("_CIBERSORT$", "", names(ciber))

# 保存结果
write.csv(ciber, "{output_dir}/cibersort_output/01.ciber_res.csv", row.names=FALSE)
cat("Saved CIBERSORT results\\n")

# EPIC
cat("Running EPIC...\\n")
epic <- deconvo_tme(eset=expr, method="epic")
names(epic) <- gsub("_EPIC$", "", names(epic))
write.csv(epic, "{output_dir}/epic_output/01.epic_res.csv", row.names=FALSE)
cat("Saved EPIC results\\n")

# ESTIMATE
cat("Running ESTIMATE...\\n")
estimate <- deconvo_tme(eset=expr, method="estimate")
names(estimate) <- gsub("_ESTIMATE$", "", names(estimate))
write.csv(estimate, "{output_dir}/estimate_output/01.estimate_res.csv", row.names=FALSE)
cat("Saved ESTIMATE results\\n")

cat("Done!\\n")
""".format(expr_file=expr_file, group_file=group_file, output_dir=output_dir)

    # 保存临时脚本
    script_file = "/tmp/run_iobr.R"
    with open(script_file, "w") as f:
        f.write(r_script)

    # 运行 R 脚本
    result = subprocess.run(["Rscript", script_file],
                          capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    return result.returncode == 0

if __name__ == "__main__":
    print("=" * 50)
    print("运行免疫浸润分析 (IOBR)")
    print("=" * 50)
    run_r_analysis()
