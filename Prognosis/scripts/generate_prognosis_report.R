#!/usr/bin/env Rscript
# Prognosis Analysis Report Generator

library(optparse)
library(officer)
library(flextable)
library(ggplot2)
library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(timeROC)
library(qs)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
package_version_safe <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("unknown")
  as.character(utils::packageVersion(pkg))
}

# .qs2 保存（按分析节点，与 MachineLearn 一致）
save_qs2 <- function(object, path) {
  qs::qsave(object, path)
  invisible(path)
}

# Helpers - Chinese: SimSun, English/Numbers: Times New Roman
cn_font <- "SimSun"
en_font <- "Times New Roman"

build_fpar <- function(text, font_size = 10.5, bold = FALSE, align = "left") {
  if (!nzchar(trimws(text))) {
    return(officer::fpar("", fp_p = officer::fp_par(text.align = align)))
  }
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  is_cn <- grepl("[一-龥]", chars)
  runs <- list()
  start <- 1L
  for (i in seq_along(chars)) {
    is_break <- i == length(chars) || is_cn[[i + 1L]] != is_cn[[start]]
    if (is_break) {
      segment <- paste(chars[start:i], collapse = "")
      font_family <- if (is_cn[[start]]) cn_font else en_font
      runs[[length(runs) + 1L]] <- officer::ftext(
        segment,
        prop = officer::fp_text(font.family = font_family, font.size = font_size, bold = bold)
      )
      start <- i + 1L
    }
  }
  para_props <- officer::fp_par(text.align = align)
  do.call(officer::fpar, c(runs, list(fp_p = para_props)))
}

mk_par <- function(doc, text, font_size = 10.5, bold = FALSE, indent = FALSE) {
  if (indent) text <- paste0("\u3000\u3000", text)
  body_add_fpar(doc, value = build_fpar(text, font_size = font_size, bold = bold, align = "left"))
}
mk_par_center <- function(doc, text, font_size = 10.5, bold = FALSE) {
  body_add_fpar(doc, value = build_fpar(text, font_size = font_size, bold = bold, align = "center"))
}
mk_blank <- function(doc) { body_add_par(doc, "") }
mk_img <- function(doc, path, w, h) {
  img_block <- fpar(
    external_img(src = path, width = w, height = h),
    fp_p = fp_par(text.align = "center")
  )
  body_add_fpar(doc, value = img_block)
}

# Parse args
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL),
  make_option(c("-o", "--output"), type = "character", default = NULL),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL),
  make_option(c("--cox-dir"), type = "character", default = NULL),
  make_option(c("--risk-dir"), type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$timestamp)) opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_root <- "/media/desk16/share/secure/Prognosis/report"
if (!dir.exists(report_root)) dir.create(report_root, recursive = TRUE)
template_path <- "/media/desk16/share/secure/Prognosis/templates/prognosis_report_template.docx"

# 报告输出强制落在 report/ 目录，不受结果目录参数影响
if (is.null(opt$output)) {
  opt$output <- file.path(report_root, paste0("Prognosis_Report_", opt$timestamp, ".docx"))
} else {
  out_name <- if (grepl("\\.docx$", basename(opt$output), ignore.case = TRUE)) {
    basename(opt$output)
  } else {
    paste0("Prognosis_Report_", opt$timestamp, ".docx")
  }
  opt$output <- file.path(report_root, out_name)
}
if (is.null(opt$input)) opt$input <- "/media/desk16/share/secure/Prognosis/result"

cat("================================================================================\n")
cat("  Prognosis Report Generator\n")
cat("  Timestamp:", opt$timestamp, "\n")
cat("================================================================================\n\n")

# Set paths
if (!is.null(opt$`cox-dir`) && !is.null(opt$`risk-dir`)) {
  cox_dir <- opt$`cox-dir`
  risk_dir <- opt$`risk-dir`
  # 即使传入目录或任意路径，也统一写到 report_root
  out_name <- if (grepl("\\.docx$", basename(opt$output), ignore.case = TRUE)) {
    basename(opt$output)
  } else {
    paste0("Prognosis_Report_", opt$timestamp, ".docx")
  }
  opt$output <- file.path(report_root, out_name)
} else {
  result_dir <- opt$input
  cox_dir <- file.path(result_dir, "Cox")
  risk_dir <- file.path(result_dir, "RiskModel")
}

# Read data
cox_all <- NULL; cox_ph <- NULL; lasso_genes <- NULL; lasso_coef <- NULL
if (dir.exists(cox_dir)) {
  if (file.exists(f1 <- file.path(cox_dir, "01.univariate_cox_all.csv"))) {
    cox_all <- read.csv(f1, stringsAsFactors = FALSE); cat("Loaded Cox all\n") }
  if (file.exists(f2 <- file.path(cox_dir, "01.univariate_cox_ph.csv"))) {
    cox_ph <- read.csv(f2, stringsAsFactors = FALSE); cat("Loaded Cox PH\n") }
  if (file.exists(f3 <- file.path(cox_dir, "03.lasso_genes.csv"))) {
    lasso_genes <- read.csv(f3, stringsAsFactors = FALSE); cat("Loaded Lasso genes\n") }
  if (file.exists(f4 <- file.path(cox_dir, "03.Lasso_Coefficients.csv"))) {
    lasso_coef <- read.csv(f4, stringsAsFactors = FALSE); cat("Loaded Lasso coef\n") }
}
risk_dat <- NULL
if (dir.exists(risk_dir) && file.exists(f <- file.path(risk_dir, "01.rs_dat.csv"))) {
  risk_dat <- read.csv(f, stringsAsFactors = FALSE); cat("Loaded RiskModel\n")
}
heatmap_gene_stats <- NULL
if (dir.exists(risk_dir) && file.exists(hf <- file.path(risk_dir, "05.heatmap_gene_stats.csv"))) {
  heatmap_gene_stats <- read.csv(hf, stringsAsFactors = FALSE); cat("Loaded HeatMap gene stats\n")
}
km_stats <- NULL
if (dir.exists(risk_dir) && file.exists(kf <- file.path(risk_dir, "01.KM_stats.csv"))) {
  km_stats <- read.csv(kf, stringsAsFactors = FALSE); cat("Loaded KM stats\n")
}
survstat_stats <- NULL
if (dir.exists(risk_dir) && file.exists(sf <- file.path(risk_dir, "03.SurvStat_stats.csv"))) {
  survstat_stats <- read.csv(sf, stringsAsFactors = FALSE); cat("Loaded SurvStat stats\n")
}
riskscore_stats <- NULL
if (dir.exists(risk_dir) && file.exists(rf <- file.path(risk_dir, "04.RiskScore_stats.csv"))) {
  riskscore_stats <- read.csv(rf, stringsAsFactors = FALSE); cat("Loaded RiskScore stats\n")
}

# 优先读取 RiskModel 输出的固定AUC文件，避免报告端二次计算偏差
auc_from_file <- c(`1` = NA_real_, `3` = NA_real_, `5` = NA_real_)
auc_file <- file.path(risk_dir, "02.timeROC_auc.csv")
if (file.exists(auc_file)) {
  auc_df <- tryCatch(read.csv(auc_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(auc_df) && nrow(auc_df) > 0 && all(c("time_year", "auc") %in% colnames(auc_df))) {
    for (yr in c(1, 3, 5)) {
      hit <- auc_df[as.integer(auc_df$time_year) == yr, "auc", drop = TRUE]
      if (length(hit) > 0) auc_from_file[as.character(yr)] <- suppressWarnings(as.numeric(hit[1]))
    }
    cat("Loaded timeROC AUCs from 02.timeROC_auc.csv\n")
  }
}

# Compute stats
km_pval <- NA; auc_1y <- NA; auc_3y <- NA; auc_5y <- NA
total_n <- NA_integer_; high_n <- NA_integer_; low_n <- NA_integer_
high_death <- NA_integer_; low_death <- NA_integer_
risk_cutoff <- NA_real_; high_death_rate <- NA_real_; low_death_rate <- NA_real_
if (!is.null(risk_dat) && nrow(risk_dat) > 0 && "OS" %in% colnames(risk_dat)) {
  total_n <- nrow(risk_dat)
  if ("RiskGroup" %in% colnames(risk_dat)) {
    high_n <- sum(risk_dat$RiskGroup == "High", na.rm = TRUE)
    low_n  <- sum(risk_dat$RiskGroup == "Low",  na.rm = TRUE)
    high_death <- sum(risk_dat$RiskGroup == "High" & risk_dat$OS == 1, na.rm = TRUE)
    low_death  <- sum(risk_dat$RiskGroup == "Low"  & risk_dat$OS == 1, na.rm = TRUE)
    high_death_rate <- if (!is.na(high_n) && high_n > 0) high_death / high_n else NA_real_
    low_death_rate  <- if (!is.na(low_n) && low_n > 0) low_death / low_n else NA_real_
  }
  if ("RiskScore" %in% colnames(risk_dat)) {
    risk_cutoff <- suppressWarnings(as.numeric(stats::median(risk_dat$RiskScore, na.rm = TRUE)))
    if (!is.finite(risk_cutoff)) risk_cutoff <- NA_real_
  }
  fit <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = risk_dat)
  km_pval <- pchisq(fit$chisq, df = 1, lower.tail = FALSE)
  if ("RiskScore" %in% colnames(risk_dat) && any(is.na(auc_from_file))) {
    roc <- tryCatch(
      timeROC(
        T = risk_dat$OS.time,
        delta = risk_dat$OS,
        marker = risk_dat$RiskScore,
        times = c(1, 3, 5) * 365.25,
        cause = 1,
        weighting = "marginal"
      ),
      error = function(e) NULL
    )
    if (!is.null(roc) && !is.null(roc$AUC) && length(roc$AUC) >= 3) {
      auc_1y <- round(roc$AUC[1], 3)
      auc_3y <- round(roc$AUC[2], 3)
      auc_5y <- round(roc$AUC[3], 3)
    }
  }
}
if (!is.na(auc_from_file["1"])) auc_1y <- round(auc_from_file["1"], 3)
if (!is.na(auc_from_file["3"])) auc_3y <- round(auc_from_file["3"], 3)
if (!is.na(auc_from_file["5"])) auc_5y <- round(auc_from_file["5"], 3)

# 保存 .qs2（按分析节点，与 MachineLearn 一致）
result_dir_for_qs <- if (!is.null(opt$`cox-dir`) && !is.null(opt$`risk-dir`)) {
  cox_parent <- dirname(normalizePath(opt$`cox-dir`, mustWork = FALSE))
  risk_parent <- dirname(normalizePath(opt$`risk-dir`, mustWork = FALSE))
  if (identical(cox_parent, risk_parent)) cox_parent else cox_parent
} else {
  opt$input
}
if (!dir.exists(result_dir_for_qs)) dir.create(result_dir_for_qs, recursive = TRUE)
# dataset_name 取自输出目录 basename，命名可由 -o 路径控制
dataset_name <- basename(result_dir_for_qs)
cox_node_path       <- file.path(result_dir_for_qs, sprintf("%s_cox_node.qs2",     dataset_name))
riskmodel_node_path <- file.path(result_dir_for_qs, sprintf("%s_riskmodel_node.qs2", dataset_name))
save_qs2(list(
  cox_all = cox_all, cox_ph = cox_ph,
  lasso_genes = lasso_genes, lasso_coef = lasso_coef
), cox_node_path); cat("Saved:", cox_node_path, "\n")
save_qs2(list(
  risk_dat = risk_dat,
  km_pval = km_pval, auc_1y = auc_1y, auc_3y = auc_3y, auc_5y = auc_5y
), riskmodel_node_path); cat("Saved:", riskmodel_node_path, "\n")

# Build document
doc <- if (file.exists(template_path)) {
  cat("Using report template:", template_path, "\n")
  read_docx(path = template_path)
} else {
  read_docx()
}
doc <- mk_blank(doc)
doc <- mk_par(doc, "预后分析报告", font_size = 16, bold = TRUE)
doc <- mk_blank(doc)

# ---- 1. Methods ----
doc <- mk_par(doc, "1. 方法", font_size = 14, bold = TRUE)
doc <- mk_blank(doc)

doc <- mk_par(doc, "1.1 方法概述", font_size = 12, bold = TRUE)
doc <- mk_blank(doc)
doc <- mk_par(doc, paste0(
  "采用R(v", as.character(getRversion()), ", https://www.r-project.org)预后建模；",
  "survival(v", package_version_safe("survival"), ")、glmnet(v", package_version_safe("glmnet"), ")、",
  "survminer(v", package_version_safe("survminer"), ")、timeROC(v", package_version_safe("timeROC"), ")",
  "（https://cran.r-project.org）用于Cox+PH检验(p>0.05)、10折Lasso(family='cox',lambda.min)、",
  "KM与1/3/5年AUC评估。"
), indent = TRUE)
doc <- mk_blank(doc)

# ---- 2. Results ----
doc <- mk_par(doc, "2. 结果", font_size = 14, bold = TRUE)
doc <- mk_blank(doc)

doc <- mk_par(doc, "2.1 汇总统计", font_size = 12, bold = TRUE)
doc <- mk_blank(doc)
if (!is.null(risk_dat) && nrow(risk_dat) > 0) {
  summary_parts <- c(
    paste0("总样本", total_n, "例，其中高风险组", high_n, "例、低风险组", low_n, "例；"),
    paste0("高风险组死亡", high_death, "例（", sprintf("%.1f", high_death / high_n * 100), "%），低风险组死亡", low_death, "例（", sprintf("%.1f", low_death / low_n * 100), "%）。")
  )
  if (!is.na(km_pval)) summary_parts <- c(summary_parts, sprintf("Log-rank检验p=%.4g。", km_pval))
  if (!is.na(auc_1y))  summary_parts <- c(summary_parts, sprintf("1/3/5年AUC分别为%.3f/%.3f/%.3f。", auc_1y, auc_3y, auc_5y))
  summary_parts <- c(summary_parts, paste0("高风险组占比", sprintf("%.1f", high_n / total_n * 100), "%，低风险组占比", sprintf("%.1f", low_n / total_n * 100), "%。"))
  doc <- mk_par(doc, paste0(summary_parts, collapse = ""), indent = TRUE)
  doc <- mk_par(doc, "相关结果表：`01.rs_dat.csv`。", indent = TRUE)
} else {
  doc <- mk_par(doc, "无风险模型结果可用。")
}
doc <- mk_blank(doc)

# 2.2 Lasso
doc <- mk_par(doc, "2.2 Lasso特征选择结果", font_size = 12, bold = TRUE)
doc <- mk_blank(doc)
if (!is.null(lasso_coef) && nrow(lasso_coef) > 0) {
  lasso_gene_names <- if ("gene" %in% colnames(lasso_coef)) paste(lasso_coef$gene, collapse = "、") else paste(lasso_coef[[1]], collapse = "、")
  doc <- mk_par(doc, paste0(
    "Lasso回归（family='cox'，10折交叉验证，lambda.min）筛选出", nrow(lasso_coef), "个关键预后基因（", lasso_gene_names, "），表1为各基因的回归系数。"
  ), indent = TRUE)
  td <- lasso_coef
  if ("gene" %in% colnames(td) && "coefficient" %in% colnames(td)) {
    td <- data.frame(Gene = td$gene, Coefficient = round(as.numeric(td$coefficient), 6))
  } else {
    colnames(td)[1:2] <- c("Gene", "Coefficient")
    td[,2] <- round(as.numeric(td[,2]), 6)
  }
  ft <- flextable(td) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all") %>%
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "all") %>%
    align(align = "left", j = 1, part = "body") %>%
    set_table_properties(layout = "autofit", width = 1)
  doc <- mk_par_center(doc, "表1：Lasso筛选基因回归系数表")
  doc <- doc %>% body_add_flextable(ft, align = "center")
  doc <- mk_par(doc, "相关结果表：`03.lasso_genes.csv`，`03.Lasso_Coefficients.csv`。", indent = TRUE)
  doc <- mk_blank(doc)
} else {
  doc <- mk_par(doc, "无Lasso筛选结果。")
  doc <- mk_blank(doc)
}

# 2.3 Cox
doc <- mk_par(doc, "2.3 单因素Cox分析结果", font_size = 12, bold = TRUE)
doc <- mk_blank(doc)
if (!is.null(cox_all) && nrow(cox_all) > 0) {
  cox_desc <- paste0("共对", nrow(cox_all), "个基因进行单因素Cox分析")
  if (!is.null(cox_ph) && nrow(cox_ph) > 0) {
    cox_desc <- paste0(cox_desc, "，其中", nrow(cox_ph), "个通过PH假定检验（cox.zph p>0.05）进入后续筛选")
  }
  cox_desc <- paste0(cox_desc, "（见表2）。")
  doc <- mk_par(doc, cox_desc, indent = TRUE)
  gene_col <- if ("Gene" %in% colnames(cox_all)) "Gene" else
               if ("gene" %in% colnames(cox_all)) "gene" else colnames(cox_all)[1]
  has_sep <- "HR" %in% colnames(cox_all) && "L95CI" %in% colnames(cox_all)
  if (has_sep) {
    pval_col <- if ("p.value" %in% colnames(cox_all)) "p.value" else grep("p", colnames(cox_all), value = TRUE)[1]
    td <- head(cox_all[, c(gene_col, "HR", "L95CI", "H95CI", pval_col)], 20)
    colnames(td) <- c("Gene", "HR", "95%CI_Lower", "95%CI_Upper", "P-value")
  } else {
    td <- head(cox_all[, 1:min(5, ncol(cox_all))], 20)
  }
  ft <- flextable(td) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    fontsize(size = 9, part = "all") %>%
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "all") %>%
    align(align = "left", j = 1, part = "body") %>%
    set_table_properties(layout = "autofit", width = 1)
  doc <- mk_par_center(doc, "表2：单因素Cox分析结果表")
  doc <- doc %>% body_add_flextable(ft, align = "center")
  doc <- mk_par(doc, "相关结果表：`01.univariate_cox_all.csv`，`01.univariate_cox_ph.csv`。", indent = TRUE)
  doc <- mk_blank(doc)
} else {
  doc <- mk_par(doc, "无Cox分析结果。")
  doc <- mk_blank(doc)
}

# 2.4 Figures
doc <- mk_par(doc, "2.4 可视化结果", font_size = 12, bold = TRUE)
doc <- mk_blank(doc)

# HeatMap 动态解释：优先使用图5专属结果表
heatmap_detail_text <- ""
if (!is.null(heatmap_gene_stats) && nrow(heatmap_gene_stats) > 0 &&
    all(c("gene", "delta_high_minus_low") %in% colnames(heatmap_gene_stats))) {
  deltas <- suppressWarnings(as.numeric(heatmap_gene_stats$delta_high_minus_low))
  valid <- is.finite(deltas)
  if (any(valid)) {
    deltas <- deltas[valid]
    genes <- as.character(heatmap_gene_stats$gene[valid])
    up_n <- sum(deltas > 0, na.rm = TRUE)
    down_n <- sum(deltas < 0, na.rm = TRUE)
    top_i <- which.max(abs(deltas))
    top_gene <- genes[top_i]
    top_diff <- deltas[top_i]
    top_dir <- ifelse(top_diff > 0, "高风险组更高", "低风险组更高")
    heatmap_detail_text <- sprintf(
      "纳入热图的%d个基因中，高风险组相对低风险组表达升高%d个、降低%d个；差异最大的基因为%s（High-Low均值差=%.4f，%s）。",
      length(deltas), up_n, down_n, top_gene, top_diff, top_dir
    )
  }
}

fig_count <- 0
plots <- list(
  c("01.KM", "Kaplan-Meier生存曲线",
    ifelse(is.na(km_pval), "",
           sprintf("X轴为随访时间（time），Y轴为累计生存概率。红/蓝曲线分别表示高风险组与低风险组，曲线短线标记为删失样本。组间差异采用Log-rank检验（p = %.4f）。", km_pval))),
  c("02.timeROC", "时间依赖ROC曲线",
    ifelse(is.na(auc_1y), "",
           sprintf("X轴为1-特异度（False Positive Rate），Y轴为敏感度（True Positive Rate）。分别给出1年、3年、5年预测曲线；对应AUC为%.3f、%.3f、%.3f。", auc_1y, auc_3y, auc_5y))),
  c("03.SurvStat", "生存状态分布图",
    "X轴为按风险评分升序排列的患者，Y轴为生存时间（OS.time）。点/三角分别表示存活/死亡状态；颜色区分高低风险组。灰色虚线为风险分组截断位置。"),
  c("04.RiskScore", "风险评分分布图",
    "X轴为按风险评分升序排列的患者，Y轴为风险评分（RiskScore）。红色虚线表示中位风险评分，并据此划分High/Low风险组。"),
  c("05.HeatMap", "预后基因表达热图",
    "列为患者样本（按风险评分排序），行为Lasso筛选预后基因。颜色表示标准化表达量Z-score（红高蓝低，截断至[-2,2]）；行方向进行层次聚类。")
)

for (p in plots) {
  pf <- file.path(risk_dir, paste0(p[1], ".png"))
  if (!file.exists(pf)) pf <- file.path(risk_dir, paste0(p[1], ".pdf"))
  if (file.exists(pf)) {
    fig_count <- fig_count + 1
    doc <- mk_img(doc, pf, 5, 4)
    doc <- mk_par_center(doc, paste0("图", fig_count, "：", p[2]))
    if (nzchar(p[3])) doc <- mk_par_center(doc, paste0("图注：", p[3]))
    # 参考列线图模块：每张图后追加可复用的客观数值解读
    figure_text <- switch(p[1],
      "01.KM" = if (!is.na(km_pval) && !is.na(high_death_rate) && !is.na(low_death_rate)) {
        sprintf("图%d结果显示，Log-rank检验p=%.4f；高风险组死亡率为%.2f%%，低风险组死亡率为%.2f%%。", fig_count, km_pval, high_death_rate * 100, low_death_rate * 100)
      } else if (!is.null(km_stats) && nrow(km_stats) > 0 &&
                 all(c("logrank_p", "high_death", "high_n", "low_death", "low_n") %in% colnames(km_stats))) {
        hrate <- suppressWarnings(as.numeric(km_stats$high_death[1]) / as.numeric(km_stats$high_n[1]) * 100)
        lrate <- suppressWarnings(as.numeric(km_stats$low_death[1]) / as.numeric(km_stats$low_n[1]) * 100)
        sprintf("图%d结果显示，Log-rank检验p=%.4f；高风险组死亡率为%.2f%%，低风险组死亡率为%.2f%%。", fig_count, as.numeric(km_stats$logrank_p[1]), hrate, lrate)
      } else {
        sprintf("图%d结果显示，高低风险组生存曲线存在可视化分离趋势。", fig_count)
      },
      "02.timeROC" = if (!is.na(auc_1y) && !is.na(auc_3y) && !is.na(auc_5y)) {
        sprintf("图%d结果显示，模型在1/3/5年的AUC分别为%.3f/%.3f/%.3f。", fig_count, auc_1y, auc_3y, auc_5y)
      } else {
        sprintf("图%d结果显示，时间依赖ROC用于评估模型在不同时间点的区分能力。", fig_count)
      },
      "03.SurvStat" = if (!is.na(total_n) && !is.na(high_n) && !is.na(low_n)) {
        sprintf("图%d结果显示，样本按风险评分排序后呈现事件分层分布（总样本%d例；High/Low=%d/%d）。", fig_count, total_n, high_n, low_n)
      } else if (!is.null(survstat_stats) && nrow(survstat_stats) > 0 &&
                 all(c("total_samples", "high_n", "low_n") %in% colnames(survstat_stats))) {
        sprintf("图%d结果显示，样本按风险评分排序后呈现事件分层分布（总样本%d例；High/Low=%d/%d）。", fig_count, as.integer(survstat_stats$total_samples[1]), as.integer(survstat_stats$high_n[1]), as.integer(survstat_stats$low_n[1]))
      } else {
        sprintf("图%d结果显示，生存状态随风险评分排序呈现分布差异。", fig_count)
      },
      "04.RiskScore" = if (!is.na(risk_cutoff)) {
        sprintf("图%d结果显示，风险评分中位数阈值为%.4f，并据此完成High/Low分组。", fig_count, risk_cutoff)
      } else if (!is.null(riskscore_stats) && nrow(riskscore_stats) > 0 && "cutoff" %in% colnames(riskscore_stats)) {
        sprintf("图%d结果显示，风险评分中位数阈值为%.4f，并据此完成High/Low分组。", fig_count, as.numeric(riskscore_stats$cutoff[1]))
      } else {
        sprintf("图%d结果显示，风险评分分布用于定义高低风险分组。", fig_count)
      },
      "05.HeatMap" = if (nzchar(heatmap_detail_text)) {
        sprintf("图%d结果显示，%s", fig_count, heatmap_detail_text)
      } else if (!is.null(lasso_coef) && nrow(lasso_coef) > 0) {
        sprintf("图%d结果显示，热图展示了%d个Lasso筛选基因在样本中的表达异质性。", fig_count, nrow(lasso_coef))
      } else {
        sprintf("图%d结果显示，热图用于观察预后相关基因表达模式。", fig_count)
      },
      ""
    )
    related_table_text <- switch(p[1],
      "01.KM" = if (file.exists(file.path(risk_dir, "01.KM_stats.csv"))) "相关结果表：`01.KM_stats.csv`。" else "",
      "02.timeROC" = if (file.exists(file.path(risk_dir, "02.timeROC_auc.csv"))) "相关结果表：`02.timeROC_auc.csv`。" else "",
      "03.SurvStat" = if (file.exists(file.path(risk_dir, "03.SurvStat_stats.csv"))) "相关结果表：`03.SurvStat_stats.csv`。" else "",
      "04.RiskScore" = if (file.exists(file.path(risk_dir, "04.RiskScore_stats.csv"))) "相关结果表：`04.RiskScore_stats.csv`。" else "",
      "05.HeatMap" = if (file.exists(file.path(risk_dir, "05.heatmap_gene_stats.csv"))) "相关结果表：`05.heatmap_gene_stats.csv`。" else "",
      ""
    )
    combined_fig_text <- paste0(figure_text, if (nzchar(related_table_text)) related_table_text else "")
    if (nzchar(combined_fig_text)) doc <- mk_par(doc, combined_fig_text, indent = TRUE)
    doc <- mk_blank(doc)
    cat("Added:", p[1], "\n")
  }
}
if (fig_count == 0) {
  doc <- mk_par(doc, "无可视化结果。")
}

# Save
print(doc, target = opt$output)
cat("\n================================================================================\n")
cat("Report generated:", opt$output, "\n")
cat("================================================================================\n")

# 保存 sessionInfo（logs/ 子目录，与 MachineLearn 一致）
log_dir <- file.path(result_dir_for_qs, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
session_path <- file.path(log_dir, sprintf("%s_sessionInfo.%s.txt", dataset_name, opt$timestamp))
writeLines(capture.output(sessionInfo()), session_path)
cat("Session info saved:", session_path, "\n")
