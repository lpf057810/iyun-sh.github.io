#!/usr/bin/env Rscript
# ==============================================================================
# Prognosis/scripts/r-RiskModel.R
# 风险模型分析脚本
# 支持命令行参数输入
# ==============================================================================

# 保存原始工作目录（需要在路径函数之前定义）
orig_wd <- getwd()

# 时间戳初始化
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
SCRIPT_NAME <- "Prognosis_RiskModel"

# 中文日志函数
log_info <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] [信息] %s", timestamp, message)
  cat(log_line, "\n")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

log_warn <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] [警告] %s", timestamp, message)
  cat(log_line, "\n")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

log_error <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] [错误] %s", timestamp, message)
  cat(log_line, "\n")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

# 解析命令行参数
library(optparse)

option_list <- list(
  make_option(c("--expr"), type = "character", help = "表达矩阵文件，或合并格式文件（含OS.time/OS/基因表达）"),
  make_option(c("--group"), type = "character", default = NULL, help = "分组文件（可选，合并格式时自动生成）"),
  make_option(c("--survival"), type = "character", default = NULL, help = "生存数据文件（可选，合并格式时自动拆分）"),
  make_option(c("--genes"), type = "character", default = NULL, help = "基因列表文件（可选，合并格式时自动提取）"),
  make_option(c("--lasso-coef"), type = "character", help = "Lasso系数文件"),
  make_option(c("--output"), type = "character", default = ".", help = "输出目录"),
  make_option(c("--script-dir"), type = "character", default = NULL, help = "脚本目录"),
  make_option(c("--report"), action = "store_true", default = FALSE, help = "是否生成Word报告")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$expr)) {
  stop("请提供表达矩阵文件: --expr")
}

# 保存参数
expr_file <- opt$expr
group_file <- opt$group
survival_file <- opt$survival
genes_file <- opt$genes
lasso_coef_file <- opt$`lasso-coef`

# 输出目录使用绝对路径（如果是相对路径则相对于工作目录）
output_dir <- if (grepl("^/", opt$output)) opt$output else file.path(orig_wd, opt$output)

# 加载依赖包
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ComplexHeatmap)
  library(survivalROC)
  library(qs)
})

# 统一 .qs2 保存（优先 qs2 包，回退 qs 包）
save_qs2 <- function(object, path) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, path)
    return(invisible(path))
  }
  qs::qsave(object, path)
  invisible(path)
}

# 获取脚本目录
SCRIPT_DIR <- opt$`script-dir`
if (is.null(SCRIPT_DIR)) {
  SCRIPT_DIR <- "/media/desk16/share/secure/Prognosis"
}

# 读取配置文件
read_config <- function(config_file) {
  if (!file.exists(config_file)) {
    return(NULL)
  }
  config <- list()
  lines <- readLines(config_file)
  current_section <- ""
  for (line in lines) {
    line <- trimws(line)
    if (grepl("^\\[.*\\]$", line)) {
      current_section <- gsub("\\[|\\]", "", line)
      config[[current_section]] <- list()
    } else if (grepl("^[^=]+=.*$", line)) {
      key <- trimws(gsub("=.*$", "", line))
      value <- trimws(gsub("^[^=]+=", "", line))
      value <- gsub('"', '', value)
      if (value == "true") value <- TRUE
      if (value == "false") value <- FALSE
      if (current_section != "") {
        config[[current_section]][[key]] <- value
      } else {
        config[[key]] <- value
      }
    }
  }
  return(config)
}

# 加载配置文件
config_file <- file.path(SCRIPT_DIR, "scripts", "config", "Prognosis.config.ini")
config <- read_config(config_file)

# 从配置中读取参数（如果配置存在）
if (!is.null(config)) {
  seed <- as.numeric(ifelse(is.null(config$analysis$seed), 123, config$analysis$seed))
} else {
  seed <- 123
}

# 初始化日志
log_dir <- file.path(SCRIPT_DIR, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, paste0(SCRIPT_NAME, ".", TIMESTAMP, ".log"))

# 写入日志的辅助函数
log_write <- function(text) {
  cat(text, "\n", file = log_file, append = TRUE)
}

# 获取运行环境信息
get_env_info <- function() {
  r_version <- R.version$minor
  os <- Sys.info()["sysname"]
  return(paste0("R ", r_version, ", ", os))
}

# 写入详细日志头
log_write("================================================================================")
log_write("风险模型分析日志 - Risk Model Analysis")
log_write("================================================================================")
log_write(sprintf("时间戳: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
log_write(sprintf("运行ID: %s", TIMESTAMP))
log_write(sprintf("脚本: %s v1.0", SCRIPT_NAME))
log_write(sprintf("运行环境: %s", get_env_info()))
log_write("")

# 输入参数部分
log_write("----------------------------------------")
log_write("输入参数")
log_write("----------------------------------------")
log_write(sprintf("随机种子: %s", seed))
log_write("")

# 输入文件部分
log_write("----------------------------------------")
log_write("输入文件")
log_write("----------------------------------------")

# 分析步骤部分
log_write("")
log_write("----------------------------------------")
log_write("分析步骤")
log_write("----------------------------------------")

# 转换路径为绝对路径
to_abs_path <- function(p) {
  if (is.null(p)) return(NULL)
  if (grepl("^/", p)) p else file.path(orig_wd, p)
}
expr_file <- to_abs_path(expr_file)
group_file <- to_abs_path(group_file)
survival_file <- to_abs_path(survival_file)
genes_file <- to_abs_path(genes_file)
if (!is.null(lasso_coef_file)) {
  lasso_coef_file <- to_abs_path(lasso_coef_file)
}

# 创建输出目录（使用绝对路径）
outdir <- output_dir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)

# 加载自定义函数
source(file.path(SCRIPT_DIR, "scripts", "ikl_function.R"))

# 读取数据
cat("[1/6] 读取数据...\n")

# ---- 自动检测合并格式（df.csv）并拆分 ----
auto_detect_merged <- function(expr_path) {
  df <- fread(expr_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "sample"

  # 检测特征：前3列名称是否符合 OS.time, OS
  first3 <- tolower(gsub("[^a-zA-Z0-9]", "", colnames(df)[1:3]))
  is_merged <- length(first3) >= 3 &&
               (grepl("ostime|ostime", first3[1], ignore.case = TRUE) || grepl("^os", first3[1], ignore.case = TRUE)) &&
               (grepl("^os$", first3[2], ignore.case = TRUE) || grepl("^os\\.", first3[2], ignore.case = TRUE)) &&
               (grepl("ostatus|osstatus|status", first3[3], ignore.case = TRUE) || grepl("^os$", first3[3], ignore.case = TRUE))

  # 备用检测：行数远大于列数，第一列是样本名，后面是数字列
  if (!is_merged && nrow(df) > ncol(df) - 3 && is.character(df[, 1])) {
    col2_num <- tryCatch(is.numeric(df[, 2]), error = function(e) FALSE)
    col3_num <- tryCatch(is.numeric(df[, 3]), error = function(e) FALSE)
    if (col2_num && col3_num && df[, 3][1] %in% c(0, 1)) {
      is_merged <- TRUE
    }
  }

  if (!is_merged) return(NULL)

  cat("  [自动检测] 检测到合并格式文件，正在拆分...\n")

  df[, 1] <- as.character(df[, 1])
  sample_col <- df[, 1]

  survival <- data.frame(
    sample = sample_col,
    OS.time = as.numeric(df[, 2]),
    OS = as.numeric(df[, 3]),
    stringsAsFactors = FALSE
  )

  expr_cols <- df[, -(1:3), drop = FALSE]
  tpm <- t(expr_cols)
  colnames(tpm) <- sample_col
  rownames(tpm) <- make.names(rownames(tpm), unique = TRUE)

  group <- data.frame(
    sample = sample_col,
    group = rep("Tumor", length(sample_col)),
    stringsAsFactors = FALSE
  )

  genes <- rownames(tpm)

  return(list(
    tpm = tpm,
    survival = survival,
    group = group,
    genes = genes,
    n_genes = length(genes),
    n_samples = length(sample_col)
  ))
}

# 首先尝试自动检测合并格式（独立于后续基因/Lasso逻辑）
merged <- auto_detect_merged(expr_file)

# 如果提供了 Lasso 系数文件，直接使用
if (!is.null(lasso_coef_file) && file.exists(lasso_coef_file)) {
  Lasso_Coef <- fread(lasso_coef_file)
  gene <- Lasso_Coef$gene
  cat("  使用Lasso系数文件中的", length(gene), "个基因\n")
} else if (!is.null(genes_file)) {
  gene <- fread(genes_file)$gene
} else {
  # 合并格式会自动提供gene
  if (is.null(merged)) {
    stop("请提供基因列表文件或Lasso系数文件")
  }
}

# 处理合并格式拆分结果
if (!is.null(merged)) {
  tpm <- merged$tpm
  survival <- merged$survival
  group <- merged$group
  gene_in_merged <- merged$genes
  if (!exists("Lasso_Coef")) {
    Lasso_Coef <- data.frame(gene = gene_in_merged, coefficient = 1)
  }
  cat(sprintf("  合并格式拆分完成: %d个样本, %d个基因\n", merged$n_samples, merged$n_genes))
}

# 仅在未自动拆分时才读取文件
if (!exists("tpm")) {
  tpm <- fread(expr_file) %>% column_to_rownames("SYMBOL")
}
if (!exists("survival")) {
  survival <- fread(survival_file, data.table = FALSE)
}
if (!exists("group")) {
  group <- fread(group_file)
}

# 筛选肿瘤样本（如果存在Tumor组）
if ("Tumor" %in% group$group) {
  tumor_sample <- group %>% filter(group == "Tumor") %>% pull("sample")
  tpm <- tpm[, tumor_sample]
}

# 合并数据
intersect_gene_exp <- tpm[gene, ]
sur_dat <- merge(survival, t(intersect_gene_exp), by.x = 1, by.y = 0)

# 如果没有提供 Lasso 系数，使用基因表达均值作为系数
if (is.null(lasso_coef_file) || !file.exists(lasso_coef_file)) {
  Lasso_Coef <- data.frame(gene = gene, coefficient = 1)
}

cat("[2/6] 计算风险评分...\n")
# 计算 RiskScore
lasso_dat <- sur_dat %>% dplyr::select(c("sample", Lasso_Coef$gene)) %>% column_to_rownames("sample")
f_dat <- sur_dat %>% dplyr::select(c("sample", "OS.time", "OS", Lasso_Coef$gene)) %>% column_to_rownames("sample")

risk_out <- lasso_dat %>% as.data.frame()
for (i in 1:nrow(risk_out)) {
  risk_out$RiskScore[i] <- sum(Lasso_Coef$coefficient * risk_out[i, ])
}
rs_dat <- merge(risk_out[, "RiskScore", drop = FALSE], f_dat, by = 0) %>% column_to_rownames("Row.names")

# 按最佳截断值分高低风险组
res.cut <- surv_cutpoint(rs_dat, time = "OS.time", event = "OS", variables = "RiskScore")
rs_dat$RiskGroup <- ifelse(rs_dat$RiskScore > res.cut$cutpoint$cutpoint, "High", "Low")
rs_dat <- rs_dat %>% dplyr::select(c("OS.time", "OS", "RiskScore", "RiskGroup"), everything())
fwrite(rs_dat %>% rownames_to_column("sample"), paste0(outdir, "/01.rs_dat.csv"))

cat("[3/6] 绘制KM生存曲线...\n")
km_dat <- rs_dat
km_dat$RiskGroup <- factor(km_dat$RiskGroup, levels = c("High", "Low"))
kmfit1 <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = km_dat %>% dplyr::select(-c("RiskScore")))
km_diff <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = km_dat %>% dplyr::select(-c("RiskScore")))
km_pval <- pchisq(km_diff$chisq, df = 1, lower.tail = FALSE)
km_stats <- data.frame(
  total_samples = nrow(km_dat),
  high_n = sum(km_dat$RiskGroup == "High", na.rm = TRUE),
  low_n = sum(km_dat$RiskGroup == "Low", na.rm = TRUE),
  high_death = sum(km_dat$RiskGroup == "High" & km_dat$OS == 1, na.rm = TRUE),
  low_death = sum(km_dat$RiskGroup == "Low" & km_dat$OS == 1, na.rm = TRUE),
  logrank_p = km_pval,
  risk_cutoff = as.numeric(res.cut$cutpoint$cutpoint),
  stringsAsFactors = FALSE
)
fwrite(km_stats, paste0(outdir, "/01.KM_stats.csv"))

s_surv <- ggsurvplot(kmfit1,
  pval = TRUE, conf.int = FALSE, legend.labs = c("High risk", "Low risk"),
  legend.title = "", title = "",
  risk.table = TRUE, palette = risk_color,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

s_surv$table <- s_surv$table + labs(x = "Time") + theme(panel.grid = element_blank(), legend.position = "none")
s_surv$plot <- s_surv$plot + labs(x = "Time", y = "Overall survival probability") +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

pdf(paste0(outdir, "/01.KM.pdf"), w = 8, h = 8)
print(s_surv)
dev.off()
png(paste0(outdir, "/01.KM.png"), w = 8*300, h = 8*300, res = 300)
print(s_surv)
dev.off()

cat("[4/6] 绘制时间ROC曲线...\n")
# 时间ROC曲线（固定显示1/3/5年）
plot_survival_roc <- function(data, years = c(1, 3, 5), colors = c("red", "green", "blue")) {
  plot_data <- data.frame()
  auc_map <- setNames(rep(NA_real_, length(years)), as.character(years))
  year_labels <- paste0(years, "-year")

  # 自动估计OS.time时间单位，确保1/3/5年在不同数据尺度下都可计算
  max_time <- suppressWarnings(max(as.numeric(data$OS.time), na.rm = TRUE))
  if (!is.finite(max_time)) max_time <- 0
  time_unit_factor <- if (max_time <= 50) {
    1          # 近似“年”
  } else if (max_time <= 600) {
    30.4375    # 近似“月”
  } else {
    365        # 近似“天”
  }

  for (year in years) {
    predict_time <- time_unit_factor * year
    roc_result <- tryCatch(
      survivalROC::survivalROC(
        Stime = data$OS.time,
        status = data$OS,
        marker = data$RiskScore,
        predict.time = predict_time,
        method = "KM"
      ),
      error = function(e) NULL
    )
    if (!is.null(roc_result) && length(roc_result$FP) > 0 && length(roc_result$TP) > 0) {
      tmp <- data.frame(
        FP = pmin(pmax(as.numeric(roc_result$FP), 0), 1),
        TP = pmin(pmax(as.numeric(roc_result$TP), 0), 1),
        year = paste0(year, "-year"),
        stringsAsFactors = FALSE
      )
      plot_data <- rbind(plot_data, tmp)
      auc_map[as.character(year)] <- as.numeric(roc_result$AUC)
    }
  }

  legend_labels <- paste0(
    years, "-year (AUC=",
    ifelse(is.na(auc_map[as.character(years)]), "NA", sprintf("%.3f", auc_map[as.character(years)])),
    ")"
  )
  names(legend_labels) <- year_labels

  # 强制保留1/3/5年的图例项（即使个别年份曲线不可计算）
  if (nrow(plot_data) == 0) {
    plot_data <- data.frame(
      FP = NA_real_, TP = NA_real_, year = factor(year_labels[1], levels = year_labels),
      stringsAsFactors = FALSE
    )
  } else {
    plot_data$year <- factor(plot_data$year, levels = year_labels)
  }

  # 若个别年份无法估计，仍保留其AUC标注；可绘制的曲线按可用数据展示
  p <- ggplot(plot_data, aes(x = FP, y = TP, color = year)) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(
      values = setNames(colors[1:length(years)], year_labels),
      breaks = year_labels,
      labels = legend_labels[year_labels],
      drop = FALSE
    ) +
    labs(title = "Time-dependent ROC Curve", x = "1 - Specificity", y = "Sensitivity", color = "Time") +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18), legend.position = "top")
  list(plot = p, auc_map = auc_map, curve_data = plot_data)
}

roc_res <- plot_survival_roc(data = rs_dat, years = c(1, 3, 5), colors = c("red", "green", "blue"))
p_roc <- roc_res$plot
ggsave(paste0(outdir, "/02.timeROC.pdf"), p_roc, w = 10, h = 9)
ggsave(paste0(outdir, "/02.timeROC.png"), p_roc, w = 10, h = 9, dpi = 300)
roc_auc_df <- data.frame(
  time_year = c(1, 3, 5),
  auc = as.numeric(roc_res$auc_map[c("1", "3", "5")]),
  stringsAsFactors = FALSE
)
fwrite(roc_auc_df, paste0(outdir, "/02.timeROC_auc.csv"))
if (!is.null(roc_res$curve_data) && nrow(roc_res$curve_data) > 0) {
  fwrite(roc_res$curve_data, paste0(outdir, "/02.timeROC_curve_points.csv"))
}

cat("[5/6] 绘制风险评分和生存状态图...\n")
# 生存状态图
sur_dat <- rs_dat[order(rs_dat$RiskScore), ]
sur_dat$Patient <- 1:nrow(sur_dat)
sur_dat$Status <- factor(sur_dat$OS, levels = c(1, 0), labels = c("Dead", "Alive"))
survstat_stats <- data.frame(
  total_samples = nrow(sur_dat),
  high_n = sum(sur_dat$RiskGroup == "High", na.rm = TRUE),
  low_n = sum(sur_dat$RiskGroup == "Low", na.rm = TRUE),
  dead_n = sum(sur_dat$Status == "Dead", na.rm = TRUE),
  alive_n = sum(sur_dat$Status == "Alive", na.rm = TRUE),
  high_dead_n = sum(sur_dat$RiskGroup == "High" & sur_dat$Status == "Dead", na.rm = TRUE),
  low_dead_n = sum(sur_dat$RiskGroup == "Low" & sur_dat$Status == "Dead", na.rm = TRUE),
  stringsAsFactors = FALSE
)
fwrite(survstat_stats, paste0(outdir, "/03.SurvStat_stats.csv"))

p_surv <- ggplot(sur_dat, aes(x = Patient, y = OS.time, color = Status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Dead" = risk_color[1], "Alive" = risk_color[2])) +
  labs(x = "Patients (increasing risk score)", y = "Time", color = NULL, title = "Survival Status Distribution") +
  theme_classic(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18), legend.position = "top")

ggsave(paste0(outdir, "/03.SurvStat.pdf"), p_surv, width = 10, height = 6)
ggsave(paste0(outdir, "/03.SurvStat.png"), p_surv, width = 10, height = 6, dpi = 300)

# 风险评分分布图
p_risk <- ggplot(sur_dat, aes(x = Patient, y = RiskScore, color = RiskGroup)) +
  geom_line(size = 1) +
  geom_vline(xintercept = sum(sur_dat$RiskGroup == "Low"), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Low" = risk_color[2], "High" = risk_color[1])) +
  labs(x = "Patients (increasing risk score)", y = "Risk score", color = NULL, title = "Risk Score Distribution") +
  theme_classic(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18), legend.position = "top")

ggsave(paste0(outdir, "/04.RiskScore.pdf"), p_risk, width = 10, height = 6)
ggsave(paste0(outdir, "/04.RiskScore.png"), p_risk, width = 10, height = 6, dpi = 300)
riskscore_stats <- data.frame(
  cutoff = as.numeric(res.cut$cutpoint$cutpoint),
  high_n = sum(sur_dat$RiskGroup == "High", na.rm = TRUE),
  low_n = sum(sur_dat$RiskGroup == "Low", na.rm = TRUE),
  high_score_mean = mean(sur_dat$RiskScore[sur_dat$RiskGroup == "High"], na.rm = TRUE),
  low_score_mean = mean(sur_dat$RiskScore[sur_dat$RiskGroup == "Low"], na.rm = TRUE),
  stringsAsFactors = FALSE
)
fwrite(riskscore_stats, paste0(outdir, "/04.RiskScore_stats.csv"))

cat("[6/6] 绘制热图...\n")
# 热图
heatmap_dat <- rs_dat[order(rs_dat$RiskScore, decreasing = TRUE), ]
heatmap_dat$RiskGroup <- factor(heatmap_dat$RiskGroup)
mat <- heatmap_dat %>% dplyr::select(-c("OS.time", "OS", "RiskScore", "RiskGroup")) %>% scale() %>% t()
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2

# 图5专属结果：热图基因在High/Low组的表达差异统计
gene_cols <- setdiff(colnames(heatmap_dat), c("OS.time", "OS", "RiskScore", "RiskGroup"))
if (length(gene_cols) > 0) {
  high_idx <- heatmap_dat$RiskGroup == "High"
  low_idx <- heatmap_dat$RiskGroup == "Low"
  heatmap_gene_stats <- data.frame(
    gene = gene_cols,
    high_mean = vapply(gene_cols, function(g) mean(suppressWarnings(as.numeric(heatmap_dat[high_idx, g, drop = TRUE])), na.rm = TRUE), numeric(1)),
    low_mean = vapply(gene_cols, function(g) mean(suppressWarnings(as.numeric(heatmap_dat[low_idx, g, drop = TRUE])), na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )
  heatmap_gene_stats$delta_high_minus_low <- heatmap_gene_stats$high_mean - heatmap_gene_stats$low_mean
  heatmap_gene_stats$direction <- ifelse(
    heatmap_gene_stats$delta_high_minus_low > 0, "High_up",
    ifelse(heatmap_gene_stats$delta_high_minus_low < 0, "Low_up", "No_change")
  )
  heatmap_gene_stats$abs_delta <- abs(heatmap_gene_stats$delta_high_minus_low)
  heatmap_gene_stats <- heatmap_gene_stats[order(-heatmap_gene_stats$abs_delta), ]
  fwrite(heatmap_gene_stats, paste0(outdir, "/05.heatmap_gene_stats.csv"))
}

col <- risk_color
col <- setNames(col, levels(heatmap_dat$RiskGroup))
col <- list(Group = col)
ht_opt$message <- FALSE

p_heat <- densityHeatmap(mat, title = "Distribution as heatmap", title_gp = gpar(fontsize = 14, fontface = "bold"), ylab = " ", height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = heatmap_dat$RiskGroup, col = col, annotation_name_gp = gpar(fontface = "bold")) %v%
  Heatmap(mat, row_names_gp = gpar(fontsize = 9), show_column_names = FALSE, show_row_names = TRUE, name = "expression", cluster_rows = TRUE, height = unit(6, "cm"), col = ramp_color)

pdf(paste0(outdir, "/05.HeatMap.pdf"), w = 7, h = 5)
print(p_heat)
dev.off()
png(paste0(outdir, "/05.HeatMap.png"), w = 7*300, h = 5*300, res = 300)
print(p_heat)
dev.off()

log_info(sprintf("输出目录: %s", outdir))

# 保存配置文件
log_info("保存运行配置文件...")
run_config <- list(
  timestamp = TIMESTAMP,
  script_name = SCRIPT_NAME,
  analysis_type = "Risk Model",
  parameters = list(
    seed = 123
  ),
  input_files = list(
    expr_file = opt$expr,
    group_file = opt$group,
    survival_file = opt$survival,
    lasso_coef_file = opt$`lasso-coef`
  ),
  output_dir = output_dir
)
save_qs2(run_config, path = paste0(outdir, "/RiskModel.config.qs2"))

# 创建追踪文件
trace_file <- paste0(outdir, "/.", TIMESTAMP)
file.create(trace_file)

# 保存按分析节点的 RiskModel 结果（与 MachineLearn 节点风格对齐）
node_root_dir <- if (basename(normalizePath(outdir, mustWork = FALSE)) %in% c("Cox", "RiskModel")) {
  dirname(normalizePath(outdir, mustWork = FALSE))
} else {
  normalizePath(outdir, mustWork = FALSE)
}
dataset_name <- basename(node_root_dir)
riskmodel_node <- list(
  timestamp = TIMESTAMP,
  script = SCRIPT_NAME,
  output_dir = normalizePath(outdir, mustWork = FALSE),
  input = list(
    expr_file = expr_file,
    group_file = group_file,
    survival_file = survival_file,
    genes_file = genes_file,
    lasso_coef_file = lasso_coef_file
  ),
  stats = list(
    sample_count = nrow(rs_dat),
    high_risk_count = sum(rs_dat$RiskGroup == "High", na.rm = TRUE),
    low_risk_count = sum(rs_dat$RiskGroup == "Low", na.rm = TRUE)
  ),
  result_files = list(
    risk_data = file.path(outdir, "01.rs_dat.csv"),
    km = file.path(outdir, "01.KM.pdf"),
    km_stats = file.path(outdir, "01.KM_stats.csv"),
    time_roc = file.path(outdir, "02.timeROC.pdf"),
    time_roc_auc = file.path(outdir, "02.timeROC_auc.csv"),
    time_roc_curve = file.path(outdir, "02.timeROC_curve_points.csv"),
    surv_status = file.path(outdir, "03.SurvStat.pdf"),
    surv_status_stats = file.path(outdir, "03.SurvStat_stats.csv"),
    risk_score = file.path(outdir, "04.RiskScore.pdf"),
    risk_score_stats = file.path(outdir, "04.RiskScore_stats.csv"),
    heatmap = file.path(outdir, "05.HeatMap.pdf"),
    heatmap_gene_stats = file.path(outdir, "05.heatmap_gene_stats.csv")
  )
)
save_qs2(riskmodel_node, file.path(node_root_dir, sprintf("%s_riskmodel_node.qs2", dataset_name)))

# 追溯信息
log_write("----------------------------------------")
log_write("追溯信息")
log_write("----------------------------------------")
log_write(sprintf("日志文件: %s.%s.log", SCRIPT_NAME, TIMESTAMP))
log_write(sprintf("配置文件: Prognosis.config.ini"))
log_write("")

# 分析完成
log_write("================================================================================")
log_write("分析成功完成")
log_write("================================================================================")

log_info("风险模型分析完成")

# 保存 sessionInfo（txt + qs2）
session_info_obj <- sessionInfo()
session_log_dir <- file.path(node_root_dir, "logs")
if (!dir.exists(session_log_dir)) dir.create(session_log_dir, recursive = TRUE)
writeLines(
  capture.output(session_info_obj),
  file.path(session_log_dir, sprintf("%s_sessionInfo.%s.txt", dataset_name, TIMESTAMP))
)
save_qs2(
  list(timestamp = TIMESTAMP, script = SCRIPT_NAME, session_info = session_info_obj),
  file.path(node_root_dir, sprintf("%s_sessionInfo.%s.qs2", dataset_name, TIMESTAMP))
)

# 自动生成报告
if (opt$report) {
  log_info("正在生成Word报告...")
  prog_dir <- "/media/desk16/share/secure/Prognosis"
  report_script <- file.path(prog_dir, "scripts", "generate_prognosis_report.R")

  # 计算路径时考虑: outdir是绝对路径，在setwd之后dirname(outdir)可能不准确
  # 先存储父目录路径
  results_dir <- if (grepl("^/", opt$output)) {
    # 如果output是绝对路径，从其父目录找Cox
    p <- opt$output
    if (!grepl("\\.docx$", p, ignore.case = TRUE) && !dir.exists(p)) {
      # output是目录，直接用
      p
    } else {
      dirname(p)
    }
  } else {
    # 相对路径，相对于orig_wd
    if (!grepl("\\.docx$", opt$output, ignore.case = TRUE)) {
      file.path(orig_wd, opt$output)
    } else {
      file.path(orig_wd, dirname(opt$output))
    }
  }

  cox_dir_candidate <- file.path(results_dir, "Cox")
  if (file.exists(cox_dir_candidate)) {
    system2("Rscript", args = c(
      report_script,
      "--cox-dir", cox_dir_candidate,
      "--risk-dir", ifelse(grepl("^/", opt$output), ifelse(!grepl("\\.docx$", opt$output, ignore.case = TRUE), opt$output, dirname(opt$output)), file.path(orig_wd, opt$output)),
      "--output", file.path(prog_dir, "report")
    ), stdout = TRUE, stderr = TRUE)
    log_info("Word报告生成完成")
  } else {
    log_warn("未找到Cox结果，跳过报告生成")
  }
}
