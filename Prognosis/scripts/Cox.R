#!/usr/bin/env Rscript
# ==============================================================================
# Prognosis/scripts/r-Cox.R
# Cox回归分析 + Lasso特征选择
# 支持命令行参数输入
# ==============================================================================

# 时间戳初始化
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
SCRIPT_NAME <- "Prognosis_Cox"

# 解析命令行参数
library(optparse)

option_list <- list(
  make_option(c("--expr"), type = "character", help = "表达矩阵文件，或合并格式文件（含OS.time/OS/基因表达）"),
  make_option(c("--group"), type = "character", default = NULL, help = "分组文件（可选，合并格式时自动生成）"),
  make_option(c("--survival"), type = "character", default = NULL, help = "生存数据文件（可选，合并格式时自动拆分）"),
  make_option(c("--genes"), type = "character", default = NULL, help = "基因列表文件（可选，合并格式时自动提取）"),
  make_option(c("--output"), type = "character", default = ".", help = "输出目录"),
  make_option(c("--script-dir"), type = "character", default = NULL, help = "脚本目录")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$expr)) {
  stop("请提供表达矩阵文件: --expr")
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
  pval_cutoff <- as.numeric(ifelse(is.null(config$analysis$pval_cutoff), 0.05, config$analysis$pval_cutoff))
  cv_folds <- as.numeric(ifelse(is.null(config$analysis$cv_folds), 10, config$analysis$cv_folds))
} else {
  seed <- 123
  pval_cutoff <- 0.05
  cv_folds <- 10
}

# 初始化目录结构
MODULE_ROOT <- SCRIPT_DIR
init_module_structure <- function(root_dir) {
  required_dirs <- c("logs", "report", "result")
  for (dir_name in required_dirs) {
    dir_path <- file.path(root_dir, dir_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
}
init_module_structure(MODULE_ROOT)

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

# 创建日志文件
log_dir <- file.path(MODULE_ROOT, "logs")
log_file <- file.path(log_dir, paste0(SCRIPT_NAME, ".", TIMESTAMP, ".log"))

# 日志输出到文件和终端
log_to_file <- function(level, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] [%s] %s\n", timestamp, toupper(level), message)
  cat(log_line)
  cat(log_line, file = log_file, append = TRUE)
}

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
log_write("Cox回归分析日志 - Cox Regression Analysis")
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
log_write(sprintf("p值阈值: %s", pval_cutoff))
log_write(sprintf("Lasso交叉验证折数: %s", cv_folds))
log_write(sprintf("随机种子: %s", seed))
log_write("")

# 保存原始工作目录
orig_wd <- getwd()

# 转换输入文件为绝对路径（如果是相对路径则相对于工作目录）
to_abs_path <- function(p) {
  if (is.null(p)) return(NULL)
  if (grepl("^/", p)) p else file.path(orig_wd, p)
}
expr_file <- to_abs_path(opt$expr)
group_file <- to_abs_path(opt$group)
survival_file <- to_abs_path(opt$survival)
genes_file <- to_abs_path(opt$genes)

# 输出目录使用绝对路径（如果是相对路径则相对于工作目录）
output_dir <- if (grepl("^/", opt$output)) opt$output else file.path(orig_wd, opt$output)

# 加载依赖包（不清理环境）
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(survival)
  library(survminer)
  library(glmnet)
  library(forestplot)
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

# 设置随机种子
set.seed(seed)

# 生成运行配置
generate_config <- function() {
  config <- list(
    timestamp = TIMESTAMP,
    script_name = SCRIPT_NAME,
    analysis_type = "Cox regression + Lasso",
    parameters = list(
      pval_cutoff = pval_cutoff,
      cv_folds = cv_folds,
      seed = seed
    ),
    input_files = list(
      expr_file = opt$expr,
      group_file = ifelse(is.null(opt$group), "auto-detected", opt$group),
      survival_file = ifelse(is.null(opt$survival), "auto-detected", opt$survival),
      genes_file = ifelse(is.null(opt$genes), "auto-detected", opt$genes)
    ),
    output_dir = output_dir
  )
  return(config)
}

# 创建输出目录（使用绝对路径）
outdir <- output_dir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)

# 加载自定义函数
source(file.path(SCRIPT_DIR, "scripts", "ikl_function.R"))

# 读取数据
cat("[1/5] 读取数据...\n")

# ---- 自动检测合并格式（df.csv）并拆分 ----
# 合并格式特征：第一列为样本名，第二列为OS.time，第三列为OS，其余为基因表达
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

  # 提取样本名
  df[, 1] <- as.character(df[, 1])
  sample_col <- df[, 1]

  # 提取生存数据（假设前3列为: sample, OS.time, OS）
  survival <- data.frame(
    sample = sample_col,
    OS.time = as.numeric(df[, 2]),
    OS = as.numeric(df[, 3]),
    stringsAsFactors = FALSE
  )

  # 提取表达矩阵（剩余列）并转置为 基因x样本 格式
  expr_cols <- df[, -(1:3), drop = FALSE]
  tpm <- t(expr_cols)
  colnames(tpm) <- sample_col
  rownames(tpm) <- make.names(rownames(tpm), unique = TRUE)

  # 生成组别（全部设为Tumor）
  group <- data.frame(
    sample = sample_col,
    group = rep("Tumor", length(sample_col)),
    stringsAsFactors = FALSE
  )

  # 提取基因列表
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

# 处理输入文件
need_merged <- is.null(opt$group) || is.null(opt$survival) || is.null(opt$genes)

if (need_merged) {
  merged <- auto_detect_merged(expr_file)
  if (!is.null(merged)) {
    tpm <- merged$tpm
    survival <- merged$survival
    group <- merged$group
    gene <- merged$genes
    cat(sprintf("  拆分完成: %d个样本, %d个基因\n", merged$n_samples, merged$n_genes))
  } else {
    cat("  [警告] 未检测到合并格式，请提供 --group, --survival, --genes 参数\n")
    cat("  尝试按标准格式读取...\n")
    gene <- fread(genes_file)$gene
    tpm <- fread(expr_file) %>% column_to_rownames("SYMBOL")
    group <- fread(group_file)
    survival <- fread(survival_file, data.table = FALSE)
  }
} else {
  gene <- fread(genes_file)$gene
  tpm <- fread(expr_file) %>% column_to_rownames("SYMBOL")
  group <- fread(group_file)
  survival <- fread(survival_file, data.table = FALSE)
}

# 记录输入文件信息
log_write("----------------------------------------")
log_write("输入文件")
log_write("----------------------------------------")
if (exists("merged") && !is.null(merged)) {
  log_write(sprintf("数据格式: 合并格式（自动拆分）"))
  log_write(sprintf("原始文件: %s", opt$expr))
} else {
  log_write(sprintf("表达矩阵: %s", opt$expr))
  log_write(sprintf("分组文件: %s", opt$group))
  log_write(sprintf("生存数据: %s", opt$survival))
  log_write(sprintf("基因列表: %s", opt$genes))
}
log_write(sprintf("样本总数: %d", ncol(tpm)))
log_write(sprintf("基因总数: %d", nrow(tpm)))
log_write("")

# 筛选肿瘤样本（如果分组文件中有Tumor组，否则使用所有样本）
if ("Tumor" %in% group$group) {
  tumor_sample <- group %>% filter(group == "Tumor") %>% pull("sample")
  tpm <- tpm[, tumor_sample]
  log_write(sprintf("筛选肿瘤样本: %d个", length(tumor_sample)))
  cat("  筛选肿瘤样本:", length(tumor_sample), "个\n")
} else {
  log_write("未找到Tumor组，使用所有样本")
  cat("  未找到Tumor组，使用所有样本\n")
}

# 合并数据
intersect_gene_exp <- tpm[gene, ]
cox_exp <- merge(survival, t(intersect_gene_exp), by.x = 1, by.y = 0)

# 记录分析步骤
log_write("----------------------------------------")
log_write("分析步骤")
log_write("----------------------------------------")

cat("[2/5] 单因素Cox分析...\n")
train_dat <- cox_exp[, c("sample", "OS.time", "OS", gene)] %>% column_to_rownames("sample")
colnames_sum <- colnames(train_dat)
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
if (length(covariates) == 0) {
  stop("未找到可用于Cox分析的基因列（请检查gene列表与表达矩阵SYMBOL是否匹配）")
}

# 单因素Cox回归
univ_formulas <- sapply(covariates, function(x) as.formula(paste0("Surv(OS.time, OS)~ `", x, "`")))
univ_models <- lapply(univ_formulas, function(x) {coxph(x, data = train_dat)})

# PH假定检验
univ_cox_zph <- lapply(univ_models, function(x) {
  y <- cox.zph(x)
  p <- y$table[-nrow(y$table),][3]
  chisq <- y$table[-nrow(y$table),][1]
  df <- y$table[-nrow(y$table),][2]
  res <- c(chisq, df, p)
  names(res) <- c("chisq", "df", "p")
  return(res)
})
univ_cox_zph_df <- as.data.frame(univ_cox_zph, check.names = FALSE) %>% t() %>% as.data.frame()
univ_cox_zph_pass <- rownames(univ_cox_zph_df)[univ_cox_zph_df$p > 0.05]

# 提取结果 - 按参考代码风格: for循环 + rbind累积，直接构建干净的结果表
uniresult <- data.frame(
  Gene = character(0),
  PH_pvalue = numeric(0),
  HR = numeric(0),
  L95CI = numeric(0),
  H95CI = numeric(0),
  "CI95%" = character(0),
  p.value = numeric(0),
  stringsAsFactors = FALSE
)
for (i in seq_along(covariates)) {
  m <- univ_models[[i]]
  gene_name <- covariates[i]
  unisum <- summary(m)

  # PH检验
  test.ph <- cox.zph(m)
  PH_p <- test.ph$table[, "p"][1]

  # 提取数值
  pvalue <- round(unisum$coefficients[, 5], 4)
  HR_val <- round(unisum$coefficients[, 2], 3)
  L95 <- round(unisum$conf.int[, 3], 3)
  H95 <- round(unisum$conf.int[, 4], 3)
  CI95 <- paste0(HR_val, "[", L95, ":", H95, "]")

  uniresult <- rbind(uniresult, data.frame(
    Gene = gene_name,
    PH_pvalue = round(PH_p, 4),
    HR = HR_val,
    L95CI = L95,
    H95CI = H95,
    "CI95%" = CI95,
    p.value = pvalue,
    stringsAsFactors = FALSE
  ))
}

# 保存所有基因Cox结果
res_mod <- uniresult
if (nrow(res_mod) > 0) {
  res_mod[, c("PH_pvalue", "HR", "L95CI", "H95CI", "p.value")] <-
    lapply(res_mod[, c("PH_pvalue", "HR", "L95CI", "H95CI", "p.value"), drop = FALSE], as.numeric)
}
write.csv(res_mod, file = paste0(outdir, "/01.univariate_cox_all.csv"), row.names = FALSE)

# PH假定通过的基因
res_mod_ph <- res_mod[res_mod$PH_pvalue > 0.05, , drop = FALSE]
write.csv(res_mod_ph, file = paste0(outdir, "/01.univariate_cox_ph.csv"), row.names = FALSE)

# 筛选显著基因 (PH>0.05 & p<pval_cutoff)
res_results_pval <- res_mod_ph[res_mod_ph$p.value < pval_cutoff, , drop = FALSE]
res_results_pval <- na.omit(res_results_pval)
write.csv(res_results_pval, file = paste0(outdir, "/02.univariate_cox_result_", pval_cutoff, ".csv"), row.names = FALSE)

cat("[3/5] 绘制森林图...\n")
# 森林图 - 读回CSV确保列名干净
forest_data <- read.csv(paste0(outdir, "/02.univariate_cox_result_", pval_cutoff, ".csv"), header = TRUE, stringsAsFactors = FALSE)
forest_data$HR <- as.numeric(forest_data$HR)
forest_data$L95CI <- as.numeric(forest_data$L95CI)
forest_data$H95CI <- as.numeric(forest_data$H95CI)
forest_data$p.value <- as.numeric(forest_data$p.value)

# 取p值最小的前20个
forest_data <- forest_data[order(forest_data$p.value), ][seq_len(min(20, nrow(forest_data))), , drop = FALSE]

tabletext <- rbind(
  c("Gene", "HR(CI95%)", "P value"),
  cbind(
    as.character(forest_data$Gene),
    as.character(forest_data$CI95.),
    ifelse(forest_data$p.value < 0.001, "< 0.001", sprintf("%.3f", forest_data$p.value))
  )
)
# All vectors for forestplot must have same length (NA for header row)
mean_vec  <- c(NA, forest_data$HR)
lower_vec <- c(NA, forest_data$L95CI)
upper_vec <- c(NA, forest_data$H95CI)

pdf(paste0(outdir, "/02.univariate_cox_", pval_cutoff, "_forest.pdf"), height = 7, width = 8)
forestplot(
  labeltext = tabletext,
  graph.pos = 3,
  mean = mean_vec,
  lower = lower_vec,
  upper = upper_vec,
  zero = 1, boxsize = 0.25,
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  fn.ci_norm = fpDrawCircleCI,
  col = fpColors(box = "red", lines = "royalblue"),
  lwd.ci = 3, ci.vertices.height = 0.1,
  title = "Univariate Cox Analysis",
  xlab = "Hazard Ratio",
  xticks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1)
)
dev.off()
png(paste0(outdir, "/02.univariate_cox_", pval_cutoff, "_forest.png"), width = 8*300, height = 7*300, res = 300)
forestplot(
  labeltext = tabletext,
  graph.pos = 3,
  mean = mean_vec,
  lower = lower_vec,
  upper = upper_vec,
  zero = 1, boxsize = 0.25,
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  fn.ci_norm = fpDrawCircleCI,
  col = fpColors(box = "red", lines = "royalblue"),
  lwd.ci = 3, ci.vertices.height = 0.1,
  title = "Univariate Cox Analysis",
  xlab = "Hazard Ratio",
  xticks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1)
)
dev.off()

cat("[4/5] Lasso特征选择...\n")
# Lasso只使用通过PH检验且p<阈值的基因
lasso_genes_input <- res_results_pval$Gene
if (length(lasso_genes_input) == 0) {
  stop("没有通过筛选的基因，无法进行Lasso分析")
}

train_data <- cox_exp[, c("sample", "OS.time", "OS", lasso_genes_input)] %>% column_to_rownames("sample")
x_all <- subset(train_data, select = -c(OS, OS.time))
y_all <- subset(train_data, select = c(OS, OS.time))

# Lasso交叉验证
cvfit <- cv.glmnet(as.matrix(x_all), Surv(y_all$OS.time, y_all$OS), nfold = 10, family = "cox")
coef.min <- coef(cvfit, s = "lambda.min")

# 提取Lasso筛选的基因
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]
write.csv(data.frame(gene = lasso_geneids), paste0(outdir, "/03.lasso_genes.csv"), row.names = FALSE)
save_qs2(cvfit, path = paste0(outdir, "/cvfit.qs2"))
save_qs2(coef.min, path = paste0(outdir, "/coef.min.qs2"))

# 绘图
png(paste0(outdir, "/03.lasso.CV.png"), width = 6*300, height = 5*300, res = 300)
plot(cvfit); dev.off()
pdf(paste0(outdir, "/03.lasso.CV.pdf"), width = 6, height = 5)
plot(cvfit); dev.off()
png(paste0(outdir, "/03.lasso.Coef.png"), width = 6*300, height = 5*300, res = 300)
plot(cvfit$glmnet.fit, xvar = "lambda"); dev.off()
pdf(paste0(outdir, "/03.lasso.Coef.pdf"), width = 6, height = 5)
plot(cvfit$glmnet.fit, xvar = "lambda"); dev.off()

# 保存系数
l.coef <- coef(cvfit, s = cvfit$lambda.min, exact = FALSE)
l.coef.mat <- as.matrix(l.coef)
df.coef <- data.frame(gene = rownames(l.coef.mat), coefficient = l.coef.mat[, 1], stringsAsFactors = FALSE)
df.coef <- df.coef[df.coef$coefficient != 0, , drop = FALSE]
fwrite(df.coef, paste0(outdir, "/03.Lasso_Coefficients.csv"))

lambda_value <- data.frame(lambda.min = cvfit$lambda.min, lambda.1se = cvfit$lambda.1se)
fwrite(lambda_value, paste0(outdir, "/03.Lasso_lambda.csv"))

# 保存配置文件
log_info("保存运行配置文件...")
run_config <- generate_config()
save_qs2(run_config, path = paste0(outdir, "/Cox.config.qs2"))

# 记录统计信息
log_write("")
log_write("----------------------------------------")
log_write("统计信息")
log_write("----------------------------------------")
log_write(sprintf("分析基因数: %d", length(covariates)))
log_write(sprintf("通过PH检验基因数: %d", length(univ_cox_zph_pass)))
log_write(sprintf("显著基因数(p<0.05): %d", nrow(res_results_pval)))
log_write(sprintf("Lasso筛选基因数: %d", length(lasso_geneids)))
log_write("")

# 记录输出文件
log_write("----------------------------------------")
log_write("输出文件")
log_write("----------------------------------------")
log_write("01.univariate_cox_all.csv - 所有基因Cox结果")
log_write("01.univariate_cox_ph.csv - 通过PH检验的基因")
log_write("02.univariate_cox_0.05_forest.pdf/png - 森林图")
log_write("03.lasso.CV.pdf/png - Lasso交叉验证图")
log_write("03.lasso.Coef.pdf/png - Lasso系数图")
log_write("03.Lasso_Coefficients.csv - Lasso筛选的基因")
log_write("")

# 创建追踪文件
trace_file <- paste0(outdir, "/.", TIMESTAMP)
file.create(trace_file)

# 保存按分析节点的 Cox 结果（与 MachineLearn 节点风格对齐）
node_root_dir <- if (basename(normalizePath(outdir, mustWork = FALSE)) %in% c("Cox", "RiskModel")) {
  dirname(normalizePath(outdir, mustWork = FALSE))
} else {
  normalizePath(outdir, mustWork = FALSE)
}
dataset_name <- basename(node_root_dir)
cox_node <- list(
  timestamp = TIMESTAMP,
  script = SCRIPT_NAME,
  output_dir = normalizePath(outdir, mustWork = FALSE),
  input = list(
    expr_file = expr_file,
    group_file = group_file,
    survival_file = survival_file,
    genes_file = genes_file
  ),
  stats = list(
    analyzed_genes = length(covariates),
    ph_pass_genes = length(univ_cox_zph_pass),
    significant_genes = nrow(res_results_pval),
    lasso_selected_genes = length(lasso_geneids)
  ),
  result_files = list(
    cox_all = file.path(outdir, "01.univariate_cox_all.csv"),
    cox_ph = file.path(outdir, "01.univariate_cox_ph.csv"),
    cox_sig = file.path(outdir, paste0("02.univariate_cox_result_", pval_cutoff, ".csv")),
    lasso_genes = file.path(outdir, "03.lasso_genes.csv"),
    lasso_coef = file.path(outdir, "03.Lasso_Coefficients.csv"),
    lasso_lambda = file.path(outdir, "03.Lasso_lambda.csv")
  )
)
save_qs2(cox_node, file.path(node_root_dir, sprintf("%s_cox_node.qs2", dataset_name)))

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
