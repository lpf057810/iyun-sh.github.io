#!/usr/bin/env Rscript
# ==============================================================================
# 免疫浸润分析 - R脚本模板
# ==============================================================================
# 使用方法:
#   Rscript script.R -e expression.csv -g group.csv -o output_dir [options]
#
# 参数说明:
#   -e, --expr       : 表达矩阵文件 (必需)
#   -g, --group      : 分组文件 (必需)
#   -o, --output     : 输出目录 [default: output]
#   --arrays         : 数据是否来自芯片 [default: FALSE]
#   --control        : 对照组名称 [default: Normal]
#   --cache          : 缓存目录 [default: cache]
#   --force          : 强制重新计算 [default: FALSE]
# ==============================================================================

# ---------- 1. 依赖包 ----------
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(rstatix)
  # 根据需要添加其他依赖
})

# 设置镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# ---------- 2. 路径处理（核心！） ----------
# 获取脚本所在目录，用于加载函数库
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)

  # 方法1: 找到 --args 之后的脚本路径
  idx <- which(args == "--args")
  if (length(idx) > 0 && idx + 1 <= length(args)) {
    script_path <- args[idx + 1]
    if (file.exists(script_path)) {
      return(dirname(script_path))
    }
  }

  # 方法2: 找到包含 .R 的参数
  r_arg <- grep("\\.R$", args, value = TRUE)
  if (length(r_arg) > 0) {
    return(dirname(r_arg[1]))
  }

  # 方法3: 使用工作目录（Python已设置）
  return(getwd())
}

# 初始化脚本目录
SCRIPT_DIR <- get_script_dir()

# 加载通用函数库（如果存在）
if (file.exists(file.path(SCRIPT_DIR, "./immune_functions.R"))) {
  source(file.path(SCRIPT_DIR, "./immune_functions.R"))
}

# ---------- 3. 标准参数定义 ----------
# 所有脚本应保持一致的标准参数
get_standard_option_list <- function() {
  list(
    make_option(c("-e", "--expr"), type = "character",
                help = "Expression matrix file path (CSV format) [required]",
                default = NULL),
    make_option(c("-g", "--group"), type = "character",
                help = "Group information file path (CSV format) [required]",
                default = NULL),
    make_option(c("-o", "--output"), type = "character",
                help = "Output directory path [default: %default]",
                default = "output"),
    make_option(c("--arrays"), type = "logical",
                help = "Whether the data is from microarray [default: %default]",
                default = FALSE),
    make_option(c("--control"), type = "character",
                help = "Control group name [default: %default]",
                default = "Normal"),
    make_option(c("--cache"), type = "character",
                help = "Cache directory [default: %default]",
                default = "cache"),
    make_option(c("--force"), type = "logical",
                help = "Force re-run calculation even if cache exists [default: %default]",
                default = FALSE)
  )
}

# ---------- 4. 错误处理工具 ----------
# 安全执行函数
safe_run <- function(expr, error_msg = "执行失败", stop_on_error = TRUE) {
  result <- tryCatch({
    eval(expr)
  }, error = function(e) {
    cat("错误:", conditionMessage(e), "\n")
    if (stop_on_error) {
      quit(status = 1)
    }
    return(NULL)
  })
  return(result)
}

# 验证必需参数
validate_required_params <- function(opt, required_params) {
  missing <- c()
  for (param in required_params) {
    if (is.null(opt[[param]]) || nchar(opt[[param]]) == 0) {
      missing <- c(missing, param)
    }
  }

  if (length(missing) > 0) {
    cat("错误: 缺少必需参数:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
  }
}

# ---------- 5. 日志工具 ----------
log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_error <- function(msg) {
  cat(sprintf("[%s] 错误: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# ---------- 6. 缓存工具 ----------
# 检查并加载缓存
load_cache <- function(cache_file, force = FALSE) {
  if (!force && file.exists(cache_file)) {
    log_info(sprintf("找到缓存文件: %s", cache_file))
    return(qs::qread(cache_file))
  }
  return(NULL)
}

# 保存缓存
save_cache <- function(data, cache_file) {
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  qs::qsave(data, cache_file)
  log_info(sprintf("已保存缓存: %s", cache_file))
}
