#!/usr/bin/env Rscript
# Nomogram.R
# 预测模型列线图（Nomogram）分析脚本
# 完全独立版本，不依赖 ikl_function.R
# Author: IKL
# Date: 2026-03-02
setwd("/media/desk16/share/secure/Nomogram")
# ============================================
# 0. 时间戳初始化（日志核心）
# ============================================
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
SCRIPT_NAME <- "Nomogram"

# 获取脚本所在目录
get_script_dir <- function() {
  # 命令行模式
  args <- commandArgs()
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("^--file=", "", file_arg)))
  }
  # RStudio 或交互模式
  if (interactive()) {
    return(getwd())
  }
  # 备用：使用当前工作目录
  return(getwd())
}

SCRIPT_DIR <- get_script_dir()

# ============================================
# 0.1 读取配置文件函数
# ============================================

#' 读取 INI 配置文件
read_config <- function(config_file) {
  if (!file.exists(config_file)) {
    return(NULL)
  }

  config <- list()
  current_section <- "default"

  lines <- readLines(config_file, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0 || substr(line, 1, 1) == "#") next

    if (substr(line, 1, 1) == "[" && substr(line, nchar(line), nchar(line)) == "]") {
      current_section <- substr(line, 2, nchar(line) - 1)
      config[[current_section]] <- list()
      next
    }

    if (grepl("=", line)) {
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      key <- trimws(parts[1])
      value <- trimws(paste(parts[-1], collapse = "="))

      if (tolower(value) == "true") value <- TRUE
      if (tolower(value) == "false") value <- FALSE
      if (!is.na(suppressWarnings(as.numeric(value))) && nchar(value) > 0) {
        value <- as.numeric(value)
      }

      config[[current_section]][[key]] <- value
    }
  }

  return(config)
} 

# 尝试从 config.ini 读取项目目录
CONFIG_FILE <- file.path(SCRIPT_DIR, "config", "config.ini")
config <- read_config(CONFIG_FILE)

# 自动推断项目根目录
infer_project_root <- function() {
  # 检查常见的项目目录结构
  possible_roots <- c(
    SCRIPT_DIR,                                        # 脚本直接在项目根目录
    dirname(SCRIPT_DIR),                               # 脚本在 scripts/ 下，项目根目录是父目录
    file.path(SCRIPT_DIR, "..")                       # 同上，规范化路径
  )

  # 常见的项目目录标识
  project_markers <- c("test_data", "data", "results", "logs", "report", "scripts")

  for (root in unique(possible_roots)) {
    # 标准化路径
    root <- normalizePath(root)

    # 检查是否存在项目目录标识
    has_marker <- any(sapply(project_markers, function(m) dir.exists(file.path(root, m))))

    if (has_marker) {
      return(root)
    }
  }

  # 默认返回脚本所在目录的上一级
  return(dirname(SCRIPT_DIR))
}

# 模块根目录 - 优先从配置文件读取，否则自动推断
# 注意：空字符串也视为无效，需要检查 nchar > 0
if (!is.null(config) && !is.null(config$output$project_dir) && nchar(config$output$project_dir) > 0) {
  # 配置文件中的路径可能是绝对路径或相对路径
  project_dir_config <- config$output$project_dir

  # 如果是绝对路径且存在，直接使用
  if (file.exists(project_dir_config)) {
    MODULE_ROOT <- project_dir_config
  } else if (file.exists(file.path(SCRIPT_DIR, "..", project_dir_config))) {
    # 可能是相对于脚本目录的路径
    MODULE_ROOT <- file.path(SCRIPT_DIR, "..", project_dir_config)
  } else if (file.exists(file.path(dirname(SCRIPT_DIR), project_dir_config))) {
    # 可能是相对于父目录的路径
    MODULE_ROOT <- file.path(dirname(SCRIPT_DIR), project_dir_config)
  } else {
    # 路径无效，自动推断
    MODULE_ROOT <- infer_project_root()
    cat("⚠ 配置文件中的项目路径无效，已自动推断:", MODULE_ROOT, "\n")
  }
} else {
  # 配置为空，使用自动推断
  MODULE_ROOT <- infer_project_root()
}

if (!is.null(config)) {
  cat("✓ 已加载配置文件:", CONFIG_FILE, "\n")
  cat("✓ 项目目录:", MODULE_ROOT, "\n")
} else {
  cat("⚠ 未找到配置文件，将使用命令行参数\n")
}

# ============================================
# 自动创建目录结构
# ============================================
init_module_structure <- function(root_dir) {
  required_dirs <- c("logs", "report", "results")

  for (dir_name in required_dirs) {
    dir_path <- file.path(root_dir, dir_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat(sprintf("✓ 创建目录: %s\n", dir_path))
    }
  }
}

# 初始化目录结构
init_module_structure(MODULE_ROOT)

# 禁用 future 包（rmda依赖）对命令行参数的解析
# 必须在加载任何可能触发 future 的包之前设置
Sys.setenv(R_FUTURE_DISABLE_PARALLEL = "TRUE")
options(
  future.globals.maxSize = 500 * 1024^2,
  future.rng.onMisuse = "ignore",
  future.fork.enable = FALSE
)

suppressPackageStartupMessages({
  library(optparse)
})

# ============================================
# 日志函数
# ============================================
#' 创建日志记录器
#' @param timestamp 时间戳
#' @param script_name 脚本名称
#' @param log_dir 日志目录
create_logger <- function(timestamp, script_name, log_dir) {
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  log_file <- file.path(log_dir, paste0(script_name, ".", timestamp, ".log"))

  # 创建日志连接
  con <- file(log_file, open = "wt")
  sink(con, split = TRUE)

  return(list(
    con = con,
    file = log_file,
    timestamp = timestamp
  ))
}

#' 关闭日志连接
close_logger <- function(logger) {
  if (!is.null(logger$con) && isOpen(logger$con)) {
    sink(type = "message")
    sink()
    close(logger$con)
  }
}

#' 写入日志头部
write_log_header <- function(logger, opt) {
  r_version <- R.version.string
  packages <- installed.packages()[, 1]
  pkg_versions <- paste0(c("data.table", "dplyr", "ggplot2", "rms", "ResourceSelection", "rmda", "pROC"),
                         sapply(c("data.table", "dplyr", "ggplot2", "rms", "ResourceSelection", "rmda", "pROC"),
                                function(p) {
                                  v <- try(packageVersion(p), silent = TRUE)
                                  if (inherits(v, "error")) "" else paste0(" (", as.character(v), ")")
                                }))

  cat("================================================================================\n")
  cat("Nomogram 预测模型日志\n")
  cat("================================================================================\n")
  cat(sprintf("时间戳: %s\n", logger$timestamp))
  cat(sprintf("运行ID: %s\n", logger$timestamp))
  cat(sprintf("脚本: %s.R\n", SCRIPT_NAME))
  cat(sprintf("运行环境: %s\n", r_version))
  cat("\n")
  cat("----------------------------------------\n")
  cat("输入参数\n")
  cat("----------------------------------------\n")
  cat(sprintf("数据文件: %s\n", opt$data_file))
  cat(sprintf("分组文件: %s\n", opt$group_file))
  cat(sprintf("基因文件: %s\n", opt$gene_file))
  cat(sprintf("比较组: %s\n", ifelse(is.null(opt$comparison), "自动", opt$comparison)))
  cat(sprintf("输出目录: %s\n", opt$output_dir))
  cat(sprintf("前缀: %s\n", opt$prefix))
  cat(sprintf("校准曲线bootstrap: %d\n", opt$calibrate_b))
  cat(sprintf("LRM最大迭代: %d\n", opt$lrm_maxit))
  cat("\n")
}

# ============================================
# 1. 命令行参数定义
# ============================================

option_list <- list(
  # 必需参数
  make_option(c("-d", "--data-file"), type = "character",
              help = "表达矩阵文件路径（CSV格式），第一列为基因名 [必需]"),
  make_option(c("-g", "--group-file"), type = "character",
              help = "分组文件路径（CSV格式），包含 sample 和 group 列 [必需]"),
  make_option(c("-G", "--gene-file"), type = "character",
              help = "基因列表文件路径（CSV格式，包含 gene 列）[必需]"),
  make_option(c("-c", "--comparison"), type = "character",
              help = "比较分组，格式: \"group1,group2\" [可选]；不指定时自动提取为 Normal,另一组"),

  # 输出参数
  make_option(c("-o", "--output-dir"), type = "character", default = ".",
              help = "输出目录 [默认: 当前目录]"),
  make_option(c("-p", "--prefix"), type = "character",
              help = "输出文件前缀（默认自动从 -d 推断，如 TCGA-LUAD）"),

  # 分析参数
  make_option(c("--calibrate-b"), type = "integer", default = 1000,
              help = "校准曲线 bootstrap 次数 [默认: 1000]"),
  make_option(c("--lrm-maxit"), type = "integer", default = 1000,
              help = "LRM 模型最大迭代次数 [默认: 1000]"),

  # 显示控制
  make_option(c("-q", "--quiet"), action = "store_true", default = FALSE,
              help = "静默模式")
)

parser <- OptionParser(
  usage = "%prog -d <data.csv> -g <group.csv> -G <genes.csv> [options]",
  option_list = option_list,
  description = "
预测模型列线图（Nomogram）分析工具

功能:
  1. Logistic 回归模型构建（rms::lrm）
  2. Nomogram 列线图绘制
  3. 校准曲线（Calibration curve）
  4. 决策曲线分析（DCA）
  5. ROC 曲线

使用示例:
  # 查看帮助
  Rscript Nomogram.R --help

  # 基本用法（自动推断比较组）
  Rscript Nomogram.R \\
    -d ../00_RawData/TCGA-LUAD.tpm.csv \\
    -g ../00_RawData/TCGA-LUAD.group.csv \\
    -G ../07_Expression/05.gene.csv \\
    -o . -p TCGA-LUAD

  # 指定比较组
  Rscript Nomogram.R \\
    -d data.csv -g group.csv -G genes.csv \\
    -c \"Normal,Tumor\" \\
    -o output/ -p MyStudy

  # 静默模式
  Rscript Nomogram.R -d data.csv -g group.csv -G genes.csv -q
"
)

opt <- parse_args(parser, convert_hyphens_to_underscores = TRUE)

# ============================================
# 1.5 合并配置文件与命令行参数
# 命令行参数优先于配置文件
# ============================================

if (!is.null(config)) {
  # 从配置文件读取参数（作为默认值）
  if (!is.null(config$input)) {
    if (is.null(opt$data_file) && !is.null(config$input$data_file)) {
      opt$data_file <- config$input$data_file
    }
    if (is.null(opt$group_file) && !is.null(config$input$group_file)) {
      opt$group_file <- config$input$group_file
    }
    if (is.null(opt$gene_file) && !is.null(config$input$gene_file)) {
      opt$gene_file <- config$input$gene_file
    }
  }

  # 注意：不再将 opt$output_dir 设置为 project_dir
  # result_dir 的设置会根据 config$directories$result 自动处理

  if (!is.null(config$analysis)) {
    if (is.null(opt$calibrate_b) || opt$calibrate_b == 1000) {
      if (!is.null(config$analysis$calibrate_b)) {
        opt$calibrate_b <- config$analysis$calibrate_b
      }
    }
    if (is.null(opt$lrm_maxit) || opt$lrm_maxit == 1000) {
      if (!is.null(config$analysis$lrm_maxit)) {
        opt$lrm_maxit <- config$analysis$lrm_maxit
      }
    }
    if (is.null(opt$comparison) && !is.null(config$analysis$comparison) && nchar(config$analysis$comparison) > 0) {
      opt$comparison <- config$analysis$comparison
    }
    if (is.null(opt$prefix) && !is.null(config$analysis$prefix) && nchar(config$analysis$prefix) > 0) {
      opt$prefix <- config$analysis$prefix
    }
  }

  # 显示配置来源
  cat("✓ 已合并配置文件参数\n")
}

# ============================================
# 2. 参数验证与路径解析
# ============================================

# 路径解析函数：自动处理相对/绝对路径
resolve_path <- function(path, project_root, script_dir) {
  if (is.null(path) || path == "") return(NULL)

  # 如果是绝对路径且存在，直接返回
  if (file.exists(path)) {
    return(path)
  }

  # 尝试相对于项目根目录
  rel_path1 <- file.path(project_root, path)
  if (file.exists(rel_path1)) {
    return(rel_path1)
  }

  # 尝试相对于脚本目录
  rel_path2 <- file.path(script_dir, "..", path)
  if (file.exists(rel_path2)) {
    return(rel_path2)
  }

  # 尝试相对于脚本父目录
  rel_path3 <- file.path(dirname(script_dir), path)
  if (file.exists(rel_path3)) {
    return(rel_path3)
  }

  # 返回原始路径（让后续的错误处理显示原始路径）
  return(path)
}

# 解析输入文件路径
if (!is.null(opt$data_file)) {
  opt$data_file <- resolve_path(opt$data_file, MODULE_ROOT, SCRIPT_DIR)
}
if (!is.null(opt$group_file)) {
  opt$group_file <- resolve_path(opt$group_file, MODULE_ROOT, SCRIPT_DIR)
}
if (!is.null(opt$gene_file)) {
  opt$gene_file <- resolve_path(opt$gene_file, MODULE_ROOT, SCRIPT_DIR)
}

if (is.null(opt$data_file) || is.null(opt$group_file) || is.null(opt$gene_file)) {
  print_help(parser)
  stop("\n错误: 必须指定 -d, -g, -G 参数", call. = FALSE)
}

if (!file.exists(opt$data_file)) {
  stop(sprintf("\n错误: 数据文件不存在: %s\n  请检查路径或使用绝对路径", opt$data_file), call. = FALSE)
}

if (!file.exists(opt$group_file)) {
  stop(sprintf("\n错误: 分组文件不存在: %s\n  请检查路径或使用绝对路径", opt$group_file), call. = FALSE)
}

if (!file.exists(opt$gene_file)) {
  stop(sprintf("\n错误: 基因文件不存在: %s\n  请检查路径或使用绝对路径", opt$gene_file), call. = FALSE)
}

# 自动推断前缀
infer_prefix_from_input <- function(data_file) {
  name <- basename(data_file)
  # 尝试匹配 GSE 或 TCGA 格式
  hit <- regmatches(name, regexpr("(GSE[0-9]+|TCGA-[A-Z0-9]+)", name, ignore.case = TRUE))
  if (length(hit) == 1 && !is.na(hit) && nzchar(hit)) {
    return(toupper(hit))
  }
  # 否则使用文件名前缀
  base <- sub("\\.[^.]*$", "", name)
  base <- sub("_(tpm|counts|expression)$", "", base, ignore.case = TRUE)
  if (nzchar(base)) {
    return(base)
  }
  "Nomogram"
}

if (is.null(opt$prefix) || trimws(opt$prefix) == "") {
  opt$prefix <- infer_prefix_from_input(opt$data_file)
}
opt$prefix <- basename(trimws(opt$prefix))

# ============================================
# 3. 加载依赖包
# ============================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)

  # Nomogram 核心包
  library(rms)
  # library(regplot)  # Using rms::nomogram instead
  library(ResourceSelection)
  library(rmda)
  library(pROC)
})

# 统一 .qs2 保存（优先 qs2 包，回退 qs 包）
save_qs2 <- function(object, path) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, path)
    return(invisible(path))
  }
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(object, path)
    return(invisible(path))
  }
  stop("缺少 qs2/qs 包，无法保存 .qs2 文件", call. = FALSE)
}

set.seed(123)

# ============================================
# 4. 内嵌函数定义
# ============================================

## 4.1 工具函数 ----

#' 创建目录（如果不存在）
dir_create <- function(x, recursive = FALSE) {
  if (!dir.exists(x)) {
    dir.create(x, recursive = recursive)
  }
}

#' 保存绘图为 PDF/PNG
save_plot <- function(filename, plot, outdir = ".", width = 7, height = 7,
                      bg = 'white', family = "sans", units = "in", res = 600) {

  is_gg <- inherits(plot, "ggplot")

  # 标准化路径（去除末尾斜杠，避免双斜杠）
  outdir <- gsub("/+$", "", outdir)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  base_name <- tools::file_path_sans_ext(filename)
  file_path_pdf <- file.path(outdir, paste0(base_name, ".pdf"))
  file_path_png <- file.path(outdir, paste0(base_name, ".png"))

  draw_plot <- function() {
    if (is_gg) {
      # ggplot2 使用 theme 设置字体
      print(plot)
    } else {
      if (is.function(plot)) {
        plot()
      } else {
        print(plot)
      }
    }
  }

  # 设置全局字体（Arial）
  par(family = family)

  # 保存 PDF - 设置全局字体
  pdf(file = file_path_pdf, width = width, height = height,
      bg = bg, family = family)
  draw_plot()
  dev.off()

  # 保存 PNG - 设置全局字体
  png(file = file_path_png, width = width, height = height,
      bg = bg, units = units, res = res, family = family)
  draw_plot()
  dev.off()

  # 重置字体设置
  par(family = "")

  cat("✓ 已保存 PDF 和 PNG:\n", file_path_pdf, "\n", file_path_png, "\n")
}

## 4.2 颜色定义 ----

group_color <- c("firebrick2", "#386CB0", "orange", "seagreen", "#BC80BD", "#17BECFFF")

## 4.3 数据处理函数 ----

#' 准备 Nomogram 分析数据
#'
#' @param df 表达矩阵（基因 × 样本）
#' @param df.group 分组信息（包含 sample 和 group 列）
#' @param genes 目标基因列表
#' @param comparison 比较组（c(control, treatment)）
#' @return list(data = 准备好的数据框, ddist = datadist 对象)
prepare_nomogram_data <- function(df, df.group, genes, comparison) {

  # 转置表达矩阵
  tdf <- as.data.frame(t(df))

  # 提取目标基因表达
  data <- tdf[df.group$sample, genes, drop = FALSE]

  # 创建二分类结局变量
  df.group$y <- ifelse(df.group$group == comparison[2], 1, 0)

  # 合并数据
  d <- merge(data, df.group, by.x = "row.names", by.y = "sample")
  rownames(d) <- d$Row.names
  d$Row.names <- NULL

  # 创建 datadist
  ddist <- rms::datadist(d)
  # 保存到全局环境以确保lrm可以访问
  assign("ddist", ddist, envir = .GlobalEnv)
  options(datadist = ddist)

  return(list(data = d, ddist = ddist))
}

## 4.4 模型构建函数 ----

#' 构建 LRM 模型
#'
#' @param data 准备好的数据框
#' @param genes 目标基因列表
#' @param maxit 最大迭代次数
#' @return lrm 模型对象
build_lrm_model <- function(data, genes, maxit = 1000) {

  formula1 <- as.formula(paste0('y ~ ', paste(genes, sep = '', collapse = ' + ')))

  # 尝试使用lrm，如果失败则使用glm作为后备
  fit <- tryCatch({
    lrm(formula1, data = data, x = TRUE, y = TRUE, maxit = maxit, se.fit = TRUE)
  }, error = function(e) {
    message("lrm failed, trying glm: ", e$message)
    # 使用glm作为后备
    glm_fit <- glm(formula1, data = data, family = binomial())
    # 转换为类似lrm的对象
    structure(list(
      coefficients = glm_fit$coefficients,
      stats = c(C = NA, R2 = NA),
      call = glm_fit$call,
      class = "glm"
    ), class = "try-error")
  })

  return(fit)
}

## 4.5 绘图函数 ----

#' 绘制 Nomogram 列线图
#'
#' @param fit lrm 模型对象
#' @param comparison 比较组
#' @param outdir 输出目录
#' @param prefix 文件前缀
plot_nomogram <- function(fit, comparison, outdir, prefix) {

  # 检查模型是否有效
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    message("Model fitting failed, skipping nomogram plot")
    # 创建一个占位图
    png_file <- file.path(outdir, "01.nomogram.png")
    png(png_file, width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, "Nomogram not available\n(Model fitting failed)", cex = 1.5)
    dev.off()
    return(invisible(NULL))
  }

  # 检查是否是lrm对象
  if (!inherits(fit, "lrm")) {
    message("Model is not lrm object, skipping nomogram")
    return(invisible(NULL))
  }

  # 创建 nomogram 对象
  nom <- tryCatch({
    nomogram(fit,
            fun = plogis,
            funlabel = paste0("Risk of ", comparison[2]),
            lp = FALSE,
            fun.at = c(.001, .999))
  }, error = function(e) {
    message("Nomogram creation failed: ", e$message)
    return(NULL)
  })

  if (is.null(nom)) {
    return(list(nomogram_obj = NULL))
  }

  # 绘图函数
  plot_func <- function() {
    plot(nom, cex.axis = 1.2, cex.var = 1.7)
  }

  # 保存图形
  save_plot(paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "01.nomogram.pdf"),
            plot_func,
            outdir = outdir,
            width = 10, height = 8)

  # 返回 nomogram 对象，供系数表计算 Points 使用
  return(list(nomogram_obj = nom))
}

#' 绘制校准曲线
#'
#' @param fit lrm 模型对象
#' @param data 数据框
#' @param outdir 输出目录
#' @param prefix 文件前缀
#' @param B bootstrap 次数
#' @return Hosmer-Lemeshow 检验 p 值
plot_calibrate <- function(fit, data, outdir, prefix, B = 1000) {

  # 捕获校准曲线计算的警告和错误
  calibrate_result <- tryCatch({
    # 使用suppressWarnings捕获警告
    suppressWarnings({
      cal1 <- calibrate(fit, cmethod = 'KM', method = 'boot', B = B)
    })
    # 检查是否有警告
    warns <- warnings()
    if (!is.null(warns) && length(warns) > 0) {
      warn_msg <- paste(names(warns), warns, sep = ": ", collapse = "\n")
      if (grepl("converge|singular", warn_msg, ignore.case = TRUE)) {
        list(success = FALSE, message = warn_msg, type = "convergence")
      } else {
        list(success = TRUE, cal = cal1, warning = warn_msg)
      }
    } else {
      list(success = TRUE, cal = cal1)
    }
  }, error = function(e) {
    list(success = FALSE, message = e$message, type = "error")
  })

  # 如果校准失败，返回NA并给出提示
  if (!calibrate_result$success) {
    cat("⚠ 校准曲线计算失败\n")
    if (calibrate_result$type == "convergence") {
      cat("  原因: 样本量不足或基因数过多，导致bootstrap无法收敛\n")
      cat("  建议: 减少基因数量或增加样本量\n")
    } else {
      cat(paste0("  错误: ", calibrate_result$message, "\n"))
    }
    cat("\n")
    # 创建一个空的或简单的校准图
    png_file <- paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "02.calibrate.png")
    png(file.path(outdir, png_file), width = 800, height = 600)
    plot.new()
    text(0.5, 0.5, "Calibration curve\nfailed to converge\n(reduce genes or add samples)",
         cex = 1.5, adj = 0.5)
    dev.off()

    pdf_file <- paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "02.calibrate.pdf")
    pdf(file.path(outdir, pdf_file), width = 8, height = 6)
    plot.new()
    text(0.5, 0.5, "Calibration curve\nfailed to converge\n(reduce genes or add samples)",
         cex = 1.5, adj = 0.5)
    dev.off()

    return(NA)
  }

  cal1 <- calibrate_result$cal

  # Hosmer-Lemeshow 检验
  hl1 <- ResourceSelection::hoslem.test(data$y, predict(fit, data), g = 10)
  hl2 <- stats::resid(fit, "gof")
  pval <- signif(hl2[5], 3)

  # 绘图函数
  plot_func <- function() {
    par(mar = c(6, 5, 2, 2))
    plot(cal1, lwd = 2, lty = 1,
         cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, cex.sub = 1.2,
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Predicted probability",
         ylab = "Actual probability",
         col = c("#00468BFF", "#ED0000FF", "#42B540FF"),
         legend = FALSE)
    lines(cal1[, c(1:3)], type = "l", lwd = 2, pch = 16, col = c("#00468BFF"))
    abline(0, 1, lty = 3, lwd = 2)
    legend(x = .6, y = .4, legend = c("Apparent", "Bias-corrected", "Ideal"),
           lty = c(1, 1, 2), lwd = 2, col = c("#00468BFF", "black", "black"), bty = "n")
    text(x = 0.2, y = 0.8, "Hosmer-Lemeshow ")
    text(x = 0.4, y = 0.8, as.expression(bquote(italic('p') == .(pval))))
  }

  # 计算校准MAE（Mean Absolute Error，偏差校正后的预测值与理想线的平均绝对误差）
  mae_val <- round(mean(abs(cal1[,"predy"] - cal1[,"calibrated.corrected"]), na.rm = TRUE), 3)

  # 保存图形
  save_plot(paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "02.calibrate.pdf"),
            plot_func,
            outdir = outdir,
            width = 8, height = 6)

  return(list(pval = pval, mae = mae_val))
}

#' 绘制决策曲线（DCA）
#'
#' @param fit lrm 模型对象
#' @param data 数据框
#' @param genes 目标基因列表
#' @param outdir 输出目录
#' @param prefix 文件前缀
plot_dca <- function(fit, data, genes, outdir, prefix) {

  # 构建模型公式
  formula1 <- as.formula(paste0('y ~ ', paste(genes, sep = '', collapse = ' + ')))

  # 完整模型 DCA
  model <- decision_curve(formula1, data = data)

  # 单基因 DCA
  dca_list <- lapply(genes, function(x) {
    formula_str <- paste0("y ~ ", x)
    gene_formula <- as.formula(formula_str)
    decision_curve(gene_formula, data = data)
  })

  # 绘图函数
  plot_func <- function() {
    plot_decision_curve(
      c(list(model), dca_list),
      curve.names = c('model', genes),
      legend.position = "bottomleft",
      cost.benefit.axis = FALSE,
      confidence.intervals = FALSE,
      standardize = FALSE,
      lty = c(1, 4:(length(dca_list) + 4), 2, 3)
    )
  }

  # 保存图形
  save_plot(paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "03.dca.pdf"),
            plot_func,
            outdir = outdir,
            width = 8, height = 7)
  # 输出图3专属统计表，供报告自动解释
  dca_stats <- data.frame(
    metric = c("sample_count", "gene_count", "threshold_min", "threshold_max", "max_net_benefit", "threshold_at_max_nb"),
    value = c(nrow(data), length(genes), NA, NA, NA, NA),
    stringsAsFactors = FALSE
  )
  if (!is.null(model$derived.data) && nrow(model$derived.data) > 0) {
    dd <- model$derived.data
    thr_col <- colnames(dd)[grepl("threshold", colnames(dd), ignore.case = TRUE)][1]
    nb_col <- colnames(dd)[grepl("^NB$|net.benefit|net_benefit", colnames(dd), ignore.case = TRUE)][1]
    if (!is.na(thr_col) && !is.na(nb_col) && thr_col %in% colnames(dd) && nb_col %in% colnames(dd)) {
      thr <- suppressWarnings(as.numeric(dd[[thr_col]]))
      nb <- suppressWarnings(as.numeric(dd[[nb_col]]))
      ok <- is.finite(thr) & is.finite(nb)
      if (any(ok)) {
        thr_ok <- thr[ok]; nb_ok <- nb[ok]
        idx <- which.max(nb_ok)
        dca_stats$value[dca_stats$metric == "threshold_min"] <- sprintf("%.4f", min(thr_ok, na.rm = TRUE))
        dca_stats$value[dca_stats$metric == "threshold_max"] <- sprintf("%.4f", max(thr_ok, na.rm = TRUE))
        dca_stats$value[dca_stats$metric == "max_net_benefit"] <- sprintf("%.4f", nb_ok[idx])
        dca_stats$value[dca_stats$metric == "threshold_at_max_nb"] <- sprintf("%.4f", thr_ok[idx])
      }
    }
    write.csv(dd, file.path(outdir, "03.dca_curve_data.csv"), row.names = FALSE)
  }
  write.csv(dca_stats, file.path(outdir, "03.dca_stats.csv"), row.names = FALSE)
  return(invisible(dca_stats))
}

#' 绘制 ROC 曲线
#'
#' @param fit lrm 模型对象
#' @param data 数据框
#' @param comparison 比较组
#' @param outdir 输出目录
#' @param prefix 文件前缀
#' @return AUC 值
plot_roc <- function(fit, data, comparison, outdir, prefix) {

  # 预测
  predicted <- predict(fit, newdata = data)

  # ROC 曲线
  roc_curve <- roc(data$group, predicted, levels = comparison)
  auc_val <- auc(roc_curve)

  # 计算最佳阈值点
  roc_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv", "youden"))
  roc_results <- data.frame(
    Metric = c("AUC", "Sensitivity", "Specificity", "PPV", "NPV", "Youden Index", "Optimal Threshold"),
    Value = c(
      round(as.numeric(auc_val), 4),
      round(roc_coords$sensitivity, 4),
      round(roc_coords$specificity, 4),
      round(roc_coords$ppv, 4),
      round(roc_coords$npv, 4),
      round(roc_coords$youden, 4),
      round(roc_coords$threshold, 4)
    )
  )
  write.csv(roc_results, file.path(outdir, "04.roc_results.csv"), row.names = FALSE)

  # 绘图函数
  plot_func <- function() {
    par(pin = c(4, 4), mar = c(6, 6, 6, 1))
    plot(roc_curve, auc.main = TRUE, legend = FALSE, color = 'darkblue',
         xlab = "1-Specificity",
         print.auc = TRUE,
         asp = 1,
         cex.axis = 1.8,
         cex.lab = 2.0,
         cex.main = 2.0,
         main = 'model',
         col = "#FF2E63",
         font.lab = 2,
         font.main = 2,
         font.sub = 2)
  }

  # 保存图形
  save_plot(paste0(ifelse(nchar(prefix) > 0, paste0(prefix, "/"), ""), "04.roc.pdf"),
            plot_func,
            outdir = outdir,
            width = 6, height = 6)

  return(as.numeric(auc_val))
}

## 4.6 Quarto 报告生成函数 (见第6.1节) ----

# ============================================
# 5. 主程序
# ============================================

# 初始化日志
log_dir <- file.path(MODULE_ROOT, "logs")
logger <- create_logger(TIMESTAMP, SCRIPT_NAME, log_dir)
write_log_header(logger, opt)

# 设置输出目录结构（按标准化规范）
# 优先级：1. -o 参数指定 2. config$directories$result 3. 默认 result/Nomogram
# 标记用户是否通过 -o 显式指定了输出目录
user_specified_output <- !is.null(opt$output_dir) && opt$output_dir != "."

if (user_specified_output) {
  # 用户通过 -o 显式指定了输出目录
  result_dir <- file.path(opt$output_dir, "Nomogram")
} else if (!is.null(config) && !is.null(config$directories$result)) {
  # 使用配置文件中的 result 目录（相对于 project_dir）
  result_dir <- file.path(MODULE_ROOT, config$directories$result)
} else {
  # 默认使用 result/Nomogram
  result_dir <- file.path(MODULE_ROOT, "results", "Nomogram")
}

# 确保目录存在
result_dir <- normalizePath(result_dir, mustWork = FALSE)
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

# 创建隐藏追踪文件
timestamp_file <- file.path(result_dir, paste0(".", TIMESTAMP))
file.create(timestamp_file)

# 保存配置文件
config_file <- file.path(result_dir, "Nomogram.config.ini")
config_lines <- c(
  "[Analysis]",
  paste0("timestamp=", TIMESTAMP),
  paste0("script=", SCRIPT_NAME, ".R"),
  "",
  "[Parameters]",
  paste0("data_file=", opt$data_file),
  paste0("group_file=", opt$group_file),
  paste0("gene_file=", opt$gene_file),
  paste0("comparison=", ifelse(is.null(opt$comparison), "auto", opt$comparison)),
  paste0("output_dir=", result_dir),
  paste0("prefix=", opt$prefix),
  paste0("calibrate_b=", opt$calibrate_b),
  paste0("lrm_maxit=", opt$lrm_maxit)
)
writeLines(config_lines, config_file)

cat(sprintf("[%s] 开始分析...\n", format(Sys.time(), "%H:%M:%S")))

# 读取数据
if (!opt$quiet) cat("\n[读取数据]\n")

# 表达矩阵
df <- fread(opt$data_file)
df <- as.data.frame(df)

# 自动检测基因列
gene_col <- NULL
if ("SYMBOL" %in% colnames(df)) {
  gene_col <- "SYMBOL"
} else if ("GeneID" %in% colnames(df)) {
  gene_col <- "GeneID"
} else {
  gene_col <- colnames(df)[1]
}
rownames(df) <- df[[gene_col]]

# 移除基因标识列
gene_id_cols <- c("SYMBOL", "GeneID", "gene_name", "GENEID", "gene_id")
df <- df[, !colnames(df) %in% gene_id_cols, drop = FALSE]

# 分组信息
df.group <- fread(opt$group_file)
df.group <- as.data.frame(df.group)

if (!"group" %in% colnames(df.group)) {
  stop("错误: 分组文件必须包含 group 列", call. = FALSE)
}
if (!"sample" %in% colnames(df.group)) {
  stop("错误: 分组文件必须包含 sample 列", call. = FALSE)
}

# 基因列表
gene_df <- fread(opt$gene_file)
if ("gene" %in% colnames(gene_df)) {
  genes <- gene_df$gene
} else if ("SYMBOL" %in% colnames(gene_df)) {
  genes <- gene_df$SYMBOL
} else {
  genes <- gene_df[[1]]
}
genes <- as.character(genes)

# 验证表达矩阵格式
# 检查第一列是否为基因名（而非样本名）
check_expr_format <- function(df, genes) {
  # 假设第一行是基因名
  first_row_values <- df[1, ]
  n_genes_in_list <- length(genes)

  # 检查有多少个基因列表中的基因匹配到行名
  matched_in_rows <- sum(genes %in% rownames(df))

  # 检查样本数：如果列数远大于基因数，可能是基因为列的格式
  n_cols <- ncol(df)
  n_rows <- nrow(df)

  # 启发式判断：如果行数接近基因数，但列数远大于基因数，可能格式错误
  if (n_cols > n_rows * 1.5 && n_rows < 100) {
    return(list(
      valid = FALSE,
      message = paste0(
        "⚠ 检测到表达矩阵可能是「基因为列、样本为行」的格式\n",
        "  当前: ", n_rows, " 行(基因) × ", n_cols, " 列(样本)\n",
        "  要求: 基因应为行，样本应为列\n",
        "  请确保输入文件格式为: 第一列为基因名(SYMBOL)，其余列为样本表达值")
    ))
  }

  # 检查是否匹配
  if (matched_in_rows == 0 && n_genes_in_list > 0) {
    # 尝试检查是否匹配到列名
    matched_in_cols <- sum(genes %in% colnames(df))
    if (matched_in_cols > 0) {
      return(list(
        valid = FALSE,
        message = paste0(
          "⚠ 表达矩阵格式可能不正确\n",
          "  基因列表在行名中匹配: ", matched_in_rows, "/", n_genes_in_list, "\n",
          "  基因列表在列名中匹配: ", matched_in_cols, "/", n_genes_in_list, "\n",
          "  当前格式: 基因为列（需要转换为基因为行）\n",
          "  要求格式: 第一列为基因名(SYMBOL)，其余列为样本表达值")
      ))
    }
  }

  return(list(valid = TRUE, message = ""))
}

# 格式检查
format_check <- check_expr_format(df, genes)
if (!format_check$valid) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                    输入文件格式错误                               ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  cat(format_check$message)
  cat("\n\n")
  stop("请检查并修正表达矩阵格式后重试", call. = FALSE)
}

# 验证基因是否存在
missing_genes <- setdiff(genes, rownames(df))
if (length(missing_genes) > 0) {
  stop(sprintf("错误: 以下基因在表达矩阵中不存在: %s", paste(missing_genes, collapse = ", ")), call. = FALSE)
}

if (!opt$quiet) {
  cat(sprintf("✓ 表达数据: %d 基因 × %d 样本\n", nrow(df), ncol(df)))
  cat(sprintf("✓ 分组数据: %d 样本\n", nrow(df.group)))
  cat(sprintf("✓ 目标基因: %d 个\n", length(genes)))
}

# 检查样本量与基因数比例
n_genes <- length(genes)
n_samples <- nrow(df.group)
ratio <- n_samples / n_genes

if (ratio < 3) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║                      样本量警告                                 ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  cat(paste0("⚠ 当前样本量与基因数比例较低: ", n_samples, " 样本 / ", n_genes, " 基因 = ", round(ratio, 2), "\n"))
  cat("  建议: 样本数至少为基因数的 3-5 倍以避免过拟合\n")
  cat("  风险: 可能导致模型过拟合、校准曲线无法收敛\n")
  cat("  建议: 减少基因数量或增加样本量\n")
  cat("\n")
}

# 解析比较分组
unique_groups <- unique(df.group$group)

if (is.null(opt$comparison) || trimws(opt$comparison) == "") {
  # 自动模式：以 Normal 为对照组
  if ("Normal" %in% unique_groups) {
    ref_group <- "Normal"
  } else {
    # 尝试常见的对照组名称
    common_refs <- c("Control", "Normal", "Healthy", "HC")
    ref_group <- intersect(common_refs, unique_groups)
    if (length(ref_group) == 0) {
      ref_group <- unique_groups[1]
    } else {
      ref_group <- ref_group[1]
    }
  }
  other_group <- setdiff(unique_groups, ref_group)
  if (length(other_group) == 0) {
    stop("错误: 无法自动确定比较组，请使用 -c 参数手动指定", call. = FALSE)
  }
  comparison <- c(ref_group, other_group[1])
  if (!opt$quiet) {
    cat(sprintf("✓ 自动提取比较分组: %s vs %s\n", comparison[1], comparison[2]))
  }
} else {
  comparison <- trimws(strsplit(opt$comparison, ",")[[1]])
  if (length(comparison) != 2) {
    stop("错误: -c 参数必须为两个分组，格式: \"group1,group2\"", call. = FALSE)
  }
  missing_groups <- setdiff(comparison, unique_groups)
  if (length(missing_groups) > 0) {
    stop(sprintf("错误: 以下分组不存在: %s", paste(missing_groups, collapse = ", ")), call. = FALSE)
  }
}

if (!opt$quiet) {
  cat(sprintf("✓ 比较分组: %s vs %s\n", comparison[1], comparison[2]))
}

# 数据质量控制
if (!opt$quiet) cat("\n[数据质控]\n")

qc_results <- list()

# 1. 检查缺失值
missing_count <- sum(is.na(df))
if (missing_count > 0) {
  qc_results$missing <- missing_count
  cat(sprintf("⚠ 表达矩阵存在 %d 个缺失值\n", missing_count))
} else {
  qc_results$missing <- 0
  cat("✓ 表达矩阵无缺失值\n")
}

# 2. 检查异常值（表达值为负或极端高值）
# 注意：基因名为df的行名，不是列名
expr_values <- as.matrix(df[genes, , drop = FALSE])
neg_count <- sum(expr_values < 0, na.rm = TRUE)
extreme_count <- sum(expr_values > 1e10, na.rm = TRUE)  # 超过1e10视为极端值

qc_results$negative <- neg_count
qc_results$extreme <- extreme_count

if (neg_count > 0) {
  cat(sprintf("⚠ 发现 %d 个负表达值（已设置为0）\n", neg_count))
  expr_values[expr_values < 0] <- 0
}
if (extreme_count > 0) {
  cat(sprintf("⚠ 发现 %d 个极端高表达值（>1e10）\n", extreme_count))
}
if (neg_count == 0 && extreme_count == 0) {
  cat("✓ 表达值范围正常\n")
}

# 3. 样本分布检查
n_group1 <- sum(df.group$group == comparison[1])
n_group2 <- sum(df.group$group == comparison[2])
qc_results$group1 <- n_group1
qc_results$group2 <- n_group2

cat(sprintf("✓ 分组样本分布: %s=%d, %s=%d\n",
            comparison[1], n_group1, comparison[2], n_group2))

# 检查样本不平衡
if (n_group1 / n_group2 < 0.5 || n_group2 / n_group1 < 0.5) {
  cat("⚠ 警告: 两组样本比例不平衡，可能影响模型稳定性\n")
  qc_results$imbalanced <- TRUE
} else {
  qc_results$imbalanced <- FALSE
}

# 4. 样本量与特征比
n_genes <- length(genes)
n_samples <- nrow(df.group)
ratio <- n_samples / n_genes
qc_results$ratio <- ratio
qc_results$overfitting_risk <- ratio < 3

if (ratio < 3) {
  cat(sprintf("⚠ 过拟合风险: 样本/基因比 = %.2f (< 3)\n", ratio))
} else {
  cat(sprintf("✓ 样本量充足: 样本/基因比 = %.2f\n", ratio))
}

# 汇总质控结果
cat("\n--- 质控摘要 ---\n")
cat(sprintf("缺失值: %d | 负值: %d | 极端值: %d\n",
            qc_results$missing, qc_results$negative, qc_results$extreme))
cat(sprintf("样本比: %.2f | 组间平衡: %s | 过拟合风险: %s\n",
            ratio,
            ifelse(qc_results$imbalanced, "否", "是"),
            ifelse(qc_results$overfitting_risk, "高", "低")))
cat("-------------------\n\n")

# 过滤分组数据
df.group <- df.group[df.group$group %in% comparison, ]

# 准备 Nomogram 数据
if (!opt$quiet) cat("\n[准备数据]\n")
prep_result <- prepare_nomogram_data(df, df.group, genes, comparison)
d <- prep_result$data

if (!opt$quiet) {
  cat(sprintf("✓ 分析样本: %d\n", nrow(d)))
  cat(sprintf("  - %s: %d\n", comparison[1], sum(d$group == comparison[1])))
  cat(sprintf("  - %s: %d\n", comparison[2], sum(d$group == comparison[2])))
}

# 构建模型
if (!opt$quiet) cat("\n[构建 LRM 模型]\n")
fit <- build_lrm_model(d, genes, maxit = opt$lrm_maxit)

if (!opt$quiet) {
  cat(sprintf("✓ 模型 C-index: %.4f\n", fit$stats["C"]))
  cat(sprintf("✓ 模型 R²: %.4f\n", fit$stats["R2"]))
}

# 模型质量评估
c_index <- fit$stats["C"]
r2 <- fit$stats["R2"]

if (c_index > 0.95) {
  cat(sprintf("⚠ 警告: C-index = %.4f > 0.95，可能存在过拟合\n", c_index))
}
if (r2 > 0.9) {
  cat(sprintf("⚠ 警告: R² = %.4f > 0.9，可能存在过拟合\n", r2))
}
if (c_index > 0.95 || r2 > 0.9) {
  cat("  建议: 增加样本量或减少特征数以提高模型泛化能力\n")
}

# 保存模型系数（完整系数，包括截距），同时计算OR和95%CI，以及变量范围和得分
# 先创建 nomogram 对象用于提取 Points 刻度
nom_result <- plot_nomogram(fit, comparison, result_dir, "")
nom_obj <- if (!is.null(nom_result)) nom_result$nomogram_obj else NULL

coef_names <- names(fit$coefficients)
coef_values <- fit$coefficients

# 计算标准误和95%CI（排除截距）
if (!is.null(fit$se)) {
  se_vec <- fit$se
} else if (!is.null(fit$var) && nrow(fit$var) > 0) {
  se_vec <- sqrt(diag(fit$var))
} else {
  se_vec <- rep(NA, length(coef_values))
  names(se_vec) <- coef_names
}
z <- qnorm(0.975)

# 为非截距变量计算OR和95%CI
is_intercept <- coef_names == "Intercept"
or_vals <- round(exp(coef_values), 4)
or_low <- rep(NA, length(coef_values))
or_high <- rep(NA, length(coef_values))
or_low[!is_intercept] <- round(exp(coef_values[!is_intercept] - z * se_vec[!is_intercept]), 4)
or_high[!is_intercept] <- round(exp(coef_values[!is_intercept] + z * se_vec[!is_intercept]), 4)

# 计算每个变量的表达范围和得分贡献
var_names <- coef_names[!is_intercept]
var_coefs <- coef_values[!is_intercept]

# 表达范围
expr_range <- sapply(var_names, function(g) {
  rng <- range(d[[g]], na.rm = TRUE)
  paste0(round(rng[1], 2), "-", round(rng[2], 2))
})

# 从列线图对象提取每个基因轴的 Points 范围
# nomogram() 返回的每个元素是 named list: CXCL8 | Xbeta | points
nom_pts_min <- sapply(var_names, function(v) NA)
nom_pts_max <- sapply(var_names, function(v) NA)
if (!is.null(nom_obj)) {
  for (v in var_names) {
    if (v %in% names(nom_obj)) {
      el <- nom_obj[[v]]
      if (is.list(el) && "points" %in% names(el)) {
        vals <- as.numeric(el[["points"]])
        nom_pts_min[v] <- min(vals, na.rm = TRUE)
        nom_pts_max[v] <- max(vals, na.rm = TRUE)
      }
    }
  }
}

# 构建 Points 范围字符串（如 "0-20"），小数值保留1位
make_pts_range <- function(min_pt, max_pt) {
  if (is.na(min_pt) || is.na(max_pt)) return(NA)
  paste0(round(min_pt, 1), "-", round(max_pt, 1))
}
points_ranges <- mapply(make_pts_range, nom_pts_min, nom_pts_max)

# 计算相对贡献：每个基因的 Points / 总 Points * 100
# 使用 nomogram 实际刻度，确保 Points 之和为 100
if (!all(is.na(nom_pts_max))) {
  total_pts <- sum(nom_pts_max, na.rm = TRUE)
  if (total_pts > 0) {
    points_contrib <- round(nom_pts_max / total_pts * 100, 1)
  } else {
    points_contrib <- rep(NA, length(var_names))
  }
} else {
  # 回退：使用系数 * 表达范围估算相对贡献
  var_range_vals <- sapply(var_names, function(g) {
    rng <- range(d[[g]], na.rm = TRUE)
    rng[2] - rng[1]
  })
  raw_points <- abs(var_coefs) * var_range_vals
  total_raw_points <- sum(raw_points)
  if (total_raw_points > 0) {
    points_contrib <- round(raw_points / total_raw_points * 100, 1)
  } else {
    points_contrib <- rep(NA, length(var_names))
  }
}

# 构建系数表
# Points 列：显示范围（0-20 格式），供表格阅读
# Points_pct 列：相对贡献百分比（数值），供报告解读
coef_df <- data.frame(
  Variable = coef_names,
  Coefficient = round(coef_values, 4),
  OR = or_vals,
  OR_lower = or_low,
  OR_upper = or_high,
  Range = c("-", expr_range),
  Points = c("-", points_ranges),
  Points_pct = c(NA, points_contrib),
  stringsAsFactors = FALSE
)

# 保存所有系数（包括截距Intercept）
write.csv(coef_df, file = file.path(result_dir, "05.model_coefficients.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(coef_df, file = file.path(result_dir, "01.nomogram_stats.csv"),
          row.names = FALSE, quote = FALSE)
if (!opt$quiet) {
  cat(sprintf("✓ 模型系数已保存: %s\n", file.path(result_dir, "05.model_coefficients.csv")))
}

# 绘制 Nomogram
if (!opt$quiet) cat("\n[绘制 Nomogram 列线图]\n")
# 列线图已在上方创建并保存
if (!opt$quiet) cat("  (列线图已生成)\n")
# 绘制校准曲线
if (!opt$quiet) cat("\n[绘制校准曲线]\n")
calib_result <- plot_calibrate(fit, d, result_dir, "", B = opt$calibrate_b)
if (is.list(calib_result)) {
  hl_pval <- calib_result$pval
  calib_mae <- calib_result$mae
} else {
  hl_pval <- calib_result
  calib_mae <- NA
}
calib_stats <- data.frame(
  metric = c("hosmer_lemeshow_p", "calibration_mae", "bootstrap_B"),
  value = c(hl_pval, calib_mae, opt$calibrate_b),
  stringsAsFactors = FALSE
)
write.csv(calib_stats, file.path(result_dir, "02.calibrate_stats.csv"), row.names = FALSE)
if (!opt$quiet) {
  cat(sprintf("✓ Hosmer-Lemeshow p = %s\n", hl_pval))
  if (!is.na(calib_mae)) {
    cat(sprintf("✓ 校准MAE = %.3f\n", calib_mae))
  }
}

# 绘制 DCA
if (!opt$quiet) cat("\n[绘制决策曲线 DCA]\n")
plot_dca(fit, d, genes, result_dir, "")

# 绘制 ROC
if (!opt$quiet) cat("\n[绘制 ROC 曲线]\n")
auc_val <- plot_roc(fit, d, comparison, result_dir, "")
if (!opt$quiet) {
  cat(sprintf("✓ AUC = %.4f\n", auc_val))
}

# 保存按分析节点拆分的 .qs2 结果（与 MachineLearn 节点风格对齐）
dataset_name <- if (!is.null(opt$prefix) && nzchar(trimws(opt$prefix))) opt$prefix else basename(result_dir)
save_qs2(
  list(
    timestamp = TIMESTAMP,
    script = SCRIPT_NAME,
    input = list(
      data_file = opt$data_file,
      group_file = opt$group_file,
      gene_file = opt$gene_file,
      comparison = comparison
    ),
    dimensions = list(
      input_samples = nrow(df.group),
      input_genes = nrow(df),
      target_genes = length(genes),
      analysis_samples = nrow(d)
    )
  ),
  file.path(result_dir, sprintf("%s_data_node.qs2", dataset_name))
)
save_qs2(
  list(
    timestamp = TIMESTAMP,
    model = fit,
    coefficients = coef_df,
    c_index = fit$stats["C"],
    r2 = fit$stats["R2"]
  ),
  file.path(result_dir, sprintf("%s_model_node.qs2", dataset_name))
)
save_qs2(
  list(
    timestamp = TIMESTAMP,
    metrics = list(
      auc = auc_val,
      c_index = fit$stats["C"],
      hl_pval = hl_pval,
      calib_mae = calib_mae
    ),
    files = list(
      roc_results = file.path(result_dir, "04.roc_results.csv"),
      model_coefficients = file.path(result_dir, "05.model_coefficients.csv")
    )
  ),
  file.path(result_dir, sprintf("%s_evaluation_node.qs2", dataset_name))
)

# ============================================
# 6. 输出汇总
# ============================================

# 收集输出文件
output_files <- list.files(result_dir, pattern = "\\.(pdf|png|csv)$", full.names = TRUE)

if (!opt$quiet) {
  cat("\n[运行摘要]\n")
  cat(sprintf("- 比较组: %s vs %s\n", comparison[1], comparison[2]))
  cat(sprintf("- 目标基因: %s\n", paste(genes, collapse = ", ")))
  cat(sprintf("- 模型 C-index: %.4f\n", fit$stats["C"]))
  cat(sprintf("- AUC: %.4f\n", auc_val))
  cat(sprintf("- Hosmer-Lemeshow p: %s\n", hl_pval))
  cat("- 输出文件:\n")
  for (f in output_files) {
    cat(sprintf("    * %s\n", basename(f)))
  }
}

# ─────────────────────────────────────────
# 保存 Config 快照（带时间戳）
# ─────────────────────────────────────────
config_snapshot_latest    <- file.path(result_dir, paste0("Nomogram_config_snapshot.ini"))
config_snapshot_stamped  <- file.path(result_dir, paste0("Nomogram_config_snapshot_", TIMESTAMP, ".ini"))

snapshot_content <- c(
  paste0("# Nomogram Config Snapshot - ", TIMESTAMP),
  paste0("# Generated: ", Sys.time()),
  paste0("# Data: ", opt$data_file),
  paste0("# Genes: ", paste(genes, collapse = ",")),
  paste0("# Comparison: ", paste(comparison, collapse = " vs ")),
  paste0(""),
  paste0("[run]"),
  paste0("timestamp = ", TIMESTAMP),
  paste0("dataset = ", basename(normalizePath(opt$data_file))),
  paste0("data_file = ", normalizePath(opt$data_file)),
  paste0("group_file = ", normalizePath(opt$group_file)),
  paste0("gene_file = ", normalizePath(opt$gene_file)),
  paste0("comparison = ", paste(comparison, collapse = ",")),
  paste0("calibrate_b = ", opt$calibrate_b),
  paste0("lrm_maxit = ", opt$lrm_maxit),
  paste0("output_dir = ", normalizePath(result_dir)),
  paste0("c_index = ", fit$stats["C"]),
  paste0("r2 = ", fit$stats["R2"]),
  paste0("auc = ", auc_val),
  paste0("hl_pval = ", hl_pval),
  paste0("calib_mae = ", calib_mae)
)
writeLines(snapshot_content, config_snapshot_latest)
writeLines(snapshot_content, config_snapshot_stamped)
if (!opt$quiet) cat(sprintf("✓ 配置快照已保存: %s\n", basename(config_snapshot_latest)))

# ─────────────────────────────────────────
# 保存 SessionInfo
# ─────────────────────────────────────────
logs_dir <- file.path(MODULE_ROOT, "logs")
sessioninfo_file <- file.path(logs_dir, paste0("Nomogram_sessionInfo_", TIMESTAMP, ".txt"))

si_capture <- paste0(
  capture.output(print(sessionInfo()), file = NULL),
  collapse = "\n"
)
sink(sessioninfo_file)
print(sessionInfo())
sink()
if (!opt$quiet) cat(sprintf("✓ SessionInfo 已保存: %s\n", basename(sessioninfo_file)))
save_qs2(
  list(timestamp = TIMESTAMP, script = SCRIPT_NAME, session_info = sessionInfo()),
  file.path(result_dir, sprintf("%s_sessionInfo.%s.qs2", dataset_name, TIMESTAMP))
)

# ─────────────────────────────────────────
# 生成报告（使用 Quarto）
# ─────────────────────────────────────────
report_dir <- file.path(MODULE_ROOT, "report")
if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
}
# 报告文件名添加 nomogram_ 前缀
report_file <- file.path(report_dir, paste0("nomogram_", TIMESTAMP, ".docx"))

# 准备报告参数
report_params <- list(
  TIMESTAMP = TIMESTAMP,
  DATE = format(Sys.Date(), "%Y-%m-%d"),
  DATA_FILE = opt$data_file,
  GROUP_FILE = opt$group_file,
  GENE_FILE = opt$gene_file,
  COMPARISON = paste(comparison[1], comparison[2], sep = " vs "),
  DATASET_NAME = basename(normalizePath(opt$data_file)),
  GENES = paste(genes, collapse = ", "),
  CALIBRATE_B = opt$calibrate_b,
  LRM_MAXIT = opt$lrm_maxit,
  C_INDEX = sprintf("%.4f", fit$stats["C"]),
  R2 = sprintf("%.4f", fit$stats["R2"]),
  HL_PVAL = hl_pval,
  CALIB_MAE = calib_mae,
  AUC_VAL = sprintf("%.4f", auc_val),
  OUTPUT_DIR = normalizePath(result_dir),
  COEFFICIENTS_TABLE = "See result files"
)

# 生成 Word 报告 (使用 generate_nomogram_report.R)
cat("正在生成报告...\n")
report_script <- file.path(MODULE_ROOT, "scripts", "generate_nomogram_report.R")

# 保存参数到临时文件供报告脚本使用
params_file <- file.path(result_dir, "_report_params.RData")
save(list = c("report_params", "comparison", "genes"), file = params_file)
save_qs2(
  list(report_params = report_params, comparison = comparison, genes = genes),
  file.path(result_dir, sprintf("%s_report_input.qs2", dataset_name))
)

if (file.exists(report_script)) {
  # 调用报告生成脚本
  report_cmd <- sprintf(
    "Rscript %s -i %s -o %s -T %s",
    report_script,
    result_dir,
    report_file,
    TIMESTAMP
  )
  report_result <- system(report_cmd)
  if (report_result == 0 && file.exists(report_file)) {
    cat(sprintf("✓ Word 报告已生成: %s\n", report_file))
  } else {
    cat("⚠ 报告生成失败\n")
  }
} else {
  cat("⚠ 报告生成脚本不存在，跳过报告生成\n")
}

if (!opt$quiet) {
  cat(sprintf("- 报告文件: %s\n", normalizePath(report_file, mustWork = FALSE)))
  cat("\n========================================\n")
  cat("✓ 分析完成！\n")
  cat("========================================\n")
}

# 写入统计信息和关闭日志
cat("\n----------------------------------------\n")
cat("统计信息\n")
cat("----------------------------------------\n")
cat(sprintf("样本总数: %d\n", nrow(d)))
cat(sprintf("  - %s: %d\n", comparison[1], sum(d$group == comparison[1])))
cat(sprintf("  - %s: %d\n", comparison[2], sum(d$group == comparison[2])))
cat(sprintf("目标基因: %d\n", length(genes)))
cat(sprintf("模型 C-index: %.4f\n", fit$stats["C"]))
cat(sprintf("AUC: %.4f\n", auc_val))
cat(sprintf("Hosmer-Lemeshow p: %s\n", hl_pval))

cat("\n----------------------------------------\n")
cat("输出文件\n")
cat("----------------------------------------\n")
for (f in output_files) {
  cat(sprintf("  - %s\n", basename(f)))
}

cat("\n----------------------------------------\n")
cat("追溯信息\n")
cat("----------------------------------------\n")
cat(sprintf("日志文件: %s.%s.log\n", SCRIPT_NAME, TIMESTAMP))
cat(sprintf("配置文件: %s\n", basename(config_file)))
cat(sprintf("报告文件: %s\n", basename(report_file)))

cat("\n================================================================================\n")
cat("分析成功完成\n")
cat("================================================================================\n")

# 关闭日志
close_logger(logger)
