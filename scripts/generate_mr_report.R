#!/usr/bin/env Rscript

# =============================================================================
# MR Analysis Report Generator - 优化版
# 参考真实报告样式：具体数据来源、软件版本、结果总结、图例说明
# =============================================================================

library(optparse)
library(officer)
library(flextable)
library(ggplot2)
library(dplyr)
library(magrittr)

# ============================================================================
# 字体与段落样式辅助函数
# ============================================================================
cn_font <- "SimSun"
en_font <- "Times New Roman"

# 逐字符混合字体 fpar：中文用 SimSun，英文用 Times New Roman
build_mixed_font_fpar <- function(text,
                                   align       = "left",
                                   font_size   = 10.5,
                                   bold        = FALSE,
                                   text_ind    = NULL,
                                   padding_top = 0,
                                   padding_bot = 0,
                                   word_style  = NULL) {
  if (!is.null(text_ind) && is.numeric(text_ind) && text_ind > 0 && nzchar(trimws(text))) {
    text <- paste0("\u3000\u3000", text)
  }
  para_args <- list(
    text.align     = align,
    padding.top    = padding_top,
    padding.bottom = padding_bot
  )
  if (!is.null(word_style) && is.character(word_style) && nzchar(word_style)) {
    para_args$word_style <- word_style
  }
  para_props <- do.call(officer::fp_par, para_args)

  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) {
    return(officer::fpar("", fp_p = para_props))
  }

  is_cn <- grepl("[一-龥]", chars)
  runs  <- list()
  start <- 1L

  for (index in seq_along(chars)) {
    is_break <- index == length(chars) || is_cn[[index + 1L]] != is_cn[[index]]
    if (is_break) {
      segment     <- paste(chars[start:index], collapse = "")
      font_family <- if (is_cn[[start]]) cn_font else en_font
      runs[[length(runs) + 1L]] <- officer::ftext(
        segment,
        prop = officer::fp_text(
          font.family = font_family,
          font.size   = font_size,
          bold        = bold
        )
      )
      start <- index + 1L
    }
  }

  do.call(officer::fpar, c(runs, list(fp_p = para_props)))
}

# 正文段落：首行缩进 2 字符（420 twips），混排字体
body_add_mixed_para <- function(doc, text, font_size = 10.5, bold = FALSE) {
  officer::body_add_fpar(
    doc,
    value = build_mixed_font_fpar(
      text,
      align     = "left",
      font_size = font_size,
      bold      = bold,
      text_ind  = 420  # 首行缩进 2 中文字符
    )
  )
}

# 章节标题：无缩进，混排字体，加粗
body_add_mixed_heading <- function(doc, text, font_size = 14, bold = TRUE) {
  officer::body_add_fpar(
    doc,
    value = build_mixed_font_fpar(
      text,
      align     = "left",
      font_size = font_size,
      bold      = bold,
      text_ind  = NULL  # 标题不缩进
    )
  )
}

# 图注段落：居中，无缩进，混排字体
body_add_mixed_figure_caption <- function(doc, text, font_size = 10.5, bold = FALSE) {
  officer::body_add_fpar(
    doc,
    value = build_mixed_font_fpar(
      text,
      align     = "center",
      font_size = font_size,
      bold      = bold,
      text_ind  = NULL  # 图注不缩进
    )
  )
}

# 空段落
body_add_empty <- function(doc) {
  officer::body_add_fpar(
    doc,
    value = officer::fpar("")
  )
}

# Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Result directory containing tables and figures"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output report file (.docx)"),
  make_option(c("-t", "--type"), type = "character", default = "eqtl",
              help = "MR analysis type: eqtl, pqtl, traditional, eqtl_pqtl, coloc"),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run")
)

opt <- parse_args(OptionParser(option_list = option_list))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- args[grep(file_arg, args)]
  if (length(hit) == 0) return(getwd())
  dirname(normalizePath(sub(file_arg, "", hit[[1]]), mustWork = FALSE))
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
report_template <- file.path(project_root, "templates", "mr_report_template.docx")

# Default timestamp
if (is.null(opt$timestamp)) {
  opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

# Default output
if (is.null(opt$output)) {
  opt$output <- paste0("MR_Report_", opt$type, "_", opt$timestamp, ".docx")
}

cat("================================================================================\n")
cat("  MR分析报告生成器（优化版）\n")
cat("  类型:", opt$type, "\n")
cat("  时间戳:", opt$timestamp, "\n")
cat("================================================================================\n\n")

# Set analysis type first (needed for column name detection)
analysis_type <- opt$type

# Set result directory from input
result_dir <- opt$input

# Set paths
tables_dir <- file.path(result_dir, "tables")
figures_dir <- file.path(result_dir, "figures")

# ============================================================================
# 读取数据
# ============================================================================
integrated_results <- NULL
is_traditional <- FALSE

csv_path <- file.path(result_dir, "00.Complete_MR_Results.csv")
if (file.exists(csv_path)) {
  integrated_results <- read.csv(csv_path, stringsAsFactors = FALSE, nrows = 1)
  if ("gene_symbol" %in% colnames(integrated_results)) {
    first_gene <- as.character(integrated_results$gene_symbol[1])
    if (grepl("^(ukb|finn)", first_gene)) {
      is_traditional <- TRUE
      cat("Detected as Traditional MR\n")
    }
  }
  integrated_results <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("Loaded results from 00.Complete_MR_Results.csv\n")
} else if (file.exists(file.path(tables_dir, "00.Complete_MR_Results.csv"))) {
  integrated_results <- read.csv(file.path(tables_dir, "00.Complete_MR_Results.csv"), stringsAsFactors = FALSE)
  cat("Loaded results from tables/00.Complete_MR_Results.csv\n")
}

# Read heterogeneity and pleiotropy from merged format files
het_results <- NULL
pleio_results <- NULL
het_file <- file.path(result_dir, "02.Heterogeneity_All_Genes.csv")
if (file.exists(het_file)) {
  het_results <- read.csv(het_file, stringsAsFactors = FALSE)
  cat("Loaded heterogeneity from:", het_file, "\n")
} else if (file.exists(file.path(tables_dir, "02.Heterogeneity_All_Genes.csv"))) {
  het_results <- read.csv(file.path(tables_dir, "02.Heterogeneity_All_Genes.csv"), stringsAsFactors = FALSE)
  cat("Loaded heterogeneity from tables/\n")
}

pleio_file <- file.path(result_dir, "03.Pleiotropy_All_Genes.csv")
if (file.exists(pleio_file)) {
  pleio_results <- read.csv(pleio_file, stringsAsFactors = FALSE)
  cat("Loaded pleiotropy from:", pleio_file, "\n")
} else if (file.exists(file.path(tables_dir, "03.Pleiotropy_All_Genes.csv"))) {
  pleio_results <- read.csv(file.path(tables_dir, "03.Pleiotropy_All_Genes.csv"), stringsAsFactors = FALSE)
  cat("Loaded pleiotropy from tables/\n")
}

# 读取01文件用于统计SNP数
mr01_results <- NULL
mr01_file <- file.path(result_dir, "01.MR_Results_All_Genes.csv")
if (file.exists(mr01_file)) {
  mr01_results <- read.csv(mr01_file, stringsAsFactors = FALSE)
  cat("Loaded MR results from:", mr01_file, "\n")
} else if (file.exists(file.path(tables_dir, "01.MR_Results_All_Genes.csv"))) {
  mr01_results <- read.csv(file.path(tables_dir, "01.MR_Results_All_Genes.csv"), stringsAsFactors = FALSE)
  cat("Loaded MR results from tables/01.MR_Results_All_Genes.csv\n")
}

# 读取三重筛选结果
three_filter_results <- NULL
tf_file <- file.path(result_dir, "06.Three_Filter_Results.csv")
if (file.exists(tf_file)) {
  three_filter_results <- read.csv(tf_file, stringsAsFactors = FALSE)
  cat("Loaded three-filter results from:", tf_file, "\n")
}

# 读取因果靶点结果
causal_results <- NULL
causal_file <- file.path(result_dir, "07.MR_res_gene.csv")
if (file.exists(causal_file)) {
  causal_results <- read.csv(causal_file, stringsAsFactors = FALSE)
  cat("Loaded causal targets from:", causal_file, "\n")
}

# ============================================================================
# 表格口径基因集：优先 07（最终因果），回退 06（三重筛选通过）
# ============================================================================
normalize_gene_ids <- function(values) {
  normalized <- toupper(trimws(as.character(values)))
  normalized[nzchar(normalized)]
}

final_table_genes <- character()
final_table_gene_source <- "all_table_rows"

if (!is.null(causal_results) && nrow(causal_results) > 0) {
  causal_id_cols <- c("gene_symbol", "exposure", "ensembl_id")[c("gene_symbol", "exposure", "ensembl_id") %in% colnames(causal_results)]
  if (length(causal_id_cols) > 0) {
    final_table_genes <- unique(unlist(lapply(causal_id_cols, function(col_name) {
      normalize_gene_ids(causal_results[[col_name]])
    }), use.names = FALSE))
    if (length(final_table_genes) > 0) final_table_gene_source <- "07.MR_res_gene.csv"
  }
}

if (length(final_table_genes) == 0 && !is.null(three_filter_results) && nrow(three_filter_results) > 0) {
  tf_cols <- colnames(three_filter_results)
  exposure_col <- tf_cols[tolower(tf_cols) == "exposure"][1]
  pass_col <- tf_cols[tolower(tf_cols) == "three_filter_pass"][1]
  if (!is.na(exposure_col) && nzchar(exposure_col)) {
    tf_candidates <- three_filter_results[[exposure_col]]
    if (!is.na(pass_col) && nzchar(pass_col)) {
      pass_vals <- toupper(trimws(as.character(three_filter_results[[pass_col]])))
      tf_candidates <- tf_candidates[pass_vals %in% c("YES", "TRUE", "1")]
    }
    final_table_genes <- unique(normalize_gene_ids(tf_candidates))
    if (length(final_table_genes) > 0) {
      final_table_gene_source <- "06.Three_Filter_Results.csv"
    }
  }
}

filter_by_final_genes <- function(df, candidate_cols) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  if (length(final_table_genes) == 0) {
    return(df[0, , drop = FALSE])
  }
  hit_cols <- candidate_cols[candidate_cols %in% colnames(df)]
  if (length(hit_cols) == 0) return(df[0, , drop = FALSE])
  keep <- rep(FALSE, nrow(df))
  for (col_name in hit_cols) {
    keep <- keep | (normalize_gene_ids(df[[col_name]]) %in% final_table_genes)
  }
  df[keep, , drop = FALSE]
}

mr01_results_for_table <- filter_by_final_genes(mr01_results, c("Exposure", "gene_symbol", "ensembl_id", "protein_id", "exposure"))
het_results_for_table <- filter_by_final_genes(het_results, c("Exposure", "gene_symbol", "ensembl_id", "protein_id", "exposure"))

if (length(final_table_genes) > 0) {
  cat(sprintf("Table genes source: %s (%d genes)\n", final_table_gene_source, length(final_table_genes)))
  if (!is.null(mr01_results) && nrow(mr01_results) > 0) {
    cat(sprintf("Table1 rows after filtering (01): %d / %d\n", nrow(mr01_results_for_table), nrow(mr01_results)))
  }
  if (!is.null(het_results) && nrow(het_results) > 0) {
    cat(sprintf("Table2 rows after filtering (02): %d / %d\n", nrow(het_results_for_table), nrow(het_results)))
  }
} else {
  cat("Table genes source: strict 07/06 filter (no genes found, Table1/Table2 will be empty)\n")
}

# 读取outcome名称
outcome_name <- "疾病"
if (!is.null(integrated_results) && nrow(integrated_results) > 0) {
  if ("outcome" %in% colnames(integrated_results)) {
    outcome_name <- as.character(integrated_results$outcome[1])
    outcome_name <- gsub("\\s*\\|\\|.*", "", outcome_name)
  }
}

# ============================================================================
# 读取分析参数配置
# ============================================================================
config_file <- file.path(result_dir, "analysis_config.ini")
params <- list(
  pval = "5×10⁻⁸",
  kb = "10000",
  r2 = "0.001",
  fstat = "10",
  eaf = "1%",
  minsnps = "3",
  exposure_type = "both",
  outcome_id = "疾病GWAS"
)

if (file.exists(config_file)) {
  config_lines <- readLines(config_file)
  for (line in config_lines) {
    if (grepl("^pval_threshold", line)) {
      val <- gsub(".*= ", "", line)
      if (val == "5e-8") params$pval <- "5×10⁻⁸"
      else if (val == "1e-6") params$pval <- "1×10⁻⁶"
      else if (val == "1e-5") params$pval <- "1×10⁻⁵"
      else if (val == "5e-6") params$pval <- "5×10⁻⁶"
      else params$pval <- val
    }
    if (grepl("^clump_kb", line)) params$kb <- gsub(".*= ", "", line)
    if (grepl("^clump_r2", line)) params$r2 <- gsub(".*= ", "", line)
    if (grepl("^fstat_threshold", line)) params$fstat <- gsub(".*= ", "", line)
    if (grepl("^eaf_threshold", line)) {
      val <- as.numeric(gsub(".*= ", "", line))
      params$eaf <- paste0(val * 100, "%")
    }
    if (grepl("^min_snps", line)) params$minsnps <- gsub(".*= ", "", line)
    if (grepl("^exposure_type", line)) params$exposure_type <- gsub(".*= ", "", line)
    if (grepl("^outcome_gwas", line)) params$outcome_id <- gsub(".*= ", "", line)
  }
  cat("Loaded config parameters from:", config_file, "\n")
} else {
  cat("No config file found, using default parameters\n")
}

# ============================================================================
# 统计结果数据（用于生成结果摘要）
# ============================================================================
# 统计基因数量
n_total_genes <- 0
n_significant_genes <- 0
significant_genes <- character()
n_total_snps <- 0
n_ivw_sig <- 0
n_het_pass <- 0
n_pleio_pass <- 0
n_steiger_pass <- 0
n_causal <- 0
causal_gene_list <- character()

if (!is.null(integrated_results) && nrow(integrated_results) > 0) {
  n_total_genes <- nrow(integrated_results)
  # 筛选显著的基因（IVW p < 0.05）
  if ("ivw_pval" %in% colnames(integrated_results)) {
    pval_vec <- suppressWarnings(as.numeric(integrated_results$ivw_pval))
    sig_idx <- which(!is.na(pval_vec) & pval_vec < 0.05)
    n_significant_genes <- length(sig_idx)
    if (n_significant_genes > 0) {
      if ("gene_symbol" %in% colnames(integrated_results)) {
        significant_genes <- as.character(integrated_results$gene_symbol[sig_idx])
      } else if ("exposure" %in% colnames(integrated_results)) {
        significant_genes <- as.character(integrated_results$exposure[sig_idx])
      }
    }
  }
  # 总SNP数
  if ("n_snps" %in% colnames(integrated_results)) {
    n_total_snps <- sum(as.numeric(integrated_results$n_snps), na.rm = TRUE)
  }
}

# 统计IVW显著的
if (!is.null(mr01_results) && nrow(mr01_results) > 0) {
  if (analysis_type == "traditional") {
    ivw_rows <- mr01_results[mr01_results$Method == "Inverse variance weighted", ]
    if ("P.value" %in% colnames(ivw_rows)) {
      ivw_sig_p <- suppressWarnings(as.numeric(ivw_rows$P.value))
      n_ivw_sig <- sum(!is.na(ivw_sig_p) & ivw_sig_p < 0.05, na.rm = TRUE)
    }
  } else {
    ivw_rows <- mr01_results[mr01_results$method == "Inverse variance weighted", ]
    if ("pval" %in% colnames(ivw_rows)) {
      ivw_sig_p <- suppressWarnings(as.numeric(ivw_rows$pval))
      n_ivw_sig <- sum(!is.na(ivw_sig_p) & ivw_sig_p < 0.05, na.rm = TRUE)
    }
  }
}

# 统计异质性通过的（Q_p > 0.05）
if (!is.null(het_results) && nrow(het_results) > 0) {
  q_col <- if ("Q_P_value" %in% colnames(het_results)) "Q_P_value" else "Q_P_value"
  q_pvals <- suppressWarnings(as.numeric(het_results[[q_col]]))
  # 只统计IVW行的异质性
  het_ivw <- het_results[grep("IVW|Inverse", het_results$Method, ignore.case = TRUE), ]
  if (nrow(het_ivw) > 0) {
    q_pvals_ivw <- suppressWarnings(as.numeric(het_ivw[[q_col]]))
    n_het_pass <- sum(!is.na(q_pvals_ivw) & q_pvals_ivw > 0.05, na.rm = TRUE)
  }
}

# 统计多效性通过的（Egger_intercept_p > 0.05）
if (!is.null(pleio_results) && nrow(pleio_results) > 0) {
  ep_col <- if ("Egger_intercept_pval" %in% colnames(pleio_results)) "Egger_intercept_pval" else "Egger_intercept_pval"
  ep_pvals <- suppressWarnings(as.numeric(pleio_results[[ep_col]]))
  n_pleio_pass <- sum(!is.na(ep_pvals) & ep_pvals > 0.05, na.rm = TRUE)
}

# 因果靶点数量
n_causal <- if (!is.null(causal_results) && nrow(causal_results) > 0) nrow(causal_results) else 0
if (n_causal > 0) {
  if ("gene_symbol" %in% colnames(causal_results)) {
    causal_gene_list <- paste(head(as.character(causal_results$gene_symbol), 3), collapse = "、")
    if (n_causal > 3) causal_gene_list <- paste0(causal_gene_list, "等")
  } else if ("exposure" %in% colnames(causal_results)) {
    causal_gene_list <- paste(head(as.character(causal_results$exposure), 3), collapse = "、")
    if (n_causal > 3) causal_gene_list <- paste0(causal_gene_list, "等")
  }
}

cat("Statistics: Total genes:", n_total_genes, ", Significant:", n_significant_genes, "\n")

# 补充统计（用于文段描述）
or_range_text <- ""
direction_text <- ""
or_vals <- numeric()
if (!is.null(integrated_results) && nrow(integrated_results) > 0) {
  # eQTL/pQTL: use ivw_or
  if ("ivw_or" %in% colnames(integrated_results)) {
    or_vals <- suppressWarnings(as.numeric(integrated_results$ivw_or))
  }
  # Traditional: use ivw_or if available (some versions have it)
  if (length(or_vals) == 0 && "or" %in% colnames(integrated_results)) {
    or_vals <- suppressWarnings(as.numeric(integrated_results$or))
  }
  or_vals <- or_vals[!is.na(or_vals)]
}
# Traditional: also try to get OR from per-exposure mr_results files
if (length(or_vals) == 0 && exists("tables_dir")) {
  mr_files <- list.files(tables_dir, pattern = "_mr_results\\.csv$", full.names = TRUE, recursive = TRUE)
  for (f in mr_files) {
    mr_dt <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(mr_dt) && nrow(mr_dt) > 0 && "or" %in% colnames(mr_dt)) {
      ivw_row <- mr_dt[grep("Inverse variance weighted", mr_dt$method, ignore.case = TRUE), ]
      if (nrow(ivw_row) > 0) {
        or_v <- suppressWarnings(as.numeric(ivw_row$or[1]))
        if (!is.na(or_v)) or_vals <- c(or_vals, or_v)
      }
    }
  }
}
if (length(or_vals) > 0) {
  or_min <- sprintf("%.3f", min(or_vals))
  or_max <- sprintf("%.3f", max(or_vals))
  if (or_min == or_max) {
    or_range_text <- paste0("所有暴露的IVW效应值（OR）均为", or_min, "。")
  } else {
    or_range_text <- paste0("所有暴露的IVW效应值（OR）范围为", or_min, "至", or_max, "。")
  }
}

# ============================================================================
# 创建Word文档（优先使用项目模板）
# ============================================================================
doc <- if (file.exists(report_template)) {
  cat("Using report template:", report_template, "\n")
  read_docx(path = report_template)
} else {
  read_docx()
}

# ========== 标题 ==========
doc <- doc %>%
  body_add_mixed_heading("孟德尔随机化分析报告", font_size = 16, bold = TRUE) %>%
  body_add_empty()

# ========== 1. 方法 ==========
doc <- doc %>%
  body_add_mixed_heading("1. 方法", font_size = 14, bold = TRUE) %>%
  body_add_empty()

# Read successful genes list
successful_genes <- NULL
success_file <- if (analysis_type == "pqtl") {
  file.path(result_dir, "09.pQTL_Successful_Genes.csv")
} else if (analysis_type == "traditional") {
  file.path(result_dir, "09.Traditional_Successful_Exposures.csv")
} else {
  file.path(result_dir, "09.eQTL_Successful_Genes.csv")
}
if (file.exists(success_file)) {
  successful_genes <- read.csv(success_file, stringsAsFactors = FALSE)
  cat("Loaded successful genes from:", success_file, "\n")
}

# 方法精简版（核心包 + 关键参数）
pkg_tsmr <- tryCatch(as.character(utils::packageVersion("TwoSampleMR")), error = function(e) "unknown")
pkg_mrpresso <- tryCatch(as.character(utils::packageVersion("MRPRESSO")), error = function(e) "unknown")
r_ver <- as.character(getRversion())
exposure_term <- if (analysis_type == "traditional") "暴露因素" else if (analysis_type == "pqtl") "pQTL特征" else "候选基因"
method_text <- paste0(
  "本研究采用R（v", r_ver, "，https://www.r-project.org）与TwoSampleMR（v", pkg_tsmr, "，https://mrcieu.github.io/TwoSampleMR/）开展双样本MR分析，以",
  exposure_term, "为暴露、疾病为结局。工具变量筛选参数为P<", params$pval, "、r2<", params$r2, "、kb=", params$kb, "、F>", params$fstat, "、MAF>", params$eaf,
  "，并进行harmonise_data协调。主分析采用IVW（P<0.05），结合Cochran Q、MR-Egger及MR-PRESSO（v", pkg_mrpresso, "，https://github.com/rondolab/MR-PRESSO）进行稳健性评估，并辅以leave-one-out和Steiger方向检验。"
)
doc <- doc %>%
  body_add_mixed_heading("1.1 方法概述", font_size = 12, bold = TRUE) %>%
  body_add_mixed_para(method_text) %>%
  body_add_empty()

# 旧版长方法段保留但不执行
if (FALSE) {
if (analysis_type == "eqtl") {
  doc <- doc %>%
    body_add_mixed_heading("1.1 GWAS数据获取", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过IEU OpenGWAS数据库（https://gwas.mrcieu.ac.uk/）获取疾病相关的GWAS数据。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.2 eQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过双样本MR分析来探讨候选基因与疾病风险之间的因果关系，将单核苷酸多态性（SNPs）定义为工具变量（IVs）。从血液eQTL的数据获取候选基因的eQTL数据（https://www.eQTLGen.com），作为暴露因素，并通过IEU OpenGWAS数据库获得疾病的GWAS数据作为结局因素。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.3 工具变量筛选", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para(paste0("首先，进行工具变量的筛选。通过R包TwoSampleMR（https://mrcieu.github.io/TwoSampleMR/）读取暴露因素并确定与暴露因素独立相关的SNP，将它们作为评估暴露因素和结局之间因果关系的工具变量（IVs）。筛选指标：包括P<", params$pval, "，寻找与暴露因素显著相关的工具变量（SNP）；使用clump=TRUE，r2=", params$r2, "，kb=", params$kb, "这三个指标去掉存在连锁不平衡（LD）的工具变量（SNP）。去掉与结局显著相关的工具变量（剔除F-statistics<", params$fstat, "的SNPs），并将暴露因素-工具变量-结局进行匹配。SNP数目尽量大于", params$minsnps, "个；保留最小等位基因频率MAF>", params$eaf, "位点；排除具有中等位基因频率的回文SNP；使用Harmonise函数协调exposure和outcome位点。")) %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.4 孟德尔随机化分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("为探讨候选基因与疾病风险之间的因果关系，使用TwoSampleMR软件包（https://mrcieu.github.io/TwoSampleMR/）进行MR分析（五种算法：Inverse variance weighted、MR Egger、Weighted median、Simple mode、Weighted mode），并基于IVW方法的结果来评估候选基因与疾病风险之间的关系（P<0.05）。我们还采用Cochran的Q统计量来检验异质性，其中p>0.05表明不存在异质性。使用MR-Egger回归和MR-PRESSO分析评估潜在的水平多效性，p>0.05表明不存在水平多效性。还进行留一法检验来判断结果具有可靠性。最后，为了验证正向分析的结果不受反向因果效应的干扰，基于IVW方法，用R包TwoSampleMR进行Steiger方向性分析（correct_causal_direction=1/TRUE，steiger_test_adj<0.05才能说明通过了反向的检验，如果计算出NA，则视为不通过。）") %>%
    body_add_empty()

} else if (analysis_type == "pqtl") {
  doc <- doc %>%
    body_add_mixed_heading("1.1 GWAS数据获取", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过IEU OpenGWAS数据库（https://gwas.mrcieu.ac.uk/）获取疾病相关的GWAS数据。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.2 pQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("从血浆蛋白pQTL的数据（https://www.pqtlgen.org/）中获取候选基因的pQTL数据作为暴露因素，通过IEU OpenGWAS数据库获得疾病的GWAS数据作为结局因素。工具变量筛选步骤同eQTL分析。最终整理两部分的MR结果取并集作为因果基因进行后续分析。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.3 工具变量筛选", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para(paste0("筛选条件：包括P<", params$pval, "，寻找与暴露因素显著相关的工具变量（SNP）；使用clump=TRUE，r2=", params$r2, "，kb=", params$kb, "去掉存在连锁不平衡（LD）的工具变量。去掉与结局显著相关的工具变量（剔除F-statistics<", params$fstat, "的SNPs）。SNP数目尽量大于", params$minsnps, "个；保留最小等位基因频率MAF>", params$eaf, "位点；使用Harmonise函数协调exposure和outcome位点。")) %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.4 孟德尔随机化分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("使用TwoSampleMR软件包（https://mrcieu.github.io/TwoSampleMR/）进行MR分析（五种算法：Inverse variance weighted、MR Egger、Weighted median、Simple mode、Weighted mode），并基于IVW方法的结果来评估候选基因与疾病风险之间的关系（P<0.05）。采用Cochran的Q统计量来检验异质性，p>0.05表明不存在异质性。使用MR-Egger回归和MR-PRESSO软件包（https://github.com/rondolab/MR-PRESSO）分析评估潜在的水平多效性，p>0.05表明不存在水平多效性。还进行留一法检验来判断结果的可靠性。最后进行Steiger方向性分析验证因果方向。") %>%
    body_add_empty()

} else if (analysis_type == "traditional") {
  doc <- doc %>%
    body_add_mixed_heading("1.1 GWAS数据获取", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("使用IEU OpenGWAS数据库（https://gwas.mrcieu.ac.uk/）检索，获取暴露因素和结局的GWAS数据。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.2 工具变量筛选", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para(paste0("通过R包TwoSampleMR（https://mrcieu.github.io/TwoSampleMR/）中的extract_instruments进行暴露因素的读取和工具变量的筛选。筛选条件：p<", params$pval, "，寻找与暴露因素显著相关的工具变量（SNP），使用clump=TRUE，r2=", params$r2, "，kb=", params$kb, "，去掉存在连锁不平衡（LD）的工具变量（SNP）。基于结局GWAS数据，以及前面筛选后的工具变量，去掉与结局显著相关的工具变量，通过R包TwoSampleMR函数harmonise_data统一效应等位与效应量，并将暴露因素-工具变量-结局进行匹配。对所有工具变量进行F统计量计算，全部F统计量大于", params$fstat, "，不存在弱工具变量。")) %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.3 孟德尔随机化分析暴露因素与结局的关系", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("使用R包TwoSampleMR的mr函数结合五种算法来进行MR分析（MR Egger，Weighted median，Inverse variance weighted，Simple mode，Weighted mode）。结果主要参考Inverse variance weighted（IVW），并参考后续的异质性分析结果的Cochran Q的评估结果，存在异质性则使用随机效应IVW，不存在异质性则使用固定效应IVW。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.4 暴露因素与结局的相关性分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("绘制散点图，判断暴露因素与结局的相关性，横坐标是SNP对暴露的效应，纵坐标是SNP对结局的效应，彩色的线表示的是MR不同算法的拟合结果。当截距为0时，表明工具变量对暴露因素和结局的效应相当。当截距不为0时，则暗示可能存在混杂因素。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.5 暴露因素对于结局的诊断效能", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("绘制森林图，对各暴露因素的SNP位点对结局的诊断效能进行判断。大于0则表明该SNP位点与结局是正向的关系，小于0表明该SNP位点是负向的关系。使用IVW的方法对SNP的整体诊断效能进行评价（红点和红线）。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.6 随机性判断", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("绘制漏斗图，根据分布情况，判断分析是否符合孟德尔第二定律。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.7 敏感性分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("（1）异质性检验（Heterogeneity test）：通过评估使用的不同工具变量得到结果的差异性，来反映选择的工具变量的有效性。具体而言，使用R package TwoSampleMR的mr_heterogeneity函数，计算Cochran's Q统计量。Q=Σwj(βj-β)²，其中βj是第j个IV得到的系数估计值，wj是对应的权重，β是利用IVW或者MR-Egger合并得到的汇总估计值（pooled estimate）。") %>%
    body_add_mixed_para("（2）多效性检测（Pleiotropy test）：孟德尔随机化的假设的其中之一是工具变量必须通过暴露因素来影响结果。如果工具变量可以不通过暴露因素直接影响结果，即检验结果存在水平多效性，就违反了孟德尔随机化的思想。使用mr_pleiotropy_test函数检验评价是否存在水平多效性，不存在水平多效性表明结果可靠。p值小于0.05，表明存在水平多效性。同时使用MR-PRESSO软件包（https://github.com/rondolab/MR-PRESSO）进行水平多效性检验。") %>%
    body_add_mixed_para("（3）逐个剔除检验（Leave-one-out sensitivity test）：通过leave-one-out函数进行留一分析以评估MR结果是否依赖于特定的SNP，逐步剔除每个SNP，计算剩余SNP的meta效应，观察剔除每个SNP后结果是否发生变化，如果剔除了某一个SNP后，结果改变很大，说明存在某一个SNP对结果影响很大。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.8 Steiger方向性分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过R包TwoSampleMR进行Steiger方向性分析，以判断方向是否正确，判断为TRUE，表明方向正确。") %>%
    body_add_empty()

} else if (analysis_type == "eqtl_pqtl") {
  doc <- doc %>%
    body_add_mixed_heading("1.1 GWAS数据获取", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过IEU OpenGWAS数据库（https://gwas.mrcieu.ac.uk/）获取疾病相关的GWAS数据。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.2 eQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过双样本MR分析来探讨候选基因与疾病风险之间的因果关系，将单核苷酸多态性（SNPs）定义为工具变量（IVs）。从血液eQTL的数据获取候选基因的eQTL数据（https://www.eQTLGen.com），作为暴露因素，并通过IEU OpenGWAS数据库获得疾病的GWAS数据作为结局因素。工具变量筛选、MR分析步骤同前所述。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.3 pQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("从血浆蛋白pQTL的数据（https://www.pqtlgen.org/）中获取候选基因的pQTL数据作为暴露因素，通过IEU OpenGWAS数据库获得疾病的GWAS数据作为结局因素。工具变量筛选、MR分析步骤同前所述。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.4 MR分析整合", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("整理两部分的MR结果取并集作为因果基因进行后续分析。") %>%
    body_add_empty()

} else if (analysis_type == "coloc") {
  doc <- doc %>%
    body_add_mixed_heading("1.1 GWAS数据获取", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过IEU OpenGWAS数据库（https://gwas.mrcieu.ac.uk/）获取疾病相关的GWAS数据。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.2 eQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("通过双样本MR分析来探讨候选基因与疾病风险之间的因果关系。MR分析步骤同前所述（详见1.2-1.5）。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.3 pQTL的MR分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("从血浆蛋白pQTL的数据（https://www.pqtlgen.org/）中获取候选基因的pQTL数据作为暴露因素进行分析。") %>%
    body_add_empty() %>%
    body_add_mixed_heading("1.4 贝叶斯共定位分析", font_size = 12, bold = TRUE) %>%
    body_add_mixed_para("为了进一步加强因果推断的可靠性，使用R包coloc（https://chr1swallace.github.io/coloc/index.html）对暴露因素（通过Steiger方向性分析的候选暴露因素）与疾病进行贝叶斯共定位分析。共定位分析共包含五种假设：") %>%
    body_add_mixed_para("H0：表型1（GWAS）和表型2（eQTL/pQTL）与某个基因组区域的所有SNP位点无显著相关；") %>%
    body_add_mixed_para("H1/H2：表型1（GWAS）或表型2（eQTL/pQTL）与某个基因组区域的SNP位点显著相关；") %>%
    body_add_mixed_para("H3：表型1（GWAS）和表型2（eQTL/pQTL）与某个基因组区域的SNP位点显著相关，但由不同的因果变异位点驱动；") %>%
    body_add_mixed_para("H4：表型1（GWAS）和表型2（eQTL/pQTL）与某个基因组区域的SNP位点显著相关，且由同一个因果变异位点驱动。") %>%
    body_add_mixed_para("PP.H4.abf值大于0.6认为共定位阳性（60%概率由同一个连锁区间的SNP突变影响了两个表型，且由同一个因果变异位点驱动）。选择PP.H4>0.6的基因记作候选关键基因进行后续分析。") %>%
    body_add_empty()
}

}

# ========== 2. 结果 ==========
doc <- doc %>%
  body_add_mixed_heading("2. 结果", font_size = 14, bold = TRUE) %>%
  body_add_empty()

# ---- 2.1 MR分析结果 ----
doc <- doc %>% body_add_mixed_heading("2.1 孟德尔随机化分析结果", font_size = 12, bold = TRUE)

# 生成结果叙述文字
if (n_total_genes > 0) {
  if (analysis_type == "traditional") {
    # 传统MR叙述
    summary_text <- paste0(
      "按照方法所述，将筛选出的", n_total_genes, "个候选靶点作为暴露变量，疾病作为结局变量，",
      "进行两样本MR分析，以探讨候选靶点与疾病风险之间的因果关系。",
      "按方法进行数据过滤后，使用harmonise_data函数将暴露数据与结局数据进行合并，",
      "合并后共有", ifelse(n_total_snps > 0, n_total_snps, "若干"), "个SNP用于MR分析。"
    )
    if (n_ivw_sig > 0) {
      summary_text <- paste0(summary_text, "MR分析发现有", n_ivw_sig, "个暴露变量与",
        "疾病风险之间有因果关系（p<0.05）。")
      if (n_causal > 0 && nchar(causal_gene_list) > 0) {
        summary_text <- paste0(summary_text, "其中",
          ifelse(n_het_pass > 0, paste0(n_het_pass, "个暴露变量不存在异质性与水平多效性（p>0.05），"), ""),
          "因此将这", n_causal, "个靶点作为因果靶点",
          ifelse(nchar(causal_gene_list) > 0, paste0("（包括", causal_gene_list, "）"), ""),
          "。")
      }
    } else {
      summary_text <- paste0(summary_text, "MR分析未发现显著的因果关系（p<0.05）。")
      if (n_het_pass > 0) {
        summary_text <- paste0(summary_text, "其中", n_het_pass, "个暴露变量不存在异质性（Q p>0.05），",
          n_pleio_pass, "个暴露变量不存在水平多效性（p>0.05）。")
      }
    }
    doc <- doc %>% body_add_mixed_para(summary_text)
  } else {
    # eQTL/pQTL叙述
    summary_text <- paste0(
      "按照方法所述，对", n_total_genes, "个候选基因进行双样本孟德尔随机化分析，",
      "以探讨候选基因与疾病风险之间的因果关系。"
    )
    if (n_ivw_sig > 0) {
      summary_text <- paste0(summary_text, "其中有", n_ivw_sig, "个基因的IVW方法P值<0.05，",
        "表明这些基因与疾病风险存在显著的因果关联。")
    } else {
      summary_text <- paste0(summary_text, "本次分析未发现显著的因果关系（IVW方法P值均≥0.05）。")
    }
    if (nchar(or_range_text) > 0) {
      summary_text <- paste0(summary_text, or_range_text)
    }
    if (nchar(direction_text) > 0) {
      summary_text <- paste0(summary_text, direction_text)
    }
    if (n_causal > 0) {
      summary_text <- paste0(summary_text, "经过异质性检验、水平多效性检验和方向性分析后，",
        "共有", n_causal, "个基因通过三重筛选",
        ifelse(nchar(causal_gene_list) > 0, paste0("（包括", causal_gene_list, "）"), ""),
        "，作为因果靶点进行后续分析。")
    } else if (n_total_genes > 0 && n_ivw_sig == 0) {
      n_het_fail <- n_total_genes - n_het_pass
      n_pleio_fail <- n_total_genes - n_pleio_pass
      het_desc <- if (n_het_fail == n_total_genes) {
        "所有基因均存在显著异质性（Q p<0.05）"
      } else if (n_het_pass == n_total_genes) {
        "所有基因均不存在显著异质性（Q p>0.05）"
      } else {
        paste0(n_het_fail, "个基因存在显著异质性，", n_het_pass, "个基因不存在显著异质性（Q p>0.05）")
      }
      pleio_desc <- if (n_pleio_fail == n_total_genes) {
        "所有基因均存在水平多效性（p<0.05）"
      } else if (n_pleio_pass == n_total_genes) {
        "所有基因均不存在水平多效性（p>0.05）"
      } else {
        paste0(n_pleio_fail, "个基因存在水平多效性，", n_pleio_pass, "个基因不存在水平多效性（p>0.05）")
      }
      summary_text <- paste0(summary_text, "异质性检验显示，", het_desc, "；水平多效性检验显示，", pleio_desc, "。")
    }
    doc <- doc %>% body_add_mixed_para(summary_text)
  }
}

# ---- Table 1 ----
doc <- doc %>% body_add_empty()
doc <- doc %>% body_add_mixed_figure_caption("表1 不同暴露因素的孟德尔随机化分析结果")

if (is_traditional) {
  table1_list <- list()
  mr_files <- list.files(tables_dir, pattern = "_mr_results\\.csv$", full.names = TRUE, recursive = TRUE)
  cat("Found", length(mr_files), "mr_results files\n")

  for (f in mr_files) {
    mr_dt <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(mr_dt) == 0) next

    exp_name <- gsub("_mr_results\\.csv$", "", basename(f))
    exp_name <- gsub("_", " ", exp_name)

    outcome_name <- ""
    if ("outcome" %in% colnames(mr_dt)) {
      outcome_name <- as.character(mr_dt$outcome[1])
    } else if ("id.outcome" %in% colnames(mr_dt)) {
      outcome_name <- as.character(mr_dt$id.outcome[1])
    }

    nsnp_col <- if ("nsnp" %in% colnames(mr_dt)) "nsnp" else "n_snps"

    for (i in 1:nrow(mr_dt)) {
      method <- as.character(mr_dt$method[i])
      n_snps_val <- mr_dt[i, nsnp_col]
      pval <- mr_dt$pval[i]
      or_val <- mr_dt$or[i]

      pval_str <- ifelse(is.na(pval), "", tryCatch(sprintf("%.4f", as.numeric(pval)), error = function(e) ""))
      or_str <- ifelse(is.na(or_val), "", tryCatch(sprintf("%.3f", as.numeric(or_val)), error = function(e) ""))

      table1_list[[length(table1_list) + 1]] <- data.frame(
        Exposure = exp_name,
        Outcome = "Disease",
        Method = method,
        SNP = n_snps_val,
        p_value = pval_str,
        OR = or_str,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(table1_list) > 0) {
    table1_df <- do.call(rbind, table1_list)

    ft <- flextable(table1_df) %>%
      bold(part = "header", bold = TRUE) %>%
      fontsize(size = 9, part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      align(align = "center", part = "all") %>%
      align(align = "left", j = 1:2) %>%
      merge_v(j = 1) %>%
      border_remove() %>%
      border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
      border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
      width(j = 1:6, width = c(1.1, 0.8, 2.4, 0.5, 0.8, 0.7))

    doc <- doc %>% body_add_flextable(ft, align = "center")
    # 表1内容解释
    doc <- doc %>% body_add_empty()
    table1_ivw_sig <- if (nrow(table1_df) > 0) {
      ivw_rows <- table1_df[grepl("Inverse variance weighted|IVW", table1_df$Method, ignore.case = TRUE), , drop = FALSE]
      if (nrow(ivw_rows) > 0) {
        pvals <- suppressWarnings(as.numeric(ivw_rows$p_value))
        sum(!is.na(pvals) & pvals < 0.05, na.rm = TRUE)
      } else {
        0
      }
    } else {
      0
    }
    if (table1_ivw_sig > 0) {
      table1_desc <- paste0("表1结果显示，在", table1_ivw_sig, "个暴露因素中，IVW方法分析显示存在显著因果关联（p<0.05）。")
    } else {
      table1_desc <- paste0("表1结果显示，所有暴露因素的IVW方法分析均未发现显著因果关联（p≥0.05）。")
    }
    if (nchar(or_range_text) > 0) {
      table1_desc <- paste0(table1_desc, " ", or_range_text)
    }
    if (nchar(direction_text) > 0) {
      table1_desc <- paste0(table1_desc, " ", direction_text)
    }
    doc <- doc %>% body_add_mixed_para(table1_desc)
  }

} else if (!is.null(mr01_results_for_table) && nrow(mr01_results_for_table) > 0) {
  table1_list <- list()
  get_method_prefix <- function(method_name) {
    m <- tolower(trimws(as.character(method_name)))
    if (grepl("inverse variance weighted|\\bivw\\b", m)) return("ivw")
    if (grepl("egger", m)) return("egger")
    if (grepl("weighted median", m)) return("wm")
    if (grepl("weighted mode", m)) return("wmode")
    if (grepl("simple mode", m)) return("smode")
    NA_character_
  }
  for (i in 1:nrow(mr01_results_for_table)) {
    row <- mr01_results_for_table[i, , drop = FALSE]
    exp_name <- if ("Exposure" %in% colnames(row)) {
      as.character(row$Exposure[1])
    } else if ("gene_symbol" %in% colnames(row)) {
      as.character(row$gene_symbol[1])
    } else if ("ensembl_id" %in% colnames(row)) {
      as.character(row$ensembl_id[1])
    } else {
      as.character(i)
    }
    method_name <- if ("Method" %in% colnames(row)) as.character(row$Method[1]) else ""
    n_snps_val <- if ("n_SNPs" %in% colnames(row)) row$n_SNPs[1] else if ("n_snps" %in% colnames(row)) row$n_snps[1] else if ("nsnp" %in% colnames(row)) row$nsnp[1] else NA
    pval <- NA_real_
    or_val <- NA_real_
    if ("pval" %in% colnames(row)) pval <- suppressWarnings(as.numeric(row$pval[1]))
    if ("or" %in% colnames(row)) or_val <- suppressWarnings(as.numeric(row$or[1]))
    if (is.na(pval) || is.na(or_val)) {
      prefix <- get_method_prefix(method_name)
      if (!is.na(prefix)) {
        p_col <- paste0(prefix, "_pval")
        or_col <- paste0(prefix, "_or")
        if (is.na(pval) && p_col %in% colnames(row)) pval <- suppressWarnings(as.numeric(row[[p_col]][1]))
        if (is.na(or_val) && or_col %in% colnames(row)) or_val <- suppressWarnings(as.numeric(row[[or_col]][1]))
      }
    }
    table1_list[[length(table1_list) + 1]] <- data.frame(
      Exposure = exp_name,
      Outcome = "Disease",
      Method = method_name,
      SNP = n_snps_val,
      p_value = ifelse(is.na(pval), "", sprintf("%.4f", pval)),
      OR = ifelse(is.na(or_val), "", sprintf("%.3f", or_val)),
      stringsAsFactors = FALSE
    )
  }

  if (length(table1_list) > 0) {
    table1_df <- do.call(rbind, table1_list)

    ft <- flextable(table1_df) %>%
      bold(part = "header", bold = TRUE) %>%
      fontsize(size = 9, part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      align(align = "center", part = "all") %>%
      align(align = "left", j = 1:2) %>%
      merge_v(j = 1) %>%
      border_remove() %>%
      border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
      border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
      border(part = "body", border.bottom = officer::fp_border(width = 1.5, color = "black"),
             i = nrow(table1_df)) %>%
      width(j = 1:6, width = c(1.0, 0.7, 1.5, 0.5, 0.8, 0.7))

    doc <- doc %>% body_add_flextable(ft, align = "center")
    # 表1内容解释
    doc <- doc %>% body_add_empty()
    table1_ivw_sig <- if ("Method" %in% colnames(mr01_results_for_table)) {
      ivw_rows <- mr01_results_for_table[grepl("Inverse variance weighted|IVW", mr01_results_for_table$Method, ignore.case = TRUE), , drop = FALSE]
      if (nrow(ivw_rows) > 0) {
        if ("pval" %in% colnames(ivw_rows)) {
          sum(suppressWarnings(as.numeric(ivw_rows$pval)) < 0.05, na.rm = TRUE)
        } else if ("ivw_pval" %in% colnames(ivw_rows)) {
          sum(suppressWarnings(as.numeric(ivw_rows$ivw_pval)) < 0.05, na.rm = TRUE)
        } else {
          0
        }
      } else {
        0
      }
    } else {
      0
    }
    if (table1_ivw_sig > 0) {
      table1_desc <- paste0("表1结果显示，在", table1_ivw_sig, "个候选基因中，IVW方法分析显示存在显著因果关联（p<0.05）。")
    } else {
      table1_desc <- paste0("表1结果显示，所有候选基因的IVW方法分析均未发现显著因果关联（p≥0.05）。")
    }
    if (nchar(or_range_text) > 0) {
      table1_desc <- paste0(table1_desc, " ", or_range_text)
    }
    if (nchar(direction_text) > 0) {
      table1_desc <- paste0(table1_desc, " ", direction_text)
    }
    doc <- doc %>% body_add_mixed_para(table1_desc)
  }
}

doc <- doc %>% body_add_empty()

# ---- 2.2 敏感性分析结果 ----
doc <- doc %>% body_add_mixed_heading("2.2 敏感性分析结果", font_size = 12, bold = TRUE)
doc <- doc %>% body_add_empty()
doc <- doc %>% body_add_mixed_figure_caption("表2 暴露与疾病之间的敏感性分析")

# 从合并格式文件构建 Table 2
has_het_data <- !is.null(het_results_for_table) && nrow(het_results_for_table) > 0

if (has_het_data) {
  sens_list <- list()

  # 收集所有 exposure
  all_exposures <- unique(het_results_for_table$Exposure)

  for (exp in all_exposures) {
    exp_het <- if (has_het_data) het_results_for_table[het_results_for_table$Exposure == exp, ] else NULL

    # heterogeneity 文件的每一行直接用原始列
    if (!is.null(exp_het) && nrow(exp_het) > 0) {
      for (r in 1:nrow(exp_het)) {
        row <- exp_het[r, ]
        safe_num <- function(x) {
          v <- suppressWarnings(as.numeric(x))
          if (length(v) == 0) return(NA_real_)
          v[1]
        }
        sens_list[[length(sens_list) + 1]] <- data.frame(
          Exposure = exp,
          Outcome = "Disease",
          Method = as.character(row$Method[1]),
          Heterogeneity_Q = ifelse(is.na(safe_num(row$Heterogeneity_Q)), "",
                                   sprintf("%.4f", safe_num(row$Heterogeneity_Q))),
          Q_df = ifelse(is.na(safe_num(row$Q_df)), "",
                        as.integer(safe_num(row$Q_df))),
          Q_P_value = ifelse(is.na(safe_num(row$Q_P_value)), "",
                             sprintf("%.4f", safe_num(row$Q_P_value))),
          Egger_intercept_pval = ifelse(is.na(safe_num(row$Egger_intercept_pval)), "",
                                        sprintf("%.6f", safe_num(row$Egger_intercept_pval))),
          MR_PRESSO_global_p = ifelse(is.na(safe_num(row[["MR_PRESSO_global_p"]])), "",
                                      sprintf("%.4f", safe_num(row[["MR_PRESSO_global_p"]])))
        )
      }
    }
  }

  if (length(sens_list) > 0) {
    sens_df <- do.call(rbind, sens_list)

    ft <- flextable(sens_df) %>%
      bold(part = "header", bold = TRUE) %>%
      fontsize(size = 9, part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      align(align = "center", part = "all") %>%
      align(align = "left", j = c(1, 2, 3)) %>%
      merge_v(j = 1) %>%
      border_remove() %>%
      border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
      border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
      border(part = "body", border.bottom = officer::fp_border(width = 1.5, color = "black"),
             i = nrow(sens_df)) %>%
      width(j = 1:8, width = c(1.0, 0.7, 1.5, 0.8, 0.5, 0.8, 0.8, 0.8))

    doc <- doc %>% body_add_flextable(ft, align = "center")
    # 表2内容解释
    doc <- doc %>% body_add_empty()
    # n_het_pass 只算 IVW 行
    het_ivw <- het_results_for_table[grep("IVW|Inverse", het_results_for_table$Method, ignore.case = TRUE), ]
    ep_pvals_het <- suppressWarnings(as.numeric(het_ivw$Q_P_value))
    n_het_pass <- sum(!is.na(ep_pvals_het) & ep_pvals_het > 0.05, na.rm = TRUE)
    # n_pleio_pass 按 exposure 统计（避免 MR-Egger/IVW 双行重复计数）
    pleio_by_exposure <- lapply(split(het_results_for_table$Egger_intercept_pval, het_results_for_table$Exposure), function(values) {
      num <- suppressWarnings(as.numeric(values))
      num <- num[!is.na(num)]
      if (length(num) == 0) NA_real_ else num[[1]]
    })
    pleio_pvals <- unlist(pleio_by_exposure, use.names = FALSE)
    n_pleio_pass <- sum(!is.na(pleio_pvals) & pleio_pvals > 0.05, na.rm = TRUE)
    table2_total_genes <- length(unique(het_results_for_table$Exposure))
    n_het_fail <- table2_total_genes - n_het_pass
    n_pleio_fail <- table2_total_genes - n_pleio_pass
    het_desc2 <- if (n_het_fail == table2_total_genes) {
      "所有暴露因素均存在显著异质性（Q p<0.05）"
    } else if (n_het_pass == table2_total_genes) {
      "所有暴露因素均不存在显著异质性（Q p>0.05）"
    } else {
      paste0(n_het_fail, "个暴露因素存在显著异质性，", n_het_pass, "个不存在异质性（Q p>0.05）")
    }
    pleio_desc2 <- if (n_pleio_fail == table2_total_genes) {
      "所有暴露因素均存在水平多效性（p<0.05）"
    } else if (n_pleio_pass == table2_total_genes) {
      "所有暴露因素均不存在水平多效性（p>0.05）"
    } else {
      paste0(n_pleio_fail, "个暴露因素存在水平多效性，", n_pleio_pass, "个不存在水平多效性（p>0.05）")
    }
    doc <- doc %>% body_add_mixed_para(paste0("表2 结果显示，异质性检验，", het_desc2, "；水平多效性检验，", pleio_desc2, "。"))
  }
}

doc <- doc %>% body_add_empty()

# ---- 2.3 可视化结果 ----
doc <- doc %>% body_add_mixed_heading("2.3 可视化结果", font_size = 12, bold = TRUE)

doc <- doc %>% body_add_empty()

# 获取展示基因（显著基因优先，否则取p值最小的1个基因）
target_genes <- character()
if (!is.null(integrated_results) && nrow(integrated_results) > 0) {
  if ("ivw_pval" %in% colnames(integrated_results)) {
    sig_idx <- which(as.numeric(integrated_results$ivw_pval) < 0.05)
    if (length(sig_idx) > 0) {
      if ("gene_symbol" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$gene_symbol[sig_idx])
      } else if ("protein_id" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$protein_id[sig_idx])
      } else {
        target_genes <- as.character(integrated_results$exposure[sig_idx])
      }
    }
  }
  # 如果没有显著的，取IVW p值最小的1个基因
  if (length(target_genes) == 0 && nrow(integrated_results) > 0) {
    if ("ivw_pval" %in% colnames(integrated_results)) {
      pvals <- suppressWarnings(as.numeric(integrated_results$ivw_pval))
      best_idx <- which.min(pvals)[1]
      if ("gene_symbol" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$gene_symbol[best_idx])
      } else if ("protein_id" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$protein_id[best_idx])
      } else {
        target_genes <- as.character(integrated_results$exposure[best_idx])
      }
    } else if (nrow(integrated_results) > 0) {
      if ("gene_symbol" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$gene_symbol[1])
      } else if ("protein_id" %in% colnames(integrated_results)) {
        target_genes <- as.character(integrated_results$protein_id[1])
      } else {
        target_genes <- as.character(integrated_results$exposure[1])
      }
    }
  }
}
target_genes <- unique(target_genes)
if (length(target_genes) > 1) {
  # 报告仅展示一个示例基因的四类图形，避免正文过长
  target_genes <- target_genes[1]
}

# 图号字母映射 (A=scatter, B=forest, C=funnel, D=leaveoneout per gene)
figure_letters <- LETTERS[1:26]
plot_types <- c("scatter", "forest", "funnel", "leaveoneout")
plot_subdirs <- c("01_scatter", "02_forest", "03_funnel", "04_leaveoneout")
plot_titles <- c("散点图（SNP效应）", "森林图（SNP效应）", "漏斗图（异质性）", "逐个剔除检验图（敏感性）")
plot_descs <- list(
  "x轴代表SNP对暴露的影响，y轴代表SNP对结局的影响。每一个点代表了一个工具变量SNP，每个点上的线实际反映的是95%置信区间。横坐标是SNP对暴露因素的效应；纵坐标是SNP对结局因素的效应，两个效应之比即暴露对结局的效应，即图中彩色线的斜率。不同颜色的线表示不同的算法，不同算法的线总体斜向上，则表示随着基因水平的升高，疾病的发病风险升高；不同算法的线总体斜向下，则表示随着基因水平的升高，疾病的发病风险降低。",
  "每一条水平实线反映的是单个SNP对结局的效应。横坐标为效应值（OR），纵坐标为SNP名称。当效应值的95%置信区间不跨越0时，表明该SNP与结局存在显著关联。看最下面的红线（IVW算法），红线完全在0左边，说明基因水平增加能降低疾病的发病风险；红线完全在0右边，说明基因水平增加能增加疾病的发病风险。",
  "可以观察SNP的异质性，主要关注IVW线的左右两边的点是否大致对称。如果有特别离群的点，说明存在离群值，分析结果的可靠性需要谨慎解读。",
  "孟德尔随机化逐个剔除检验。逐步剔除每个SNP后，重新计算IVW效应值（图中每一条水平线代表剔除对应SNP后的IVW效应值和95%置信区间，最右端的菱形代表纳入全部SNP时的IVW合并效应）。若剔除某个SNP后结果发生显著改变，说明该SNP对结果影响较大，可能存在异常值。"
)

find_plots_recursive <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("\\.png$", files)]
  if (length(files) > 0) return(files[1])
  return(NULL)
}

figure_count <- 0

for (g_idx in seq_along(target_genes)) {
  gene_name <- target_genes[g_idx]
  gene_pattern <- gsub("-", "_", gene_name)
  gene_pattern <- sub("^eqtl-a-", "", gene_pattern)
  gene_pattern <- sub("^pqtl-a-", "", gene_pattern)

  for (p_idx in seq_along(plot_types)) {
    plot_type <- plot_types[p_idx]
    plot_subdir <- plot_subdirs[p_idx]

    # 尝试多种文件名模式（同时匹配 scatter.png 和 scatter_plot.png）
    base_patterns <- c(
      paste0(gene_pattern, "\\.", plot_type, "(_plot)?\\.png$"),
      paste0("[0-9]{2}\\..*", gene_pattern, ".*\\.", plot_type, "(_plot)?\\.png$"),
      paste0("[0-9]{2}\\.", gene_pattern, "\\.", plot_type, "(_plot)?\\.png$"),
      paste0(plot_type, "(_plot)?\\.png$")
    )

    found_file <- NULL
    for (p in base_patterns) {
      subdir_path <- file.path(figures_dir, plot_subdir)
      if (dir.exists(subdir_path)) {
        files <- list.files(subdir_path, pattern = p, full.names = TRUE)
        if (length(files) > 0) {
          found_file <- files[1]
          break
        }
      }
    }

    if (is.null(found_file)) {
      # 在整个figures目录递归搜索
      for (p in base_patterns) {
        found_file <- find_plots_recursive(figures_dir, p)
        if (!is.null(found_file)) break
      }
    }

    if (!is.null(found_file)) {
      figure_count <<- figure_count + 1
      fig_num <- figure_count

      img_block <- fpar(
        external_img(src = found_file, width = 5, height = 4),
        fp_p = fp_par(text.align = "center")
      )
      doc <- doc %>% body_add_fpar(value = img_block)
      doc <- doc %>% body_add_mixed_figure_caption(paste0("图", fig_num, "：", plot_titles[p_idx]))
      doc <- doc %>% body_add_mixed_figure_caption(paste0("图注：", plot_descs[[p_idx]]))

      cat("Added plot:", gene_name, "-", plot_type, "(图", fig_num, ")\n")
    }
  }
}

if (figure_count == 0) {
  doc <- doc %>% body_add_mixed_para("无可视化结果可用。")
}

# Save document
print(doc, target = opt$output)

cat("\n================================================================================\n")
cat("报告已生成：", opt$output, "\n")
cat("================================================================================\n")
