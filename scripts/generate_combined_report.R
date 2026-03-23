#!/usr/bin/env Rscript

# =============================================================================
# Combined eQTL + pQTL MR Report Generator
# =============================================================================
# Generates Word (.docx) report with both eQTL and pQTL results
# =============================================================================

library(optparse)
library(officer)
library(flextable)
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
  officer::body_add_fpar(doc, value = officer::fpar(""))
}

option_list <- list(
  make_option(c("--eqtl-dir"), type = "character", default = NULL,
              help = "eQTL results directory", dest = "eqtl_dir"),
  make_option(c("--pqtl-dir"), type = "character", default = NULL,
              help = "pQTL results directory", dest = "pqtl_dir"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output report file (.docx)"),
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
report_template <- file.path(project_root, "templates", "mr_combined_report_template.docx")

if (is.null(opt$timestamp)) {
  opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

if (is.null(opt$output)) {
  opt$output <- paste0("MR_eqtl_pqtl_", opt$timestamp, ".docx")
}

cat("================================================================================\n")
cat("  MR eQTL+pQTL 联合报告生成器\n")
cat("  时间戳:", opt$timestamp, "\n")
cat("================================================================================\n\n")

# Read eQTL results
load_mr_results <- function(result_dir) {
  candidate_files <- c(
    file.path(result_dir, "00.Complete_MR_Results.csv"),
    file.path(result_dir, "tables", "00.Complete_MR_Results.csv")
  )
  hit <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  read.csv(hit, stringsAsFactors = FALSE)
}

load_common_genes <- function(eqtl_dir, pqtl_dir) {
  roots <- unique(c(normalizePath(file.path(eqtl_dir, ".."), mustWork = FALSE),
                    normalizePath(file.path(pqtl_dir, ".."), mustWork = FALSE),
                    normalizePath(file.path(eqtl_dir, "..", ".."), mustWork = FALSE)))
  candidates <- unique(c(
    file.path(roots, "common_significant_genes_concise.csv"),
    file.path(eqtl_dir, "..", "common_significant_genes_concise.csv")
  ))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  read.csv(hit, stringsAsFactors = FALSE)
}

eqtl_results <- load_mr_results(opt$eqtl_dir)

# Read pQTL results
pqtl_results <- load_mr_results(opt$pqtl_dir)

# Read common genes
common_genes <- load_common_genes(opt$eqtl_dir, opt$pqtl_dir)

# Create Word document (use project template when available)
doc <- if (file.exists(report_template)) {
  cat("Using report template:", report_template, "\n")
  read_docx(path = report_template)
} else {
  read_docx()
}

# ========== Title ==========
doc <- doc %>%
  body_add_mixed_heading("eQTL + pQTL 孟德尔随机化联合分析报告", font_size = 16, bold = TRUE) %>%
  body_add_empty()

# ========== Methods Section ==========
doc <- doc %>%
  body_add_mixed_heading("1. 方法", font_size = 14, bold = TRUE) %>%
  body_add_empty() %>%
  body_add_mixed_heading("1.1 孟德尔随机化分析", font_size = 12, bold = TRUE) %>%
  body_add_mixed_para(
    "孟德尔随机化（MR）是一种利用遗传变异作为工具变量来推断暴露与结局之间因果关系的方法。本分析同时进行eQTL（表达数量性状位点）和pQTL（蛋白数量性状位点）两种MR分析，以识别与疾病相关的基因和蛋白。"
  ) %>%
  body_add_empty() %>%
  body_add_mixed_heading("1.2 统计方法", font_size = 12, bold = TRUE) %>%
  body_add_mixed_para(
    "1) 逆方差加权法（IVW）- 主要MR方法\n2) MR-Egger - 对水平多效性具有稳健性\n3) 加权中位数法 - 对无效工具变量具有稳健性"
  ) %>%
  body_add_mixed_para(
    "关键参数：显著性判定以IVW方法P值<0.05为阈值；报告同时汇总eQTL与pQTL结果并进行共有基因交集展示。"
  ) %>%
  body_add_empty()

# ========== Results Section ==========
doc <- doc %>%
  body_add_mixed_heading("2. 结果", font_size = 14, bold = TRUE) %>%
  body_add_empty()

# eQTL Summary
doc <- doc %>% body_add_mixed_heading("2.1 eQTL分析结果", font_size = 12, bold = TRUE)

if (!is.null(eqtl_results) && nrow(eqtl_results) > 0) {
  total_eqtl <- nrow(eqtl_results)
  sig_eqtl <- sum(!is.na(eqtl_results$ivw_pval) & as.numeric(eqtl_results$ivw_pval) < 0.05, na.rm = TRUE)

  doc <- doc %>%
    body_add_mixed_para(paste0("客观统计：分析基因总数为", total_eqtl, "，其中显著基因数（IVW p < 0.05）为", sig_eqtl, "。")) %>%
    body_add_mixed_para("关键参数：统计口径采用IVW方法结果并以P<0.05定义显著性。") %>%
    body_add_mixed_para("相关结果表：`integrated_mr_results_summary.csv`。")
} else {
  doc <- doc %>% body_add_mixed_para("无eQTL结果可用。")
}

doc <- doc %>% body_add_empty()

# pQTL Summary
doc <- doc %>% body_add_mixed_heading("2.2 pQTL分析结果", font_size = 12, bold = TRUE)

if (!is.null(pqtl_results) && nrow(pqtl_results) > 0) {
  total_pqtl <- nrow(pqtl_results)
  sig_pqtl <- sum(!is.na(pqtl_results$ivw_pval) & as.numeric(pqtl_results$ivw_pval) < 0.05, na.rm = TRUE)

  doc <- doc %>%
    body_add_mixed_para(paste0("客观统计：分析蛋白总数为", total_pqtl, "，其中显著蛋白数（IVW p < 0.05）为", sig_pqtl, "。")) %>%
    body_add_mixed_para("关键参数：统计口径采用IVW方法结果并以P<0.05定义显著性。") %>%
    body_add_mixed_para("相关结果表：`integrated_mr_results_summary.csv`。")
} else {
  doc <- doc %>% body_add_mixed_para("无pQTL结果可用。")
}

doc <- doc %>% body_add_empty()

# Common genes
doc <- doc %>% body_add_mixed_heading("2.3 eQTL与pQTL共有基因", font_size = 12, bold = TRUE)

if (!is.null(common_genes) && nrow(common_genes) > 0) {
  doc <- doc %>%
    body_add_mixed_para(paste0("客观统计：eQTL与pQTL共有显著基因数量为", nrow(common_genes), "。")) %>%
    body_add_mixed_para("关键参数：共有基因基于两类分析的显著结果集合交集定义。")

  # Create common genes table
  cols_map <- c(
    Gene = "Gene",
    eQTL_n_SNPs = "eQTL_n_SNPs",
    eQTL_OR = "eQTL_OR",
    eQTL_Pval = "eQTL_Pval",
    pQTL_n_SNPs = "pQTL_n_SNPs",
    pQTL_OR = "pQTL_OR",
    pQTL_Pval = "pQTL_Pval"
  )

  cols_exist <- names(cols_map)[names(cols_map) %in% colnames(common_genes)]

  if (length(cols_exist) > 0) {
    table_data <- common_genes[, cols_exist, drop = FALSE]
    colnames(table_data) <- cols_map[cols_exist]

    ft <- flextable(table_data) %>%
      bold(part = "header", bold = TRUE) %>%
      fontsize(size = 9, part = "all") %>%
      font(fontname = "Times New Roman", part = "all") %>%
      align(align = "center", part = "all") %>%
      align(align = "left", part = "body", j = 1) %>%
      border_remove() %>%
      border_outer(border = officer::fp_border(width = 1.5, color = "black"), part = "header") %>%
      border_inner_h(border = officer::fp_border(width = 0.5, color = "black"), part = "header") %>%
      border_outer(border = officer::fp_border(width = 1.5, color = "black"), part = "body") %>%
      set_table_properties(layout = "autofit", width = 1)

    doc <- doc %>%
      body_add_mixed_para("共有基因详细结果：") %>%
      body_add_flextable(ft) %>%
      body_add_mixed_para("相关结果表：`common_significant_genes_concise.csv`。")
  }
} else {
  doc <- doc %>% body_add_mixed_para("无共有显著基因。")
}

# ========== Visualization ==========
doc <- doc %>% body_add_empty() %>%
  body_add_mixed_heading("2.4 可视化", font_size = 12, bold = TRUE)

# Function to add plots from a results directory
add_plots_to_doc <- function(doc, results_dir, results_data, label_prefix) {
  if (is.null(results_data) || nrow(results_data) == 0) return(doc)

  # Get first gene/exposure
  first_item <- NA_character_
  for (candidate_col in c("gene_symbol", "exposure", "protein_id", "ensembl_id")) {
    if (candidate_col %in% colnames(results_data)) {
      first_item <- as.character(results_data[[candidate_col]][1])
      break
    }
  }
  if (is.null(first_item) || is.na(first_item)) {
    return(doc)
  }

  figures_dir <- file.path(results_dir, "figures")

  if (!dir.exists(figures_dir)) return(doc)

  gene_pattern <- gsub("-", "_", first_item)
  gene_pattern <- sub("^eqtl-a-", "", gene_pattern)
  gene_pattern <- sub("^pqtl-a-", "", gene_pattern)

  plot_types <- c("forest_plot", "funnel_plot", "scatter_plot", "leaveoneout_plot")
  plot_captions <- c(
    forest_plot = "展示各工具变量效应值及其置信区间，整体估计用于评估暴露与结局关系方向和强度。",
    funnel_plot = "评估效应值分布对称性，用于辅助判断发表偏倚及潜在异质性。",
    scatter_plot = "横轴为SNP对暴露的效应，纵轴为SNP对结局的效应；不同拟合线代表不同MR方法估计。",
    leaveoneout_plot = "逐个剔除工具变量后重复估计总体效应，用于评估结果稳健性和异常SNP影响。"
  )

  figure_count <- 0

  for (plot_type in plot_types) {
    patterns <- c(
      paste0(gene_pattern, "_", plot_type, "\\.png$"),
      paste0(first_item, "_", plot_type, "\\.png$"),
      paste0(".*", plot_type, "\\.png$")
    )

    found_file <- NULL
    for (p in patterns) {
      files <- list.files(figures_dir, pattern = p, full.names = TRUE)
      if (length(files) > 0) {
        found_file <- files[1]
        break
      }
    }

    if (!is.null(found_file)) {
      figure_count <- figure_count + 1

      plot_title <- switch(plot_type,
        "forest_plot" = "森林图（各SNP效应大小）",
        "funnel_plot" = "漏斗图（发表偏倚评估）",
        "scatter_plot" = "散点图（SNP对暴露与结局的效应）",
        "leaveoneout_plot" = "留一法图（敏感性分析）"
      )

      doc <- doc %>%
        body_add_mixed_figure_caption(paste0(label_prefix, figure_count, "：", plot_title)) %>%
        body_add_img(src = found_file, width = 5, height = 4, style = "Normal") %>%
        body_add_mixed_figure_caption(paste0("图注：", plot_captions[[plot_type]]))

      cat("Added", label_prefix, "plot:", plot_type, "\n")
    }
  }

  return(doc)
}

# Add eQTL plots
doc <- add_plots_to_doc(doc, opt$eqtl_dir, eqtl_results, "图eQTL")

# Add pQTL plots
doc <- add_plots_to_doc(doc, opt$pqtl_dir, pqtl_results, "图pQTL")

# Save document
print(doc, target = opt$output)

cat("\n================================================================================\n")
cat("联合报告已生成：", opt$output, "\n")
cat("================================================================================\n")
