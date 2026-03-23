#!/usr/bin/env Rscript
# ==============================================================================
# ESTIMATE 免疫浸润分析
# ==============================================================================
# 功能：使用 ESTIMATE 算法计算肿瘤微环境评分
# 输入：表达矩阵、分组信息
# 输出：ImmuneScore、StromalScore、ESTIMATEScore
#
# 作者: [项目团队]
# 创建: 2026-02-03
# 更新: 2026-02-15
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(IOBR)
  library(glue)
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(qs2)
  library(ggplot2)
  library(rstatix)
  library(RColorBrewer)
})

# 设置镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 加载函数库
# Python脚本已设置cwd为script_dir，所以直接用getwd()
# 获取脚本所在目录
script_dir <- if (interactive()) {
  getwd()
} else {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(normalizePath(script_path))
}
source(file.path(script_dir, "./immune_functions.R"))

option_list <- list(
  make_option(c("-e", "--expr"), type="character", 
              help="Expression matrix file path (CSV format with gene symbols as first column) [default: %default]", 
              default="expr_tpm.csv"),
  make_option(c("-g", "--group"), type="character", 
              help="Group information file path (CSV format with sample and group columns) [default: %default]", 
              default="group.csv"),
  make_option(c("-o", "--output"), type="character", 
              help="Output directory path [default: %default]", 
              default="ESTIMATE_output"),
  make_option(c("--arrays"), type="logical", 
              help="Whether the data is from microarray [default: %default]", 
              default=FALSE),
  make_option(c("--control"), type="character", 
              help="Control group name [default: %default]", 
              default="Normal"),
  make_option(c("--cache"), type="character", 
              help="Cache directory for intermediate files [default: %default]", 
              default="cache"),
  make_option(c("--force"), type="logical", 
              help="Force re-run calculation even if cache exists [default: %default]", 
              default=FALSE)
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

cat("\n=== ESTIMATE Immune Infiltration Analysis ===\n")
cat(sprintf("Expression file: %s\n", opt$expr))
cat(sprintf("Group file: %s\n", opt$group))
cat(sprintf("Output directory: %s\n", opt$output))
cat("\n")

formatted_time()

dat <- fread(opt$expr) %>% column_to_rownames("SYMBOL")
cat(sprintf("Loaded expression matrix: %d genes x %d samples\n", nrow(dat), ncol(dat)))

df.group <- fread(opt$group)
colnames(df.group)[1] <- "ID"
cat(sprintf("Loaded group information: %d samples\n", nrow(df.group)))

temp_outdir <- opt$output
if (!dir.exists(temp_outdir)) {
  dir.create(temp_outdir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", temp_outdir))
}

cache_dir <- opt$cache
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

cache_file <- file.path(cache_dir, "estimate.qs2")

formatted_time()
cat(sprintf("Checking cache file: %s\n", cache_file))

if (file.exists(cache_file) && !opt$force) {
  cat("Found cached result, loading...\n")
  formatted_time()
  estimate_res <- qs_read(cache_file)
  cat("Loaded cached ESTIMATE result\n")
} else {
  cat("Running ESTIMATE...\n")
  formatted_time()
  estimate_res <- deconvo_tme(eset = dat, method = "estimate", arrays = opt$arrays)
  
  cat("Saving to cache...\n")
  qs_save(estimate_res, cache_file)
  cat(sprintf("Saved to: %s\n", cache_file))
  formatted_time()
}

names(estimate_res) <- gsub("_estimate$", "", names(estimate_res))
fwrite(estimate_res, file.path(temp_outdir, "01.estimate_res.csv"))
cat(sprintf("Saved ESTIMATE results: %s\n", file.path(temp_outdir, "01.estimate_res.csv")))

cell_cols <- setdiff(colnames(estimate_res), c("ID", "TumorPurity"))
if (length(cell_cols) > 0) {
  mypalette <- colorRampPalette(brewer.pal(min(length(cell_cols), 8), "Set1"))

  p <- estimate_res %>%
    gather(cell_type, fraction, -ID, -TumorPurity) %>%
    merge(df.group, by = 'ID') %>%
    ggplot(aes(x = ID, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Score', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    scale_fill_manual(values = mypalette(length(cell_cols))) +
    facet_grid(~ group, scales = "free", space = "free")

  save_plot('01.estimate.stacked_bar', p, outdir = temp_outdir, width = 10, height = 6, both = TRUE)
  cat("Saved stacked bar plot\n")

  tryCatch({
    bar_data <- estimate_res %>%
      gather(cell_type, fraction, -ID, -TumorPurity) %>%
      merge(df.group, by = "ID")
    bar_summary <- bar_data %>%
      group_by(group, cell_type) %>%
      summarise(mean_score = mean(fraction, na.rm = TRUE),
                sd_score = sd(fraction, na.rm = TRUE),
                n = n(), .groups = "drop") %>%
      arrange(group, desc(mean_score))
    fwrite(bar_summary, file.path(temp_outdir, "01.estimate_bar_stats.csv"))
    cat("Saved bar stats CSV\n")
  }, error = function(e) cat("Warning: bar stats generation failed:", e$message, "\n"))
}

res_estimate2 <- estimate_res %>%
  dplyr::select(-c("TumorPurity")) %>%
  column_to_rownames(var = 'ID') %>%
  t() %>%
  as.data.frame()

dat_estimate <- res_estimate2 %>%
  t %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")
dat_estimate <- merge(df.group, dat_estimate, by = "ID")
dat_estimate2 <- tidyr::gather(dat_estimate, ImmuneCell, Score, -c("ID", "group"))

# Filter out NA values and cell types with no variance
dat_estimate2 <- dat_estimate2 %>% filter(!is.na(Score) & !is.nan(Score) & Score != Inf & Score != -Inf)

# Filter to only include cell types that have data in both groups
cell_types_valid <- dat_estimate2 %>%
  group_by(ImmuneCell) %>%
  summarise(n_groups = n_distinct(group), .groups = 'drop') %>%
  filter(n_groups >= 2) %>%
  pull(ImmuneCell)

dat_estimate2 <- dat_estimate2 %>% filter(ImmuneCell %in% cell_types_valid)

if (nrow(dat_estimate2) > 0 && length(unique(dat_estimate2$group)) >= 2) {
  stat_estimate <- dat_estimate2 %>%
    group_by(ImmuneCell) %>%
    wilcox_test(Score ~ group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p")
} else {
  stat_estimate <- data.frame()
}

fwrite(stat_estimate, file.path(temp_outdir, '02.stat.estimate.csv'))
cat("Saved statistical analysis\n")

DE_estimate <- stat_estimate

violin_estimate <- dat_estimate2

if (nrow(violin_estimate) > 0 && nrow(DE_estimate) > 0) {
  p1 <- ggplot(violin_estimate, aes(x = ImmuneCell, y = Score, fill = group)) +
    geom_boxplot(width = 0.5,
                 alpha = 0.8,
                 position = position_dodge(0.9),
                 outlier.shape = NA) +
    scale_fill_manual(values = group_color[1:2]) +
    annotate(geom = "text", x = DE_estimate$ImmuneCell,
             y = max(violin_estimate$Score), size = 3,
             label = as.character(DE_estimate$p.signif)) +
    theme_classic() +
    theme(legend.position = "top") +
    theme(axis.title.x = element_text(size = 15, color = 'black'),
          axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5, color = 'black'),
          axis.title.y = element_text(size = 15, color = 'black'),
          axis.text.y = element_text(size = 15, color = 'black')) +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14))

  save_plot('03.DE.estimate.boxplot', p1, outdir = temp_outdir, width = 12, height = 7, both = TRUE)
  cat("Saved boxplot\n")
} else {
  cat("Skipped boxplot - no valid data for statistical comparison\n")
}

formatted_time()
cat("ESTIMATE analysis completed!\n\n")
