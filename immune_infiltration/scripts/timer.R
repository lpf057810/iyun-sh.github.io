#!/usr/bin/env Rscript
# ==============================================================================
# TIMER 免疫浸润分析
# ==============================================================================
# 功能：使用 TIMER 算法计算免疫细胞浸润水平
# 输入：表达矩阵、分组信息、肿瘤类型
# 输出：免疫细胞浸润水平、统计比较、相关性分析
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
  make_option(c("--genes"), type="character", 
              help="Gene list file path (CSV format with gene column) [default: %default]", 
              default=NULL),
  make_option(c("-o", "--output"), type="character", 
              help="Output directory path [default: %default]", 
              default="TIMER_output"),
  make_option(c("--arrays"), type="logical", 
              help="Whether the data is from microarray [default: %default]", 
              default=FALSE),
  make_option(c("--tissue"), type="character", 
              help="Tissue type for TIMER (e.g., LUAD, LUSC, GBM, etc.) [default: %default]", 
              default="LUAD"),
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

cat("\n=== TIMER Immune Infiltration Analysis ===\n")
cat(sprintf("Expression file: %s\n", opt$expr))
cat(sprintf("Group file: %s\n", opt$group))
if (!is.null(opt$genes)) {
  cat(sprintf("Gene list file: %s\n", opt$genes))
}
cat(sprintf("Output directory: %s\n", opt$output))
cat(sprintf("Tissue type: %s\n", opt$tissue))
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

cache_file <- file.path(cache_dir, sprintf("timer_%s.qs2", opt$tissue))

formatted_time()
cat(sprintf("Checking cache file: %s\n", cache_file))

if (file.exists(cache_file) && !opt$force) {
  cat("Found cached result, loading...\n")
  formatted_time()
  timer_res <- qs_read(cache_file)
  cat("Loaded cached TIMER result\n")
} else {
  cat("Running TIMER...\n")
  formatted_time()
  group_list <- rep(opt$tissue, dim(dat)[2])
  timer_res <- deconvo_tme(eset = dat, method = "timer", 
                          arrays = opt$arrays, 
                          group_list = group_list)
  
  cat("Saving to cache...\n")
  qs_save(timer_res, cache_file)
  cat(sprintf("Saved to: %s\n", cache_file))
  formatted_time()
}

names(timer_res) <- gsub("_TIMER$", "", names(timer_res))
fwrite(timer_res, file.path(temp_outdir, "01.timer_res.csv"))
cat(sprintf("Saved TIMER results: %s\n", file.path(temp_outdir, "01.timer_res.csv")))

cell_cols <- setdiff(colnames(timer_res), c("ID", "TumorPurity"))
if (length(cell_cols) > 0) {
  mypalette <- colorRampPalette(brewer.pal(min(length(cell_cols), 8), "Dark2"))

  # Handle case where TumorPurity doesn't exist
  if ("TumorPurity" %in% colnames(timer_res)) {
    p_data <- timer_res %>% gather(cell_type, fraction, -ID, -TumorPurity)
  } else {
    p_data <- timer_res %>% gather(cell_type, fraction, -ID)
  }

  p <- p_data %>%
    merge(df.group, by = 'ID') %>%
    ggplot(aes(x = ID, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Relative Percent', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    scale_fill_manual(values = mypalette(length(cell_cols))) +
    facet_grid(~ group, scales = "free", space = "free")

  save_plot('01.timer.stacked_bar', p, outdir = temp_outdir, width = 10, height = 6, both = TRUE)
  cat("Saved stacked bar plot\n")
}

# Handle case where TumorPurity doesn't exist
if ("TumorPurity" %in% colnames(timer_res)) {
  res_timer2 <- timer_res %>% dplyr::select(-c("TumorPurity"))
} else {
  res_timer2 <- timer_res
}
res_timer2 <- res_timer2 %>%
  column_to_rownames(var = 'ID') %>%
  t() %>%
  as.data.frame()
res_timer2 <- res_timer2[rowSums(res_timer2) > 0, ]

dat_timer <- res_timer2 %>%
  t %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")
dat_timer <- merge(df.group, dat_timer, by = "ID")
dat_timer2 <- tidyr::gather(dat_timer, ImmuneCell, Score, -c("ID", "group"))

# Filter out NA values and cell types with no variance
dat_timer2 <- dat_timer2 %>% filter(!is.na(Score) & !is.nan(Score) & Score != Inf & Score != -Inf)

# Filter to only include cell types that have data in both groups
cell_types_valid <- dat_timer2 %>%
  group_by(ImmuneCell) %>%
  summarise(n_groups = n_distinct(group), .groups = 'drop') %>%
  filter(n_groups >= 2) %>%
  pull(ImmuneCell)

dat_timer2 <- dat_timer2 %>% filter(ImmuneCell %in% cell_types_valid)

if (nrow(dat_timer2) > 0 && length(unique(dat_timer2$group)) >= 2) {
  stat_timer <- dat_timer2 %>%
    group_by(ImmuneCell) %>%
    wilcox_test(Score ~ group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p")
} else {
  stat_timer <- data.frame()
}

fwrite(stat_timer, file.path(temp_outdir, '02.stat.timer.csv'))
cat("Saved statistical analysis\n")

DE_timer <- stat_timer

violin_timer <- dat_timer2

if (nrow(violin_timer) > 0 && nrow(DE_timer) > 0) {
  p1 <- ggplot(violin_timer, aes(x = ImmuneCell, y = Score, fill = group)) +
    geom_boxplot(width = 0.5,
                 alpha = 0.8,
                 position = position_dodge(0.9),
                 outlier.shape = NA) +
    scale_fill_manual(values = group_color[1:2]) +
    annotate(geom = "text", x = DE_timer$ImmuneCell,
             y = max(violin_timer$Score), size = 3,
             label = as.character(DE_timer$p.signif)) +
  theme_classic() +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5, color = 'black'),
        axis.title.y = element_text(size = 15, color = 'black'),
        axis.text.y = element_text(size = 15, color = 'black')) +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14))

save_plot('03.DE.timer.boxplot', p1, outdir = temp_outdir, width = 12, height = 7, both = TRUE)
  cat("Saved boxplot\n")
} else {
  cat("Skipped boxplot - no valid data for statistical comparison\n")
}

if (!is.null(opt$genes)) {
  gene <- fread(opt$genes)$gene
  cat(sprintf("Loaded %d genes for correlation analysis\n", length(gene)))
  
  hub_exp <- dat[gene, ]
  cortiic <- res_timer2
  cortiic <- cortiic[rowSums(cortiic) > 0, ]
  DE_timer <- stat_timer[which(stat_timer$p < 0.05), ]
  
  if (nrow(DE_timer) > 0) {
    cortiic <- cortiic[DE_timer$ImmuneCell, ]
    if (ncol(cortiic) > 0) {
      cortiic <- cortiic[, colnames(hub_exp)]
      
      if (ncol(cortiic) > 0 && nrow(cortiic) > 0) {
        res <- get_corr(t(hub_exp), t(cortiic), method = 'spearman', cor_val = 0.3)
        
        colnames(res$res)[1:2] <- c("gene", "ImmuneCell")
        fwrite(res$res, file.path(temp_outdir, '04.Correlation.gene.csv'))
        cat("Saved correlation analysis\n")
        
        res$p <- res$p + theme(axis.text = element_text(size = 10))
        save_plot('04.Correlation.gene', res$p, outdir = temp_outdir, width = 12, height = 7, both = TRUE)
        cat("Saved correlation plot\n")
      }
    }
  }
}

formatted_time()
cat("TIMER analysis completed!\n\n")
