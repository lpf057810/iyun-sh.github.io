#!/usr/bin/env Rscript
# ==============================================================================
# MCPcounter 免疫浸润分析
# ==============================================================================
# 功能：使用 MCPcounter 算法计算免疫细胞丰度评分
# 输入：表达矩阵、分组信息
# 输出：免疫细胞丰度评分、统计比较、相关性分
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
script_dir <- if (interactive()) {
  getwd()
} else {
  # 获取脚本所在目录
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("^--file=", args)])
  dirname(normalizePath(script_path))
}
source(file.path(script_dir, "immune_functions.R"))

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
              default="MCPcounter_output"),
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

cat("\n=== MCPcounter Immune Infiltration Analysis ===\n")
cat(sprintf("Expression file: %s\n", opt$expr))
cat(sprintf("Group file: %s\n", opt$group))
if (!is.null(opt$genes)) {
  cat(sprintf("Gene list file: %s\n", opt$genes))
}
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

cache_file <- file.path(cache_dir, "mcpcounter.qs2")

formatted_time()
cat(sprintf("Checking cache file: %s\n", cache_file))

if (file.exists(cache_file) && !opt$force) {
  cat("Found cached result, loading...\n")
  formatted_time()
  mcpcounter_res <- qs_read(cache_file)
  cat("Loaded cached MCPcounter result\n")
} else {
  cat("Running MCPcounter...\n")
  formatted_time()
  mcpcounter_res <- deconvo_tme(eset = dat, method = "mcpcounter", arrays = opt$arrays)
  
  cat("Saving to cache...\n")
  qs_save(mcpcounter_res, cache_file)
  cat(sprintf("Saved to: %s\n", cache_file))
  formatted_time()
}

names(mcpcounter_res) <- gsub("_MCPcounter$", "", names(mcpcounter_res))
fwrite(mcpcounter_res, file.path(temp_outdir, "01.mcpcounter_res.csv"))
cat(sprintf("Saved MCPcounter results: %s\n", file.path(temp_outdir, "01.mcpcounter_res.csv")))

cell_cols <- setdiff(colnames(mcpcounter_res), "ID")
if (length(cell_cols) > 0) {
  mypalette <- colorRampPalette(brewer.pal(min(length(cell_cols), 9), "Set2"))

  p <- mcpcounter_res %>%
    gather(cell_type, fraction, -ID) %>%
    merge(df.group, by = 'ID') %>%
    ggplot(aes(x = ID, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Abundance Score', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top') +
    scale_fill_manual(values = mypalette(length(cell_cols))) +
    facet_grid(~ group, scales = "free", space = "free")

  save_plot(file.path(temp_outdir, '01.mcpcounter.stacked_bar.pdf'), p, width = 10, height = 6, both = TRUE)
  cat("Saved stacked bar plot\n")

  tryCatch({
    bar_data <- mcpcounter_res %>%
      gather(cell_type, fraction, -ID) %>%
      merge(df.group, by = "ID")
    bar_summary <- bar_data %>%
      group_by(group, cell_type) %>%
      summarise(mean_score = mean(fraction, na.rm = TRUE),
                sd_score = sd(fraction, na.rm = TRUE),
                n = n(), .groups = "drop") %>%
      arrange(group, desc(mean_score))
    fwrite(bar_summary, file.path(temp_outdir, "01.mcpcounter_bar_stats.csv"))
    cat("Saved bar stats CSV\n")
  }, error = function(e) cat("Warning: bar stats generation failed:", e$message, "\n"))
}

res_mcpcounter2 <- mcpcounter_res %>%
  column_to_rownames(var = 'ID') %>%
  t() %>%
  as.data.frame()
res_mcpcounter2 <- res_mcpcounter2[rowSums(res_mcpcounter2) > 0, ]

dat_mcpcounter <- res_mcpcounter2 %>%
  t %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")
dat_mcpcounter <- merge(df.group, dat_mcpcounter, by = "ID")
dat_mcpcounter2 <- tidyr::gather(dat_mcpcounter, ImmuneCell, Score, -c("ID", "group"))

stat_mcpcounter <- dat_mcpcounter2 %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

fwrite(stat_mcpcounter, file.path(temp_outdir, '02.stat.mcpcounter.csv'))
cat("Saved statistical analysis\n")

DE_mcpcounter <- stat_mcpcounter

violin_mcpcounter <- dat_mcpcounter2

p1 <- ggplot(violin_mcpcounter, aes(x = ImmuneCell, y = Score, fill = group)) +
  geom_boxplot(width = 0.5,
               alpha = 0.8,
               position = position_dodge(0.9),
               outlier.shape = NA) +
  scale_fill_manual(values = group_color[1:2]) +
  annotate(geom = "text", x = DE_mcpcounter$ImmuneCell,
           y = max(violin_mcpcounter$Score), size = 3, 
           label = as.character(DE_mcpcounter$p.signif)) +
  theme_classic() +
  theme(legend.position = "top") +
  theme(axis.title.x = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5, color = 'black'),
        axis.title.y = element_text(size = 15, color = 'black'),
        axis.text.y = element_text(size = 15, color = 'black')) +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14))

save_plot(file.path(temp_outdir, '03.DE.mcpcounter.boxplot.pdf'), p1, width = 12, height = 7, both = TRUE)
cat("Saved boxplot\n")

if (!is.null(opt$genes)) {
  gene <- fread(opt$genes)$gene
  cat(sprintf("Loaded %d genes for correlation analysis\n", length(gene)))
  
  hub_exp <- dat[gene, ]
  cortiic <- res_mcpcounter2
  cortiic <- cortiic[rowSums(cortiic) > 0, ]
  DE_mcpcounter <- stat_mcpcounter[which(stat_mcpcounter$p < 0.05), ]
  
  if (nrow(DE_mcpcounter) > 0) {
    cortiic <- cortiic[DE_mcpcounter$ImmuneCell, ]
    if (ncol(cortiic) > 0) {
      cortiic <- cortiic[, colnames(hub_exp)]
      
      if (ncol(cortiic) > 0 && nrow(cortiic) > 0) {
        res <- get_corr(t(hub_exp), t(cortiic), method = 'spearman', cor_val = 0.3)
        
        colnames(res$res)[1:2] <- c("gene", "ImmuneCell")
        fwrite(res$res, file.path(temp_outdir, '04.Correlation.gene.csv'))
        cat("Saved correlation analysis\n")
        
        res$p <- res$p + theme(axis.text = element_text(size = 10))
        save_plot(file.path(temp_outdir, '04.Correlation.gene.pdf'), res$p, width = 12, height = 7, both = TRUE)
        cat("Saved correlation plot\n")
      }
    }
  }
}

formatted_time()
cat("MCPcounter analysis completed!\n\n")
