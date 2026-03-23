#!/usr/bin/env Rscript

# =============================================================================
# Find Common Genes between eQTL and pQTL
# =============================================================================
# Takes significant results from both analyses and finds common genes
# =============================================================================

library(optparse)

option_list <- list(
  make_option(c("--eqtl-dir"), type = "character", default = NULL,
              help = "eQTL results directory", dest = "eqtl_dir"),
  make_option(c("--pqtl-dir"), type = "character", default = NULL,
              help = "pQTL results directory", dest = "pqtl_dir"),
  make_option(c("-o", "--output"), type = "character", default = ".",
              help = "Output directory [default %default]"),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Default timestamp
if (is.null(opt$timestamp)) {
  opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

cat("================================================================================\n")
cat("  寻找eQTL和pQTL共有基因\n")
cat("  时间戳:", opt$timestamp, "\n")
cat("================================================================================\n\n")

load_significant_set <- function(result_dir, label) {
  file_07 <- file.path(result_dir, "07.MR_res_gene.csv")
  file_06 <- file.path(result_dir, "06.Three_Filter_Results.csv")
  file_00 <- file.path(result_dir, "00.Complete_MR_Results.csv")

  if (file.exists(file_07)) {
    df <- read.csv(file_07, stringsAsFactors = FALSE)
    if (!"gene_symbol" %in% colnames(df) && "Exposure" %in% colnames(df)) df$gene_symbol <- df$Exposure
    if (!"n_snps" %in% colnames(df)) df$n_snps <- NA
    if (!"ivw_or" %in% colnames(df)) df$ivw_or <- NA
    if (!"ivw_pval" %in% colnames(df)) df$ivw_pval <- NA
    if (!"ivw_direction" %in% colnames(df)) {
      df$ivw_direction <- ifelse(suppressWarnings(as.numeric(df$ivw_or)) > 1, "Risk", "Protective")
    }
    cat(label, "显著基因来源: 07.MR_res_gene.csv\n")
    return(df[, c("gene_symbol", "n_snps", "ivw_or", "ivw_pval", "ivw_direction"), drop = FALSE])
  }

  if (file.exists(file_06) && file.exists(file_00)) {
    tf <- read.csv(file_06, stringsAsFactors = FALSE)
    all <- read.csv(file_00, stringsAsFactors = FALSE)
    exp_col <- if ("Exposure" %in% colnames(tf)) "Exposure" else if ("gene_symbol" %in% colnames(tf)) "gene_symbol" else NA_character_
    pass_col <- if ("Three_filter_pass" %in% colnames(tf)) "Three_filter_pass" else NA_character_
    if (!is.na(exp_col) && !is.na(pass_col)) {
      tf_pass <- tf[toupper(trimws(as.character(tf[[pass_col]]))) %in% c("YES", "TRUE", "1"), , drop = FALSE]
      keep <- unique(as.character(tf_pass[[exp_col]]))
      all <- all[as.character(all$gene_symbol) %in% keep, , drop = FALSE]
      if (!"ivw_direction" %in% colnames(all)) {
        all$ivw_direction <- ifelse(suppressWarnings(as.numeric(all$ivw_or)) > 1, "Risk", "Protective")
      }
      cat(label, "显著基因来源: 06.Three_Filter_Results.csv + 00.Complete_MR_Results.csv\n")
      return(all[, c("gene_symbol", "n_snps", "ivw_or", "ivw_pval", "ivw_direction"), drop = FALSE])
    }
  }

  cat(label, "显著结果文件不存在（需 07 或 06+00）\n")
  data.frame(gene_symbol = character(0), stringsAsFactors = FALSE)
}

# Read eQTL/pQTL significant results
eqtl_sig <- load_significant_set(opt$eqtl_dir, "eQTL")
pqtl_sig <- load_significant_set(opt$pqtl_dir, "pQTL")
eqtl_genes <- unique(eqtl_sig$gene_symbol)
pqtl_genes <- unique(pqtl_sig$gene_symbol)
cat("eQTL显著基因数量:", length(eqtl_genes), "\n")
cat("pQTL显著基因数量:", length(pqtl_genes), "\n")

# Find common genes
common_genes <- intersect(eqtl_genes, pqtl_genes)
cat("\n共有显著基因数量:", length(common_genes), "\n")

if (length(common_genes) > 0) {
  cat("共有显著基因:\n")
  for (g in common_genes) {
    cat("  -", g, "\n")
  }

  # Create combined results
  if (exists("eqtl_sig") && exists("pqtl_sig")) {
    # Extract eQTL results for common genes
    eqtl_common <- eqtl_sig[eqtl_sig$gene_symbol %in% common_genes, ]
    eqtl_common$analysis_type <- "eQTL"

    # Extract pQTL results for common genes
    pqtl_common <- pqtl_sig[pqtl_sig$gene_symbol %in% common_genes, ]
    pqtl_common$analysis_type <- "pQTL"

    # Combine results
    combined_results <- rbind(eqtl_common, pqtl_common)

    # Save full results
    output_file <- file.path(opt$output, "common_significant_genes.csv")
    write.csv(combined_results, output_file, row.names = FALSE)
    cat("\n完整结果已保存到:", output_file, "\n")

    # Create concise version
    concise_results <- data.frame(
      Gene = common_genes,
      eQTL_n_SNPs = sapply(common_genes, function(g) {
        eqtl_row <- eqtl_sig[eqtl_sig$gene_symbol == g, ]
        if (nrow(eqtl_row) > 0) eqtl_row$n_snps[1] else NA
      }),
      eQTL_OR = sapply(common_genes, function(g) {
        eqtl_row <- eqtl_sig[eqtl_sig$gene_symbol == g, ]
        if (nrow(eqtl_row) > 0) eqtl_row$ivw_or[1] else NA
      }),
      eQTL_Pval = sapply(common_genes, function(g) {
        eqtl_row <- eqtl_sig[eqtl_sig$gene_symbol == g, ]
        if (nrow(eqtl_row) > 0) eqtl_row$ivw_pval[1] else NA
      }),
      eQTL_Direction = sapply(common_genes, function(g) {
        eqtl_row <- eqtl_sig[eqtl_sig$gene_symbol == g, ]
        if (nrow(eqtl_row) > 0) eqtl_row$ivw_direction[1] else NA
      }),
      pQTL_n_SNPs = sapply(common_genes, function(g) {
        pqtl_row <- pqtl_sig[pqtl_sig$gene_symbol == g, ]
        if (nrow(pqtl_row) > 0) pqtl_row$n_snps[1] else NA
      }),
      pQTL_OR = sapply(common_genes, function(g) {
        pqtl_row <- pqtl_sig[pqtl_sig$gene_symbol == g, ]
        if (nrow(pqtl_row) > 0) pqtl_row$ivw_or[1] else NA
      }),
      pQTL_Pval = sapply(common_genes, function(g) {
        pqtl_row <- pqtl_sig[pqtl_sig$gene_symbol == g, ]
        if (nrow(pqtl_row) > 0) pqtl_row$ivw_pval[1] else NA
      }),
      pQTL_Direction = sapply(common_genes, function(g) {
        pqtl_row <- pqtl_sig[pqtl_sig$gene_symbol == g, ]
        if (nrow(pqtl_row) > 0) pqtl_row$ivw_direction[1] else NA
      }),
      stringsAsFactors = FALSE
    )

    concise_file <- file.path(opt$output, "common_significant_genes_concise.csv")
    write.csv(concise_results, concise_file, row.names = FALSE)
    cat("简洁结果已保存到:", concise_file, "\n")

    # Print concise results
    cat("\n共有显著基因简洁结果:\n")
    print(concise_results)
  }

  # Save gene list
  genes_df <- data.frame(gene_symbol = common_genes)
  write.csv(genes_df, file.path(opt$output, "eqtl_pqtl_common_genes.csv"), row.names = FALSE)
  cat("\n基因列表已保存到: eqtl_pqtl_common_genes.csv\n")
} else {
  cat("没有找到共有显著基因\n")
}

cat("\n✓ 分析完成!\n")

# 保存qs2与sessionInfo（标准化产物）
if (requireNamespace("qs", quietly = TRUE)) {
  qs_file <- file.path(opt$output, "common_genes_node.qs2")
  qs::qsave(list(eqtl_sig = eqtl_sig, pqtl_sig = pqtl_sig, common_genes = common_genes), qs_file)
  cat("QS2快照已保存:", qs_file, "\n")
}
log_dir <- file.path(opt$output, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
session_file <- file.path(log_dir, sprintf("common_genes_sessionInfo.%s.txt", opt$timestamp))
writeLines(capture.output(sessionInfo()), session_file)
cat("SessionInfo已保存:", session_file, "\n")
