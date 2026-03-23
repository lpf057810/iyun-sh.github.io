# =============================================================================
# 共定位分析 (Colocalization Analysis) 自动化脚本 - 带p值筛选修正版
# =============================================================================
# 基于MR分析结果，自动筛选基因进行共定位分析
# 流程：正向MR显著 → 反向MR排除反向因果 → 共定位分析
# 关键修改：添加eQTL/GWAS p值筛选步骤，减少噪声
# =============================================================================

# ============================================================================
# 1. 参数配置（命令行）
# ============================================================================
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = ".",
              help = "MR结果目录（eqtl结果目录，或包含eqtl子目录的组合目录）"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "共定位输出目录"),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL,
              help = "运行时间戳（用于流程兼容，可选）"),
  make_option(c("--eqtl-file"), type = "character",
              default = "/media/desk16/iyunlpf/TSG/TWAS/data/eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
              help = "本地eQTL汇总文件路径"),
  make_option(c("-d", "--outcome"), type = "character", default = "ebi-a-GCST90018894",
              help = "结局GWAS ID"),
  make_option(c("--region-radius-kb"), type = "integer", default = 50,
              help = "基因区域半径（kb）"),
  make_option(c("-p", "--p-thresh"), type = "double", default = 0.05,
              help = "SNP筛选p值阈值"),
  make_option(c("--pp4-threshold"), type = "double", default = 0.6,
              help = "共定位显著阈值"),
  make_option(c("--gwas-use-proxies"), type = "integer", default = 0,
              help = "GWAS是否使用proxy（0/1）"),
  make_option(c("--use-snp-filter"), type = "integer", default = 1,
              help = "是否启用SNP筛选（0/1）"),
  make_option(c("--case-prevalence"), type = "double", default = 0.4,
              help = "病例对照患病率参数")
)

opt <- parse_args(OptionParser(option_list = option_list))

MR_DIR <- normalizePath(opt$input, mustWork = FALSE)
if (!file.exists(file.path(MR_DIR, "00.Complete_MR_Results.csv")) &&
    dir.exists(file.path(MR_DIR, "eqtl")) &&
    file.exists(file.path(MR_DIR, "eqtl", "00.Complete_MR_Results.csv"))) {
  MR_DIR <- file.path(MR_DIR, "eqtl")
}

COLOC_DIR <- if (!is.null(opt$output) && nzchar(opt$output)) {
  opt$output
} else {
  file.path(dirname(MR_DIR), "coloc")
}
dir.create(COLOC_DIR, showWarnings = FALSE, recursive = TRUE)

EQTL_FILE <- opt$`eqtl-file`
REGION_RADIUS_KB <- as.integer(opt$`region-radius-kb`)
OUTCOME_ID <- opt$outcome
P_THRESHOLD <- as.numeric(opt$`p-thresh`)
PP4_THRESHOLD <- as.numeric(opt$`pp4-threshold`)
GWAS_USE_PROXIES <- as.integer(opt$`gwas-use-proxies`)
USE_SNP_FILTER <- as.integer(opt$`use-snp-filter`) == 1L
CASE_PREVALENCE <- as.numeric(opt$`case-prevalence`)

# ============================================================================
# 2. 加载必要的R包
# ============================================================================

library(tidyverse)
library(coloc)
library(ieugwasr)
library(TwoSampleMR)
library(data.table)

cat("============================================================\n")
cat("共定位分析 - 带p值筛选修正版\n")
cat("============================================================\n\n")
cat("MR目录:", MR_DIR, "\n")
cat("输出目录:", COLOC_DIR, "\n")

if (!file.exists(EQTL_FILE)) {
  stop("未找到eQTL文件: ", EQTL_FILE)
}

# ============================================================================
# 3. 读取并筛选基因
# ============================================================================

cat("步骤1: 读取MR分析结果并筛选基因\n")
cat("------------------------------------------------------------\n")

# 3.1 读取正向MR结果（优先07，回退00）
mr_res_file <- if (file.exists(file.path(MR_DIR, "07.MR_res_gene.csv"))) {
  file.path(MR_DIR, "07.MR_res_gene.csv")
} else {
  file.path(MR_DIR, "00.Complete_MR_Results.csv")
}

forward_results <- read.csv(mr_res_file, stringsAsFactors = FALSE)
cat("  使用结果文件:", mr_res_file, "\n")

cat("  正向MR分析基因数:", nrow(forward_results), "\n")

# 兼容列名
if (!"gene_symbol" %in% colnames(forward_results) && "Exposure" %in% colnames(forward_results)) {
  forward_results$gene_symbol <- forward_results$Exposure
}
if (!"ensembl_id" %in% colnames(forward_results) && "gene_symbol" %in% colnames(forward_results)) {
  forward_results$ensembl_id <- forward_results$gene_symbol
}
if (!"exposure_id" %in% colnames(forward_results) && "ensembl_id" %in% colnames(forward_results)) {
  forward_results$exposure_id <- paste0("eqtl-a-", forward_results$ensembl_id)
}

# 筛选显著基因 (IVW p < 0.05)；07文件通常已是筛选后结果
if ("ivw_pval" %in% colnames(forward_results)) {
  forward_results$ivw_pval_numeric <- suppressWarnings(as.numeric(forward_results$ivw_pval))
  significant_genes <- forward_results[!is.na(forward_results$ivw_pval_numeric) & forward_results$ivw_pval_numeric < 0.05, ]
} else {
  significant_genes <- forward_results
}
cat("  正向MR显著基因数 (p < 0.05):", nrow(significant_genes), "\n")

# 3.2 筛选通过Steiger因果方向检验的基因（若有该列）
cat("\n  Steiger因果方向检验筛选...\n")
if (all(c("steiger_correct_direction", "steiger_pval") %in% colnames(significant_genes))) {
  significant_genes$steiger_pval_numeric <- suppressWarnings(as.numeric(gsub("<", "", significant_genes$steiger_pval)))
  cat("  - steiger_correct_direction为TRUE:",
      sum(significant_genes$steiger_correct_direction == TRUE, na.rm = TRUE), "\n")
  cat("  - steiger_pval非NA:",
      sum(!is.na(significant_genes$steiger_pval_numeric)), "\n")
  final_genes <- significant_genes[
    !is.na(significant_genes$steiger_correct_direction) &
      significant_genes$steiger_correct_direction == TRUE &
      !is.na(significant_genes$steiger_pval_numeric) &
      significant_genes$steiger_pval_numeric < 0.05,
  ]
} else {
  final_genes <- significant_genes
  cat("  - 未检测到Steiger列，使用当前显著基因集继续共定位。\n")
}

cat("  通过Steiger检验的显著基因数:", nrow(final_genes), "\n")

cat("\n  最终用于共定位分析的基因数:", nrow(final_genes), "\n")

# 测试模式：只运行前2个基因
TEST_MODE <- FALSE
TEST_N <- 2
if (TEST_MODE) {
  final_genes <- final_genes[1:min(TEST_N, nrow(final_genes)), ]
  cat("  [测试模式] 只运行前", nrow(final_genes), "个基因\n")
}

if (nrow(final_genes) == 0) {
  cat("\n  警告: 没有通过筛选的基因!\n")
  cat("  尝试使用所有正向MR显著基因...\n")
  
  # 如果没有交集，使用所有正向显著基因
  final_genes <- significant_genes
  cat("  使用所有正向MR显著基因数:", nrow(final_genes), "\n")
}

# 显示最终选择的基因
if (nrow(final_genes) > 0) {
  cat("\n  将分析的基因:\n")
  for (i in 1:min(10, nrow(final_genes))) {  # 只显示前10个
    cat("   ", i, ".", final_genes$gene_symbol[i],
        "(", final_genes$ensembl_id[i], ")\n")
  }
  if (nrow(final_genes) > 10) {
    cat("   ... 还有", nrow(final_genes) - 10, "个基因\n")
  }
}

# ============================================================================
# 4. 运行共定位分析
# ============================================================================

if (nrow(final_genes) == 0) {
  cat("\n  错误: 没有基因可用于共定位分析。\n")
  cat("  请检查筛选条件或数据。\n")
} else {
  cat("\n\n步骤2: 运行共定位分析\n")
  cat("------------------------------------------------------------\n")
  cat("  筛选设置: 保留eQTL p < ", P_THRESHOLD, " 或 GWAS p < ", P_THRESHOLD, " 的SNP\n", sep = "")
  cat("  显著阈值: PP.H4 > ", PP4_THRESHOLD, "\n", sep = "")
  
  # 存储结果
  coloc_results <- data.frame()
  coloc_details_list <- list()  # 存储详细结果
  
  # 对每个基因运行coloc
  for (i in 1:nrow(final_genes)) {
    gene_symbol <- final_genes$gene_symbol[i]
    ensembl_id <- final_genes$ensembl_id[i]
    exposure_id <- if ("exposure_id" %in% colnames(final_genes)) final_genes$exposure_id[i] else paste0("eqtl-a-", ensembl_id)
    
    cat("\n  ", i, "/", nrow(final_genes), "分析基因:", gene_symbol, "(", ensembl_id, ")\n")
    
    # 使用本地eQTL文件和API获取GWAS数据进行coloc分析
    tryCatch({
      # 设置region半径（基因上下游50kb = 50000bp）
      region_radius_bp <- REGION_RADIUS_KB * 1000
      
      cat("    1) 读取本地eQTL数据\n")
      
      # 懒加载eQTL数据（首次运行时读取），使用data.table加速
      if (!exists("eqtl_full_data")) {
        cat("      加载eQTL汇总数据文件(仅需列)...\n")
        # 只读取需要的列以加速加载，保留AssessedAllele和OtherAllele用于等位基因方向校正
        eqtl_full_data <<- fread(EQTL_FILE, select=c("GeneSymbol", "SNP", "SNPPos", "GenePos", "Zscore", "NrSamples", "Pvalue", "AssessedAllele", "OtherAllele", "Beta", "SE"))
        cat("      eQTL数据加载完成，总行数:", nrow(eqtl_full_data), "\n")
      }
      
      # 筛选目标基因的eQTL数据
      gene_eqtl <- eqtl_full_data[GeneSymbol == gene_symbol]
      
      if (nrow(gene_eqtl) == 0) {
        cat("      警告: 基因", gene_symbol, "不在eQTL数据中\n")
        next
      }
      
      cat("      基因eQTL SNP总数:", nrow(gene_eqtl), "\n")
      
      # 筛选基因区域内的SNP (GenePos ± REGION_RADIUS_KB)
      region_eqtl <- gene_eqtl[abs(SNPPos - GenePos) <= region_radius_bp]
      
      if (nrow(region_eqtl) < 10) {
        cat("      警告: 基因区域SNP数太少:", nrow(region_eqtl), "\n")
        next
      }
      
      cat("      基因区域(±", REGION_RADIUS_KB, "kb) eQTL SNP数:", nrow(region_eqtl), "\n")

      # 如果有原始beta和se，直接使用；否则从Zscore计算
      # 原始beta更准确，避免反推引入误差
      if ("Beta" %in% colnames(region_eqtl) && "SE" %in% colnames(region_eqtl)) {
        cat("      使用原始eQTL Beta和SE\n")
        region_eqtl[, beta_eqtl := Beta]
        region_eqtl[, varbeta_eqtl := SE^2]
      } else {
        cat("      从Zscore反推Beta和SE\n")
        # beta = Zscore / sqrt(N)
        # varbeta = 1 / N
        region_eqtl[, beta_eqtl := Zscore / sqrt(NrSamples)]
        region_eqtl[, varbeta_eqtl := 1 / NrSamples]
      }
      
      # 获取这些 SNP 在 GWAS 中的关联数据
      cat("    2) 获取GWAS数据:", OUTCOME_ID, "\n")
      gwas_dat <- ieugwasr::associations(variant = region_eqtl$SNP, id = OUTCOME_ID, proxies = GWAS_USE_PROXIES)
      
      if (is.null(gwas_dat) || nrow(gwas_dat) == 0) {
        cat("      警告: 无法获取GWAS数据\n")
        next
      }
      
      cat("      GWAS数据SNP数:", nrow(gwas_dat), "\n")
      
      # 合并 eQTL 和 GWAS 数据
      merged_dat <- merge(
        region_eqtl[, c("SNP", "beta_eqtl", "varbeta_eqtl", "NrSamples", "Pvalue")],
        gwas_dat[, c("rsid", "beta", "se", "p", "n", "eaf")],
        by.x = "SNP", by.y = "rsid"
      )
      
      # 重命名列
      colnames(merged_dat)[colnames(merged_dat) == "beta"] <- "beta_gwas"
      colnames(merged_dat)[colnames(merged_dat) == "se"] <- "se_gwas"
      colnames(merged_dat)[colnames(merged_dat) == "p"] <- "p_gwas"
      colnames(merged_dat)[colnames(merged_dat) == "n"] <- "n_gwas"

      # 转换eaf为真正的MAF (minor allele frequency)
      # MAF = min(eaf, 1 - eaf)
      if ("eaf" %in% colnames(merged_dat)) {
        merged_dat$MAF <- pmin(merged_dat$eaf, 1 - merged_dat$eaf)
      } else {
        merged_dat$MAF <- 0.5  # 默认值
      }
      
      if (nrow(merged_dat) < 10) {
        cat("      警告: 原始合并SNP数太少:", nrow(merged_dat), "\n")
        next
      }
      
      cat("      原始合并SNP数:", nrow(merged_dat), "\n")

      # ================================================================
      # 基于p值的SNP筛选步骤 (可开关)
      # ================================================================
      if (USE_SNP_FILTER) {
        cat("    3) 基于p值筛选SNP (阈值 = ", P_THRESHOLD, ")\n", sep = "")

        # 筛选：保留在eQTL中显著 或 在GWAS中显著的SNP
        filtered_dat <- merged_dat[merged_dat$Pvalue < P_THRESHOLD | merged_dat$p_gwas < P_THRESHOLD, ]

        # 统计筛选结果
        eqtl_sig <- sum(merged_dat$Pvalue < P_THRESHOLD, na.rm = TRUE)
        gwas_sig <- sum(merged_dat$p_gwas < P_THRESHOLD, na.rm = TRUE)
        both_sig <- sum(merged_dat$Pvalue < P_THRESHOLD & merged_dat$p_gwas < P_THRESHOLD, na.rm = TRUE)

        cat("      筛选统计:\n")
        cat("        - eQTL显著SNP: ", eqtl_sig, "\n", sep = "")
        cat("        - GWAS显著SNP: ", gwas_sig, "\n", sep = "")
        cat("        - 两者都显著SNP: ", both_sig, "\n", sep = "")
        cat("        - 筛选后SNP数: ", nrow(filtered_dat), "\n", sep = "")

        # 计算筛选前后比例
        filter_rate <- round(nrow(filtered_dat) / nrow(merged_dat) * 100, 1)
        cat("      筛选保留比例: ", filter_rate, "%\n", sep = "")
      } else {
        cat("    3) SNP筛选已关闭，使用全部SNP\n")
        filtered_dat <- merged_dat
        filter_rate <- 100
        eqtl_sig <- NA
        gwas_sig <- NA
        both_sig <- NA
      }

      # 保存筛选前的统计信息
      pre_filter_stats <- data.frame(
        gene_symbol = gene_symbol,
        original_snps = nrow(merged_dat),
        filtered_snps = nrow(filtered_dat),
        filter_percent = filter_rate,
        eqtl_sig_snps = eqtl_sig,
        gwas_sig_snps = gwas_sig,
        both_sig_snps = both_sig
      )
      
      # 获取样本量
      n_exp <- unique(na.omit(filtered_dat$NrSamples))
      n_out <- unique(na.omit(filtered_dat$n_gwas))
      
      # 运行coloc分析
      # 使用 beta + varbeta + MAF + N 进行分析
      # eQTL是quantitative traits (quant)，GWAS是case-control (cc)
      # MAF用GWAS的eaf近似
      cat("    4) 运行共定位分析 (使用", nrow(filtered_dat), "个SNP)\n", sep = "")
      
      coloc_result <- coloc::coloc.abf(
        dataset1 = list(
          beta = filtered_dat$beta_eqtl,
          varbeta = filtered_dat$varbeta_eqtl,
          MAF = filtered_dat$MAF,
          snp = filtered_dat$SNP,
          type = "quant",
          N = ifelse(length(n_exp) > 0, n_exp[1], 1000)
        ),
        dataset2 = list(
          beta = filtered_dat$beta_gwas,
          varbeta = filtered_dat$se_gwas^2,
          MAF = filtered_dat$MAF,
          snp = filtered_dat$SNP,
          type = "cc",
          s = CASE_PREVALENCE,  # 使用真实患病率
          N = ifelse(length(n_out) > 0, n_out[1], 20000)
        )
      )
      
      # 提取结果
      summary_result <- data.frame(
        gene_symbol = gene_symbol,
        ensembl_id = ensembl_id,
        exposure_id = exposure_id,
        outcome_id = OUTCOME_ID,
        n_snps = nrow(filtered_dat),
        original_snps = nrow(merged_dat),
        filter_percent = filter_rate,
        PP.H0 = coloc_result$summary["PP.H0.abf"],
        PP.H1 = coloc_result$summary["PP.H1.abf"],
        PP.H2 = coloc_result$summary["PP.H2.abf"],
        PP.H3 = coloc_result$summary["PP.H3.abf"],
        PP.H4 = coloc_result$summary["PP.H4.abf"],
        stringsAsFactors = FALSE
      )
      
      # 保存详细结果
      coloc_details_list[[gene_symbol]] <- list(
        summary = summary_result,
        pre_filter_stats = pre_filter_stats,
        coloc_result = coloc_result,
        filtered_data = filtered_dat
      )
      
      coloc_results <- rbind(coloc_results, summary_result)
      
      cat("      共定位结果: PP.H4 =", round(coloc_result$summary["PP.H4.abf"], 4), "\n")
      
    }, error = function(e) {
      cat("      错误:", e$message, "\n")
    })
  }
  
  # ============================================================================
  # 5. 保存结果
  # ============================================================================
  
  cat("\n\n步骤3: 保存结果\n")
  cat("------------------------------------------------------------\n")
  
  if (nrow(coloc_results) > 0) {
    # 添加筛选结果标记
    coloc_results$is_significant <- coloc_results$PP.H4 > PP4_THRESHOLD
    coloc_results$filter_threshold <- P_THRESHOLD
    
    # 按PP.H4降序排序
    coloc_results <- coloc_results[order(-coloc_results$PP.H4), ]
    
    # 保存完整结果
    output_file <- file.path(COLOC_DIR, "coloc_analysis_results.csv")
    write.csv(coloc_results, output_file, row.names = FALSE)
    cat("  结果已保存:", output_file, "\n")
    
    # 保存详细结果（R数据格式）
    details_file <- file.path(COLOC_DIR, "coloc_analysis_details.RData")
    save(coloc_details_list, file = details_file)
    cat("  详细结果已保存:", details_file, "\n")
    
    # 保存筛选统计摘要
    if (length(coloc_details_list) > 0) {
      filter_stats <- do.call(rbind, lapply(coloc_details_list, function(x) x$pre_filter_stats))
      filter_stats_file <- file.path(COLOC_DIR, "snp_filtering_stats.csv")
      write.csv(filter_stats, filter_stats_file, row.names = FALSE)
      cat("  SNP筛选统计已保存:", filter_stats_file, "\n")
    }
    
    # 打印结果摘要
    cat("\n  共定位分析结果摘要:\n")
    print(coloc_results[, c("gene_symbol", "original_snps", "n_snps", "filter_percent", 
                            "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "is_significant")])
    
    # 筛选显著共定位 (PP.H4 > 0.8)
    cat("\n  显著共定位基因 (PP.H4 > ", PP4_THRESHOLD, "):\n", sep = "")
    significant_coloc <- coloc_results[coloc_results$PP.H4 > PP4_THRESHOLD, ]
    if (nrow(significant_coloc) > 0) {
      print(significant_coloc[, c("gene_symbol", "PP.H4", "n_snps", "filter_percent")])
      cat("\n  显著基因数: ", nrow(significant_coloc), "/", nrow(coloc_results), "\n", sep = "")
    } else {
      cat("    无显著共定位结果 (PP.H4 > ", PP4_THRESHOLD, ")\n", sep = "")
    }
    
    # 筛选中等证据 (PP.H4 > 0.5)
    cat("\n  潜在共定位基因 (PP.H4 > 0.5):\n")
    potential_coloc <- coloc_results[coloc_results$PP.H4 > 0.5, ]
    if (nrow(potential_coloc) > 0) {
      print(potential_coloc[, c("gene_symbol", "PP.H4", "n_snps", "filter_percent")])
      cat("  潜在基因数: ", nrow(potential_coloc), "/", nrow(coloc_results), "\n", sep = "")
    } else {
      cat("    无潜在共定位结果 (PP.H4 > 0.5)\n")
    }
    
    # 筛选结果统计
    cat("\n  SNP筛选效果统计:\n")
    if (exists("filter_stats") && nrow(filter_stats) > 0) {
      cat("    平均保留比例: ", round(mean(filter_stats$filter_percent, na.rm = TRUE), 1), "%\n", sep = "")
      cat("    中位数保留比例: ", round(median(filter_stats$filter_percent, na.rm = TRUE), 1), "%\n", sep = "")
      cat("    范围: ", round(min(filter_stats$filter_percent, na.rm = TRUE), 1), 
          "% - ", round(max(filter_stats$filter_percent, na.rm = TRUE), 1), "%\n", sep = "")
    }
    
  } else {
    cat("  无有效共定位结果\n")
  }
}

cat("\n============================================================\n")
cat("共定位分析完成!\n")
cat("============================================================\n")

# 输出参数设置摘要
cat("\n参数设置摘要:\n")
cat("  - 基因区域半径: ±", REGION_RADIUS_KB, "kb\n", sep = "")
cat("  - SNP筛选阈值: p < ", P_THRESHOLD, "\n", sep = "")
cat("  - 共定位显著阈值: PP.H4 > ", PP4_THRESHOLD, "\n", sep = "")

# 保存qs2与sessionInfo（标准化产物）
if (requireNamespace("qs", quietly = TRUE) && exists("coloc_results")) {
  qs_file <- file.path(COLOC_DIR, "coloc_node.qs2")
  qs::qsave(list(coloc_results = coloc_results, details = if (exists("coloc_details_list")) coloc_details_list else NULL), qs_file)
  cat("  - QS2快照:", qs_file, "\n")
}
log_dir <- file.path(COLOC_DIR, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
session_file <- file.path(log_dir, sprintf("coloc_sessionInfo.%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
writeLines(capture.output(sessionInfo()), session_file)
cat("  - SessionInfo:", session_file, "\n")
