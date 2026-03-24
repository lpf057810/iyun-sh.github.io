# 0. 辅助函数：解析可能带 "<" 前缀的 p 值
safe_parse_pval <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  x <- as.character(x)
  x <- gsub("^\\s*<\\s*", "", x)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0 || is.na(v)) return(NA_real_)
  v[1]
}

# 0b. Ensembl -> Gene Symbol（用于结果展示名）
map_ensembl_to_symbol <- function(ensembl_ids) {
  ids <- unique(sub("\\..*", "", as.character(ensembl_ids)))
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0) {
    return(setNames(character(0), character(0)))
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(setNames(rep(NA_character_, length(ids)), ids))
  }

  mapped <- tryCatch(
    AnnotationDbi::mapIds(
      x = org.Hs.eg.db::org.Hs.eg.db,
      keys = ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ),
    error = function(e) {
      cat("  Ensembl->Symbol mapping failed:", e$message, "\n")
      setNames(rep(NA_character_, length(ids)), ids)
    }
  )

  mapped_chr <- as.character(mapped)
  names(mapped_chr) <- names(mapped)
  mapped_chr
}

resolve_display_ids <- function(ids) {
  raw_ids <- as.character(ids)
  clean_ids <- sub("\\..*", "", raw_ids)
  is_ensembl <- grepl("^ENSG[0-9]+$", toupper(clean_ids))
  out <- raw_ids
  if (any(is_ensembl)) {
    symbol_map <- map_ensembl_to_symbol(clean_ids[is_ensembl])
    mapped <- unname(symbol_map[clean_ids[is_ensembl]])
    mapped <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, raw_ids[is_ensembl])
    out[is_ensembl] <- mapped
  }
  out
}

save_qs2_snapshot <- function(object, output_path) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    cat("Warning: package 'qs' not installed, skip qs2 snapshot.\n")
    return(invisible(NULL))
  }
  qs::qsave(object, output_path)
  cat("✓ QS2 snapshot saved:", output_path, "\n")
  invisible(output_path)
}

save_session_info_snapshot <- function(log_dir, dataset_name, run_tag) {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  session_path <- file.path(log_dir, sprintf("%s_sessionInfo.%s.txt", dataset_name, run_tag))
  writeLines(capture.output(sessionInfo()), session_path)
  cat("✓ Session info saved:", session_path, "\n")
  invisible(session_path)
}

# 0. 读取环境变量
# 如果有测试数据需求，可以设置 WORK_DIR 环境变量
work_dir <- Sys.getenv("WORK_DIR", getwd())
setwd(work_dir)

# Token loading handled after option parsing
# ============================================================================
# pQTL孟德尔随机化分析流程 - 修正版（使用rsID）
# ============================================================================

# 1. 加载必要的R包
library(optparse)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(dplyr)
library(data.table)
library(MRPRESSO)
library(ggsci)
library(showtext)
showtext_auto()
showtext_opts(dpi = 300)
font_add("Liberation Sans", regular = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", bold = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", italic = "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf", bolditalic = "/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf")
library(clusterProfiler)

# Config reading function (sync with MR_eqtl.R)
read_mr_config <- function() {
  config_paths <- c(
    "scripts/config/MR.config.ini",
    "../scripts/config/MR.config.ini",
    "/media/desk16/share/secure/MR/scripts/config/MR.config.ini"
  )

  config_file <- NULL
  for (path in config_paths) {
    if (file.exists(path)) {
      config_file <- path
      break
    }
  }

  if (!is.null(config_file) && requireNamespace("configr", quietly = TRUE)) {
    cat("[Config] Loading config from:", config_file, "\n")
    return(configr::read.config(config_file))
  }

  cat("[Config] No config file found or package 'configr' unavailable, using defaults\n")
  NULL
}

# Timestamp function
timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

# =============================================================================
# Token Management Functions
# =============================================================================

#' Test if token is valid by making a lightweight API call
test_token_validity <- function(token) {
  tryCatch({
    result <- ieugwasr::gwasinfo("ieu-a-2")
    if (!is.null(result)) {
      return(list(valid = TRUE, message = "Token is valid"))
    }
    return(list(valid = FALSE, message = "API returned NULL"))
  }, error = function(e) {
    err_msg <- e$message
    if (grepl("token|unauthorized|401|Please generate|forbidden", err_msg, ignore.case = TRUE)) {
      return(list(valid = FALSE, message = "TOKEN EXPIRED - Please update token"))
    }
    return(list(valid = FALSE, message = paste("API error:", err_msg)))
  })
}

#' Check if error is due to token expiry
is_token_expired_error <- function(error_msg) {
  grepl("token|unauthorized|401|Please generate|forbidden", error_msg, ignore.case = TRUE)
}

#' Retry wrapper with token expiry detection
retry_with_token_check <- function(expr, max_retries = 3, on_token_expired = NULL) {
  attempt <- 0
  last_error <- NULL

  while (attempt < max_retries) {
    attempt <- attempt + 1
    result <- tryCatch({
      expr
    }, error = function(e) {
      last_error <<- e
      NULL
    })

    if (!is.null(result)) {
      return(result)
    }

    if (!is.null(last_error) && is_token_expired_error(last_error$message)) {
      cat("[RETRY] Token expired! No more retries for expired token.\n")
      if (!is.null(on_token_expired)) {
        on_token_expired()
      }
      stop(last_error)
    }

    if (attempt < max_retries) {
      wait_time <- 2^(attempt - 1)
      cat("[RETRY] Attempt", attempt, "failed:", last_error$message, "- Retrying in", wait_time, "s...\n")
      Sys.sleep(wait_time)
    }
  }

  cat("[ERROR] All", max_retries, "attempts failed. Last error:", last_error$message, "\n")
  return(NULL)
}

#' Load token with priority: command line > env > .Renviron
load_token <- function(cmd_token = NULL) {
  if (!is.null(cmd_token) && nchar(cmd_token) > 0) {
    cat("[Token] Using token from command line\n")
    Sys.setenv(OPENGWAS_JWT = cmd_token)
    return(cmd_token)
  }

  token <- Sys.getenv("OPENGWAS_JWT")
  if (nchar(token) > 0) {
    cat("[Token] Using token from environment variable\n")
    return(token)
  }

  global_renv <- file.path(Sys.getenv("HOME"), ".Renviron")
  if (file.exists(global_renv)) {
    readRenviron(global_renv)
    token <- Sys.getenv("OPENGWAS_JWT")
    if (nchar(token) > 0) {
      cat("[Token] Using token from ~/.Renviron\n")
      return(token)
    }
  }

  cat("[Token] WARNING: No token found!\n")
  return(NULL)
}

# Safe extractor for first element/column (sync with MR_eqtl.R)
extract_safe <- function(df, col_name) {
  if (is.null(df) || nrow(df) == 0 || is.null(col_name) || !(col_name %in% colnames(df))) {
    return(NA)
  }
  val <- df[[col_name]][1]
  if (is.null(val) || length(val) == 0) return(NA)
  val
}

# Command line argument parsing
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input gene list file (CSV)"),
  make_option(c("-o", "--output"), type = "character", default = "./output",
              help = "Output directory [default %default]"),
  make_option(c("-p", "--pval"), type = "character", default = "5e-8",
              help = "P-value threshold [default %default]"),
  make_option(c("-k", "--kb"), type = "character", default = "10000",
              help = "Distance in kb for clumping [default %default]"),
  make_option(c("-r", "--r2"), type = "character", default = "0.001",
              help = "LD r2 threshold [default %default]"),
  make_option(c("-f", "--fstat"), type = "character", default = "10",
              help = "F-statistic threshold [default %default]"),
  make_option(c("-e", "--eaf"), type = "character", default = "0.01",
              help = "EAF threshold [default %default]"),
  make_option(c("-d", "--outcome"), type = "character", default = "ebi-a-GCST90091033",
              help = "Outcome GWAS ID [default %default]"),
  make_option(c("-q", "--pqtl-dir"), type = "character", default = "/media/desk16/iyunlpf/NAF/pqtl/data/cis_clump",
              help = "pQTL data directory [default %default]", dest = "pqtl_dir"),
  make_option(c("-t", "--token"), type = "character", default = NULL,
              help = "OpenGWAS JWT token (optional, overrides env/file token)"),
  make_option(c("--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run [default auto-generated]")
)

# Read config file if exists (kept for parity with MR_eqtl.R)
mr_config <- read_mr_config()

opt <- parse_args(OptionParser(option_list = option_list))

# Load token with priority
token <- load_token(opt$token)
token_loaded <- !is.null(token) && nchar(token) > 0

# Print header
cat("\n")
cat("================================================================================\n")
cat("  MR pQTL Analysis\n")
cat("  Started:", timestamp(), "\n")
cat("================================================================================\n\n")

# Set Arial as default font
par(family = "Liberation Sans")

# Test token validity at startup
cat("[Token Check] Validating OpenGWAS token...\n")
token_test <- test_token_validity(token)
if (token_test$valid) {
  cat("[Token Check] OK - Token is valid\n")
} else {
  cat("[Token Check] FAILED -", token_test$message, "\n")
  cat("[Token Check] To update token, run: ./update_opengwas_token.sh <new_token>\n")
}

base_dir <- opt$output

# Create output directory early
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# Save analysis config to result directory
config_file <- file.path(base_dir, "analysis_config.ini")
config_lines <- c(
  "[Analysis]",
  paste0("type = pqtl"),
  paste0("timestamp = ", opt$timestamp),
  "",
  "[Parameters]",
  paste0("pval_threshold = ", opt$pval),
  paste0("clump_kb = ", opt$kb),
  paste0("clump_r2 = ", opt$r2),
  paste0("fstat_threshold = ", opt$fstat),
  paste0("eaf_threshold = ", opt$eaf),
  paste0("outcome_gwas = ", opt$outcome),
  paste0("min_snps = 3")
)
writeLines(config_lines, config_file)
cat("[Config] Parameters saved to:", config_file, "\n")

# ============================================================================
# 步骤1: 读取候选基因列表（智能列识别）
# ============================================================================
cat("STEP 1: Loading candidate gene list...\n")
if (is.null(opt$input)) {
  stop("Error: Input file is required. Use -i or --input to specify the gene list file.")
}
genes <- read.csv(opt$input, check.names = FALSE)
cat("  Original dimensions:", nrow(genes), "rows,", ncol(genes), "columns\n")
cat("  Original column names:", paste(colnames(genes), collapse = ", "), "\n")
cat("  First few rows:\n")
print(head(genes))

# 智能识别基因列
gene_col <- NULL
id_col   <- NULL
gene_col_names <- c("gene", "symbol", "hgnc_symbol", "gene_symbol",
                    "entrez", "entrezid", "entrez_gene_id",
                    "ensembl", "ensembl_id", "ensembl_gene_id")
id_col_names <- c("id", "snp_id", "snpid", "rsid", "variant_id")

for (nm in gene_col_names) {
  hit <- which(tolower(colnames(genes)) == nm)
  if (length(hit) > 0) { gene_col <- colnames(genes)[hit[1]]; break }
}
for (nm in id_col_names) {
  hit <- which(tolower(colnames(genes)) == nm)
  if (length(hit) > 0) { id_col <- colnames(genes)[hit[1]]; break }
}

if (is.null(gene_col)) {
  # 备用：用内容验证排除序号列
  for (i in seq_len(ncol(genes))) {
    col_vals <- as.character(genes[[i]])
    col_vals <- na.omit(col_vals)
    if (length(col_vals) == 0) next
    is_numeric <- all(grepl("^\\d+(\\.\\d+)?$", col_vals))
    is_too_short <- all(nchar(col_vals) <= 2)
    is_gene_like <- any(grepl("^[A-Z][A-Z0-9]{2,}$", toupper(col_vals)))
    if (!is_numeric && !is_too_short && is_gene_like) {
      gene_col <- colnames(genes)[i]
      cat("  ✓ Gene column inferred by content validation: '", gene_col, "'\n", sep = "")
      break
    }
  }
  if (is.null(gene_col)) {
    gene_col <- colnames(genes)[1]
    if (ncol(genes) > 1) id_col <- colnames(genes)[2]
    cat("  [WARN] No recognized gene column found. Using first column: '", gene_col, "'\n", sep = "")
  }
} else {
  cat("  ✓ Gene column detected: '", gene_col, "'\n", sep = "")
}
if (!is.null(id_col)) {
  cat("  ✓ ID column detected: '", id_col, "'\n", sep = "")
}

# 标准化列名
if (gene_col != "gene") {
  colnames(genes)[which(colnames(genes) == gene_col)] <- "gene"
  gene_col <- "gene"
}
if (!is.null(id_col) && id_col != "id") {
  colnames(genes)[which(colnames(genes) == id_col)] <- "id"
  id_col <- "id"
} else if (is.null(id_col) && ncol(genes) > 1) {
  genes$id <- 1:nrow(genes)
  cat("  ✓ Added auto-generated id column\n")
}

candidate_genes <- unique(na.omit(genes$gene))
cat("  Total candidate genes:", length(candidate_genes), "\n")
cat("  Sample genes:", paste(head(candidate_genes, min(5, length(candidate_genes))), collapse = ", "), "\n")

# ============================================================================
# 步骤2: 扫描pQTL文件目录并匹配基因
# ============================================================================
cat("\nSTEP 2: Scanning pQTL directory and matching genes...\n")

# pQTL目录路径
pqtl <- opt$pqtl_dir
cat("  pQTL directory:", pqtl, "\n")

if (!dir.exists(pqtl)) {
  stop("ERROR: pQTL directory not found: ", pqtl)
}

# 列出所有txt文件
pqtl_files <- list.files(pqtl, pattern = "\\.txt$", full.names = TRUE)
cat("  Found", length(pqtl_files), "pQTL files\n")

# 从文件名中提取基因名
extract_gene_from_filename <- function(filename) {
  basename <- basename(filename)
  parts <- strsplit(basename, "_")[[1]]
  return(parts[1])
}

# 创建文件信息数据框
file_info <- data.frame(
  file_path = pqtl_files,
  file_name = basename(pqtl_files),
  gene_from_filename = sapply(pqtl_files, extract_gene_from_filename),
  stringsAsFactors = FALSE
)

# 大小写不敏感匹配
cat("\n  Performing case-insensitive matching...\n")
matched_files <- c()
matched_genes <- c()

# 转换为小写以便不区分大小写匹配
candidate_lower <- tolower(candidate_genes)
file_genes_lower <- tolower(file_info$gene_from_filename)

for (i in 1:length(candidate_lower)) {
  candidate_gene_lower <- candidate_lower[i]
  candidate_gene_original <- candidate_genes[i]
  
  # 查找匹配的文件
  match_indices <- which(file_genes_lower == candidate_gene_lower)
  
  if (length(match_indices) > 0) {
    # 找到匹配，记录第一个匹配的文件
    matched_files <- c(matched_files, file_info$file_path[match_indices[1]])
    matched_genes <- c(matched_genes, candidate_gene_original)
    cat("    Found match:", candidate_gene_original, "->", 
        file_info$gene_from_filename[match_indices[1]], "\n")
  }
}

cat("\n  Total matches found:", length(matched_genes), "\n")

if (length(matched_genes) == 0) {
  stop("No overlapping genes found. Please check gene naming conventions.")
}

# 分析键 -> 展示名（仅用于输出，不影响内部分析逻辑）
analysis_display_map <- setNames(resolve_display_ids(matched_genes), matched_genes)

# 创建匹配信息数据框
matched_info <- data.frame(
  candidate_gene = matched_genes,
  display_gene = unname(analysis_display_map[matched_genes]),
  pqtl_file_gene = file_info$gene_from_filename[match(basename(matched_files), file_info$file_name)],
  file_name = basename(matched_files),
  file_path = matched_files,
  stringsAsFactors = FALSE
)

write.csv(matched_info,
          file = file.path(base_dir, "08.pQTL_ID_Mapping.csv"),
          row.names = FALSE)
cat("  ✓ Matching information saved to '08.pQTL_ID_Mapping.csv'\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 1-2 COMPLETE. Found", length(matched_genes), "pQTL exposures.\n")
cat(rep("=", 60), "\n\n")

# ============================================================================
# 步骤3: 构建暴露列表并运行MR分析（使用rsID）
# ============================================================================

# 准备结局ID
outcome_id <- opt$outcome

# 工具变量筛选参数
clump_params <- list(p1 = as.numeric(opt$pval), r2 = as.numeric(opt$r2), kb = as.numeric(opt$kb))
f_stat_threshold <- as.numeric(opt$fstat)
eaf_threshold <- as.numeric(opt$eaf)

# 定义输出目录结构

tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")

# 创建目录
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "01_scatter"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "02_forest"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "03_funnel"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "04_leaveoneout"), showWarnings = FALSE, recursive = TRUE)

# 创建主日志文件
log_file <- file.path(base_dir, "MR_pQTL_Analysis_Log.txt")
sink(log_file, split = TRUE)
cat("=== pQTL Mendelian Randomization Analysis Log ===\n")
cat("Start Time:", format(Sys.time()), "\n")
cat("Analysis Type: pQTL (Protein Quantitative Trait Loci)\n")
cat("Total pQTL exposures to process:", length(matched_genes), "\n")
cat("Outcomes:", paste(outcome_id, collapse = ", "), "\n")
cat("Using rsIDs from 'rsids2' column for SNP matching\n")
cat("=================================================\n\n")

# 初始化结果列表
results_summary_pqtl <- data.frame()
all_analyses <- list()

# 主分析循环
for (i in seq_along(matched_files)) {
  file_path <- matched_files[i]
  gene_name <- matched_genes[i]
  file_basename <- basename(file_path)
  safe_name <- gsub("[^[:alnum:]]", "_", file_basename)
  
  cat("\n", rep("=", 60), "\n", sep="")
  cat("ANALYSIS [", i, "/", length(matched_files), "]: ", gene_name, "\n", sep="")
  cat("File: ", file_basename, "\n", sep="")
  cat(rep("=", 60), "\n\n", sep="")
  
  # 3.1 读取并处理pQTL工具变量文件（使用rsID）
  cat("STEP 1: Reading and processing pQTL file (using rsIDs)...\n")
  
  exposure_dat <- tryCatch({
    # 读取文件
    pqtl_data <- fread(file_path)
    
    # 验证必要的列是否存在
    required_cols <- c("ID", "ALLELE1", "ALLELE0", "BETA", "SE", "LOG10P", "A1FREQ", "rsids2")
    missing_cols <- setdiff(required_cols, colnames(pqtl_data))
    
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    cat("  File loaded successfully:", nrow(pqtl_data), "rows\n")
    
    # 计算p值（从LOG10P转换）
    pqtl_data$pval <- 10^(-pqtl_data$LOG10P)
    
    # 筛选显著SNP (p < 5e-8)
    pqtl_data <- pqtl_data[pqtl_data$pval < clump_params$p1, ]
    
    if (nrow(pqtl_data) == 0) {
      stop("No significant SNPs found (p < ", clump_params$p1, ")")
    }
    cat("  Significant SNPs (p < 5e-8):", nrow(pqtl_data), "\n")
    
    # 从rsids2列提取第一个rsID（有些可能有多个用逗号分隔）
    cat("  Extracting rsIDs from 'rsids2' column...\n")
    pqtl_data$rsid <- sapply(strsplit(as.character(pqtl_data$rsids2), ","), function(x) {
      if (length(x) > 0 && !is.na(x[1]) && x[1] != "" && x[1] != "NA") {
        return(trimws(x[1]))  # 取第一个rsID并去除空格
      } else {
        return(NA)
      }
    })
    
    # 移除没有rsID的行
    pqtl_data <- pqtl_data[!is.na(pqtl_data$rsid) & pqtl_data$rsid != "", ]
    cat("  SNPs with valid rsIDs:", nrow(pqtl_data), "\n")
    
    if (nrow(pqtl_data) == 0) {
      stop("No SNPs with valid rsIDs found")
    }
    
    # 检查rsID格式
    cat("  Sample rsIDs:", paste(head(pqtl_data$rsid, 3), collapse = ", "), "\n")
    
    # 转换为TwoSampleMR格式，使用rsid作为SNP标识符
    exp_dat <- data.frame(
      SNP = pqtl_data$rsid,  # 关键修改：使用rsid而不是ID
      effect_allele.exposure = pqtl_data$ALLELE1,
      other_allele.exposure = pqtl_data$ALLELE0,
      beta.exposure = pqtl_data$BETA,
      se.exposure = pqtl_data$SE,
      pval.exposure = pqtl_data$pval,
      eaf.exposure = pqtl_data$A1FREQ,
      exposure = gene_name,
      id.exposure = gene_name,
      stringsAsFactors = FALSE
    )
    
    # 添加其他可用信息
    if ("N" %in% colnames(pqtl_data)) {
      exp_dat$samplesize.exposure <- pqtl_data$N[1]
    }
    
    # 质量控制
    exp_dat <- exp_dat %>% 
      filter(eaf.exposure > eaf_threshold & eaf.exposure < (1 - eaf_threshold))
    cat("  After EAF filter (0.01-0.99):", nrow(exp_dat), "IVs\n")
    
    # 计算F统计量
    if (!"samplesize.exposure" %in% colnames(exp_dat)) {
      exp_dat$samplesize.exposure <- 33995
    }
    
    exp_dat$f_statistic <- (exp_dat$beta.exposure^2) / (exp_dat$se.exposure^2)
    exp_dat <- exp_dat %>% filter(f_statistic > f_stat_threshold)
    cat("  After F-statistic filter (>", f_stat_threshold, "):", nrow(exp_dat), "IVs\n")
    
    if (nrow(exp_dat) < 3) {
      cat("  WARNING: Fewer than 3 IVs remaining\n")
    }
    
    # 去重（基于SNP）
    exp_dat <- exp_dat[!duplicated(exp_dat$SNP), ]
    cat("  After removing duplicates:", nrow(exp_dat), "unique IVs\n")
    
    exp_dat
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
    results_summary_pqtl <- rbind(results_summary_pqtl, 
                                  data.frame(exposure = gene_name, status = "Failed_IV_Extraction"))
    next
  }
  
  # 3.2 提取结局数据
  cat("\nSTEP 2: Extracting outcome data...\n")
  outcome_dat <- tryCatch({
    cat("  Extracting data for", nrow(exposure_dat), "SNP(s) from outcome(s):", paste(outcome_id, collapse = ", "), "\n")

    # 尝试从第一个结局开始提取 (with retry)
    out_dat <- tryCatch({
      retry_with_token_check(
        extract_outcome_data(
          snps = exposure_dat$SNP,
          outcomes = outcome_id
        ),
        max_retries = 3,
        on_token_expired = function() {
          cat("  FATAL: Token expired during outcome extraction!\n")
        }
      )
    }, error = function(e1) {
      cat("  First attempt failed:", e1$message, "\n")
      cat("  Trying to extract outcomes separately...\n")
      
      # 如果合并提取失败，尝试分别提取
      out_list <- list()
      for (outcome in outcome_id) {
        cat("  Extracting from", outcome, "...\n")
        out_single <- tryCatch({
          retry_with_token_check(
            extract_outcome_data(
              snps = exposure_dat$SNP,
              outcomes = outcome
            ),
            max_retries = 2
          )
        }, error = function(e2) {
          cat("    Failed to extract from", outcome, ":", e2$message, "\n")
          return(NULL)
        })
        
        if (!is.null(out_single) && nrow(out_single) > 0) {
          out_list[[outcome]] <- out_single
        }
      }
      
      if (length(out_list) == 0) {
        stop("No outcome data could be extracted from any outcome")
      }
      
      # 合并所有结局数据
      do.call(rbind, out_list)
    })
    
    if (is.null(out_dat) || nrow(out_dat) == 0) {
      stop("No matching SNPs in outcome")
    }
    cat("  Found", nrow(out_dat), "SNPs in outcome(s)\n")
    out_dat
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    cat("  This may be due to:\n")
    cat("  1. rsIDs not found in OpenGWAS database\n")
    cat("  2. Network connection issues\n")
    cat("  3. Outcome GWAS IDs may be incorrect\n")
    cat("  4. Try checking if outcome IDs are valid in OpenGWAS\n")
    return(NULL)
  })
  
  if (is.null(outcome_dat)) {
    results_summary_pqtl <- rbind(results_summary_pqtl,
                                  data.frame(exposure = gene_name, status = "Failed_Outcome_Extraction"))
    next
  }
  
  # 3.3 协调数据
  cat("\nSTEP 3: Harmonising data...\n")
  dat <- tryCatch({
    h_dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
    if (nrow(h_dat) < 2) {
      stop("Not enough SNPs after harmonisation (need at least 2)")
    }
    cat("  Successfully harmonised", nrow(h_dat), "SNPs\n")
    h_dat
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(dat)) {
    results_summary_pqtl <- rbind(results_summary_pqtl,
                                  data.frame(exposure = gene_name, status = "Failed_Harmonisation"))
    next
  }

  # 3.3.1 剔除与结局显著相关的SNP（排除水平多效性，p < 1e-5）
  if ("pval.outcome" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$pval.outcome >= 1e-5, ]
    n_after <- nrow(dat)
    if (n_before > n_after) {
      cat("  Removed", n_before - n_after, "SNPs associated with outcome (p < 1e-5)\n")
    }
  }

  # 检查过滤后是否仍有足够SNP
  if (nrow(dat) < 3) {
    cat("  WARNING: Fewer than 3 SNPs after outcome-associated filter\n")
  }

  # 3.4 异质性检验 (先运行以决定使用哪种IVW模型)
  cat("\nSTEP 4: Testing heterogeneity to select MR method...\n")
  het_res <- tryCatch({
    het <- mr_heterogeneity(dat)
    if (!is.null(het)) {
      het$exposure <- gene_name
      cat("  Heterogeneity methods:", paste(unique(het$method), collapse = ", "), "\n")
    }
    het
  }, error = function(e) {
    cat("  Heterogeneity test failed:", e$message, "\n")
    return(NULL)
  })

  # 根据异质性选择MR方法
  # 如果Q_pval < 0.05，存在显著异质性，使用随机效应模型
  # 否则使用固定效应模型
  het_qvalue <- NA
  if (!is.null(het_res) && nrow(het_res) > 0) {
    cat("  Heterogeneity results:\n")
    print(het_res[, c("method", "Q", "Q_pval")])
    # 查找包含"Inverse variance"的方法
    ivw_idx <- grep("Inverse variance", het_res$method, ignore.case = TRUE)
    if (length(ivw_idx) > 0) {
      het_qvalue <- het_res$Q_pval[ivw_idx[1]]
      cat("  Selected IVW Q_pval:", het_qvalue, "\n")
    } else {
      cat("  WARNING: Could not find IVW method in heterogeneity results\n")
    }
  }

  # 根据异质性选择方法列表
  if (!is.na(het_qvalue) && het_qvalue < 0.05) {
    cat("  Heterogeneity detected (Q_pval =", format.pval(het_qvalue, digits = 3),
        ") → Using RANDOM EFFECTS IVW\n")
    method_list <- c(
      "mr_ivw_mre",      # Random effects
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_simple_mode"
    )
  } else {
    cat("  No significant heterogeneity (Q_pval =",
        ifelse(is.na(het_qvalue), "NA", format.pval(het_qvalue, digits = 3)),
        ") → Using FIXED EFFECTS IVW\n")
    method_list <- c(
      "mr_ivw_fe",       # Fixed effects
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_simple_mode"
    )
  }

  # 3.5 运行主MR分析
  cat("\nSTEP 5: Running MR analysis...\n")
  mr_res <- tryCatch({
    res <- mr(dat, method_list = method_list)

    if (nrow(res) == 0) {
      stop("MR analysis returned no results")
    }

    cat("  Methods used:", paste(unique(res$method), collapse = ", "), "\n")

    # 计算OR
    res_or <- generate_odds_ratios(res)
    res_or$exposure <- gene_name
    res_or$outcome <- outcome_id
    res_or$n_snps <- nrow(dat)
    res_or$pqtl_file <- file_basename

    res_or
  }, error = function(e) {
    cat("  ERROR in MR analysis:", e$message, "\n")
    return(NULL)
  })

  if (is.null(mr_res)) {
    results_summary_pqtl <- rbind(results_summary_pqtl,
                                  data.frame(exposure = gene_name, status = "Failed_MR_Analysis"))
    next
  }

  # 3.6 多效性检验
  cat("\nSTEP 6: Pleiotropy test...\n")
  pleio_res <- tryCatch({
    pleio <- mr_pleiotropy_test(dat)
    if (!is.null(pleio)) {
      pleio$exposure <- gene_name
    }
    pleio
  }, error = function(e) {
    cat("  Pleiotropy test failed:", e$message, "\n")
    return(NULL)
  })
  
  # 3.7 Steiger方向性分析
  cat("\nSTEP 7: Steiger directionality test...\n")
  steiger_res <- tryCatch({
    if (!"samplesize.exposure" %in% colnames(dat)) {
      dat$samplesize.exposure <- 33995
    }
    if (!"samplesize.outcome" %in% colnames(dat)) {
      dat$samplesize.outcome <- 10000
    }
    
    has_r_exposure <- "r.exposure" %in% colnames(dat)
    has_r_outcome <- "r.outcome" %in% colnames(dat)
    
    if (!has_r_exposure || !has_r_outcome) {
      cat("  Adding r values for Steiger test...\n")
      
      ncase <- 10000
      ncontrol <- 10000
      prevalence <- 0.1
      
      for (j in 1:nrow(dat)) {
        if (!has_r_exposure) {
          beta <- dat$beta.exposure[j]
          eaf <- ifelse(is.na(dat$eaf.exposure[j]), 0.5, dat$eaf.exposure[j])
          if (!is.na(beta) && !is.na(eaf)) {
            dat$r.exposure[j] <- get_r_from_lor(
              lor = beta,
              af = eaf,
              ncase = ncase,
              ncontrol = ncontrol,
              prevalence = prevalence
            )
          } else {
            dat$r.exposure[j] <- 0.1
          }
        }
        
        if (!has_r_outcome) {
          beta <- dat$beta.outcome[j]
          eaf <- ifelse(is.na(dat$eaf.outcome[j]), 0.5, dat$eaf.outcome[j])
          if (!is.na(beta) && !is.na(eaf)) {
            dat$r.outcome[j] <- get_r_from_lor(
              lor = beta,
              af = eaf,
              ncase = ncase,
              ncontrol = ncontrol,
              prevalence = prevalence
            )
          } else {
            dat$r.outcome[j] <- 0.1
          }
        }
      }
    }
    
    steiger <- directionality_test(dat)
    if (!is.null(steiger) && nrow(steiger) > 0) {
      steiger$exposure <- gene_name
      cat("  Steiger test completed successfully\n")
    } else {
      cat("  WARNING: Steiger test returned empty result\n")
      steiger <- data.frame(
        exposure = gene_name,
        snp_r2.exposure = NA,
        snp_r2.outcome = NA,
        correct_causal_direction = FALSE,
        steiger_pval = NA
      )
    }
    steiger
  }, error = function(e) {
    cat("  Steiger test failed:", e$message, "\n")
    data.frame(
      exposure = gene_name,
      snp_r2.exposure = NA,
      snp_r2.outcome = NA,
      correct_causal_direction = FALSE,
      steiger_pval = NA
    )
  })
  
  # 3.8 MR-PRESSO分析
  cat("\nSTEP 8: MR-PRESSO analysis...\n")
  presso_res <- tryCatch({
    if (nrow(dat) >= 4) {
      presso <- MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        data = dat,
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )

      # 提取结果
      presso_summary <- data.frame(
        exposure = gene_name,
        outcome = outcome_id,
        global_test_p = NA,
        n_outliers = 0,
        outlier_snps = NA,
        outlier_corrected_estimate = NA,
        outlier_corrected_se = NA,
        outlier_corrected_p = NA
      )

      if (!is.null(presso) && !is.null(presso$`MR-PRESSO results`)) {
        gt <- presso$`MR-PRESSO results`$`Global Test`
        if (!is.null(gt) && length(gt) > 0) {
          if (!is.null(gt$Pvalue)) {
            presso_summary$global_test_p <- safe_parse_pval(gt$Pvalue)
          }
        }

        if (!is.null(presso$`MR-PRESSO results`$`Outlier Test`)) {
          outlier_test <- presso$`MR-PRESSO results`$`Outlier Test`
          if (!is.null(outlier_test) && nrow(outlier_test) > 0) {
            outlier_names <- outlier_test$Outliers
            if (!is.null(outlier_names) && length(outlier_names) > 0) {
              presso_summary$n_outliers <- nrow(outlier_test)
              presso_summary$outlier_snps <- paste(outlier_names, collapse = ", ")
            }
          }
        }

        if (!is.null(presso$`Main MR results`) && nrow(presso$`Main MR results`) >= 2) {
          presso_summary$outlier_corrected_estimate <- presso$`Main MR results`[2, "Causal Estimate"]
          presso_summary$outlier_corrected_se <- presso$`Main MR results`[2, "Sd"]
          presso_summary$outlier_corrected_p <- safe_parse_pval(presso$`Main MR results`[2, "P-value"])
        }
      }
      
      presso_summary
    } else {
      cat("  SKIP: Not enough SNPs for MR-PRESSO (need at least 4)\n")
      data.frame(
        exposure = gene_name,
        outcome = outcome_id,
        global_test_p = NA,
        n_outliers = 0,
        outlier_corrected_estimate = NA,
        outlier_corrected_se = NA,
        outlier_corrected_p = NA
      )
    }
  }, error = function(e) {
    cat("  MR-PRESSO failed:", e$message, "\n")
    return(NULL)
  })

  # 3.9 可视化
  cat("\nSTEP 9: Generating plots...\n")

  # 定义通用主题 (Arial字体)
  theme_mr <- function(base_size = 12, legend.position = "right") {
    theme_bw(base_size = base_size) +
      theme(
        text = element_text(face = "bold", family = "Liberation Sans"),
        axis.title = element_text(face = "bold", size = base_size + 2, family = "Liberation Sans"),
        axis.text = element_text(size = base_size, family = "Liberation Sans"),
        legend.title = element_text(face = "bold", size = base_size, family = "Liberation Sans"),
        legend.text = element_text(size = base_size - 1, family = "Liberation Sans"),
        legend.position = legend.position,
        plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0.5, family = "Liberation Sans"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95")
      )
  }

  # 散点图 - legend在顶部
  tryCatch({
    p <- mr_scatter_plot(mr_res, dat)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12, legend.position = "top") +
        ggtitle("MR Scatter Plot")

      num_prefix <- sprintf("%02d", i)
      showtext_begin()
      ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", gene_name, ".scatter.png")),
             plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
      showtext_end()
      showtext_begin()
      ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", gene_name, ".scatter.pdf")),
             plot = p1, width = 10, height = 9, device = "pdf")
      showtext_end()
      cat("  ✓ Scatter plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Scatter plot failed:", e$message, "\n")
  })

  # 森林图 - 增加宽度
  tryCatch({
    res_single <- mr_singlesnp(dat)
    p <- mr_forest_plot(res_single)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12) +
        ggtitle("MR Forest Plot")

      num_prefix <- sprintf("%02d", i)
      showtext_begin()
      ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_name, ".forest.png")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.5), dpi = 300, bg = "white")
      showtext_end()
      showtext_begin()
      ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_name, ".forest.pdf")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.5), device = "pdf")
      showtext_end()
      # 保存rsID.csv
      if (!is.null(res_single) && nrow(res_single) > 0 && "SNP" %in% colnames(res_single)) {
        rsid_df <- data.frame(SNP = res_single$SNP, stringsAsFactors = FALSE)
        write.csv(rsid_df, file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_name, ".rsID.csv")),
                  row.names = FALSE)
      }
      cat("  ✓ Forest plot saved (PNG + PDF + rsID.csv)\n")
    }
  }, error = function(e) {
    cat("  Forest plot failed:", e$message, "\n")
  })

  # 漏斗图 - legend在顶部
  tryCatch({
    res_single <- mr_singlesnp(dat)
    p <- mr_funnel_plot(res_single)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12, legend.position = "top") +
        ggtitle("MR Funnel Plot")

      num_prefix <- sprintf("%02d", i)
      showtext_begin()
      ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", gene_name, ".funnel.png")),
             plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
      showtext_end()
      showtext_begin()
      ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", gene_name, ".funnel.pdf")),
             plot = p1, width = 10, height = 9, device = "pdf")
      showtext_end()
      cat("  ✓ Funnel plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Funnel plot failed:", e$message, "\n")
  })

  # 留一法图 - 增加宽度
  tryCatch({
    res_loo <- mr_leaveoneout(dat)
    p <- mr_leaveoneout_plot(res_loo)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12) +
        ggtitle("MR Leave-One-Out Plot")

      num_prefix <- sprintf("%02d", i)
      showtext_begin()
      ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", gene_name, ".leaveoneout.png")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.4), dpi = 300, bg = "white")
      showtext_end()
      showtext_begin()
      ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", gene_name, ".leaveoneout.pdf")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.4), device = "pdf")
      showtext_end()
      cat("  ✓ Leave-one-out plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Leave-one-out plot failed:", e$message, "\n")
  })

  # 3.10 保存结果
  cat("\nSTEP 10: Saving results...\n")
  
  # 保存协调数据
  write.csv(dat, file.path(tables_dir, paste0(gene_name, "_harmonised.csv")),
            row.names = FALSE)

  # 保存MR结果
  write.csv(mr_res, file.path(tables_dir, paste0(gene_name, "_mr_results.csv")),
            row.names = FALSE)

  # 保存异质性结果
  if (!is.null(het_res)) {
    write.csv(het_res, file.path(tables_dir, paste0(gene_name, "_heterogeneity.csv")),
              row.names = FALSE)
  }

  # 保存多效性结果
  if (!is.null(pleio_res)) {
    write.csv(pleio_res, file.path(tables_dir, paste0(gene_name, "_pleiotropy.csv")),
              row.names = FALSE)
  }

  # 保存Steiger结果
  if (!is.null(steiger_res)) {
    write.csv(steiger_res, file.path(tables_dir, paste0(gene_name, "_steiger.csv")),
              row.names = FALSE)
  }

  # 保存MR-PRESSO结果
  if (!is.null(presso_res)) {
    write.csv(presso_res, file.path(tables_dir, paste0(gene_name, "_mr_presso.csv")),
              row.names = FALSE)
  }
  
  # 存储到汇总列表
  all_analyses[[gene_name]] <- list(
    harmonised_data = dat,
    mr_results = mr_res,
    heterogeneity = het_res,
    pleiotropy = pleio_res,
    steiger = steiger_res,
    mr_presso = presso_res
  )
  
  # 记录成功
  results_summary_pqtl <- rbind(results_summary_pqtl,
                                data.frame(exposure = gene_name, status = "Success"))
  
  cat("\n✓ Analysis for", gene_name, "completed successfully!\n")
}

# 4. 保存成功分析的基因列表
cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 4: SAVING SUCCESSFULLY ANALYZED PROTEINS LIST\n")
cat(rep("=", 60), "\n")

# 提取成功分析的基因
successful_genes <- results_summary_pqtl %>% 
  filter(status == "Success") %>% 
  dplyr::select(exposure) %>% 
  distinct()
if (nrow(successful_genes) > 0) {
  successful_genes$exposure <- ifelse(
    successful_genes$exposure %in% names(analysis_display_map),
    unname(analysis_display_map[successful_genes$exposure]),
    successful_genes$exposure
  )
}

cat("  Total proteins successfully analyzed:", nrow(successful_genes), "\n")

if (nrow(successful_genes) > 0) {
  # 保存成功分析的基因列表
  successful_file <- file.path(base_dir, "09.pQTL_Successful_Genes.csv")
  write.csv(successful_genes, successful_file, row.names = FALSE)
  cat("  ✓ Successfully analyzed proteins list saved to:", successful_file, "\n")
  
  # 显示成功基因
  cat("\n  Successfully analyzed proteins (", nrow(successful_genes), "):\n")
  for (i in 1:nrow(successful_genes)) {
    cat("  ", i, ". ", successful_genes$exposure[i], "\n", sep = "")
  }
} else {
  cat("  ✗ No proteins successfully analyzed.\n")
  
  # 保存空列表
  successful_genes <- data.frame(exposure = character(), stringsAsFactors = FALSE)
  successful_file <- file.path(base_dir, "09.pQTL_Successful_Genes.csv")
  write.csv(successful_genes, successful_file, row.names = FALSE)
  cat("  ✓ Empty successful genes list saved to:", successful_file, "\n")
}

# 5. 生成综合报告
cat("\n", rep("=", 60), "\n", sep="")
cat("COMPREHENSIVE pQTL MR ANALYSIS REPORT\n")
cat(rep("=", 60), "\n\n", sep="")

cat("STUDY OVERVIEW:\n")
cat("• Analysis Type: pQTL (Protein Quantitative Trait Loci)\n")
cat("• Candidate genes processed:", length(candidate_genes), "\n")
cat("• pQTL files found:", length(pqtl_files), "\n")
cat("• Initially matched proteins:", length(matched_genes), "\n")
cat("• Successfully analyzed:", sum(results_summary_pqtl$status == "Success"), "\n")
cat("• Failed:", sum(results_summary_pqtl$status != "Success"), "\n")
cat("• Outcomes:", paste(outcome_id, collapse = ", "), "\n\n")

if (nrow(results_summary_pqtl) > 0) {
  cat("ANALYSIS STATUS SUMMARY:\n")
  status_table <- table(results_summary_pqtl$status)
  print(status_table)
  cat("\n")
  
  # 汇总主要结果 - 只包括成功分析的部分
  if (exists("all_analyses") && length(all_analyses) > 0) {
    cat("MAIN MR RESULTS (IVW method) - SUCCESSFUL ANALYSES:\n")
    all_ivw_results <- data.frame()
    for (gene_name in names(all_analyses)) {
      analysis <- all_analyses[[gene_name]]
      mr_res <- analysis$mr_results
      
      if (!is.null(mr_res) && nrow(mr_res) > 0) {
        ivw_rows <- mr_res[grepl("Inverse variance weighted", mr_res$method, ignore.case = TRUE), ]
        
        if (nrow(ivw_rows) > 0) {
          ivw_row <- ivw_rows[1, ]
          
          all_ivw_results <- rbind(all_ivw_results, 
                                   data.frame(
                                     protein = gene_name,
                                     n_snps = ivw_row$nsnp,
                                     beta = round(ivw_row$b, 4),
                                     se = round(ivw_row$se, 4),
                                     p_value = ifelse(!is.na(ivw_row$pval), format.pval(ivw_row$pval, digits = 3), "NA"),
                                     OR = ifelse(!is.na(ivw_row$or), round(ivw_row$or, 3), NA),
                                     CI_95 = ifelse(!is.na(ivw_row$or_lci95) && !is.na(ivw_row$or_uci95), 
                                                    paste0("[", round(ivw_row$or_lci95, 3), ", ", round(ivw_row$or_uci95, 3), "]"), "NA"),
                                     method_used = ivw_row$method,
                                     stringsAsFactors = FALSE
                                   )
          )
        }
      }
    }
    
    if (nrow(all_ivw_results) > 0) {
      print(all_ivw_results)
      
      # 保存主要结果汇总
      write.csv(all_ivw_results, 
                file = file.path(tables_dir, "pqtl_main_results_summary.csv"),
                row.names = FALSE)
      
      # 统计显著结果
      p_values_numeric <- suppressWarnings(as.numeric(all_ivw_results$p_value))
      n_sig <- sum(p_values_numeric < 0.05, na.rm = TRUE)
      cat("\nSignificant results (p < 0.05):", n_sig, "/", nrow(all_ivw_results), "\n")
      
      if (n_sig > 0) {
        cat("Significant proteins:\n")
        sig_indices <- which(p_values_numeric < 0.05)
        sig_proteins <- all_ivw_results[sig_indices, ]
        print(sig_proteins[, c("protein", "OR", "CI_95", "p_value")])
        
        # 保存显著结果
        write.csv(sig_proteins, 
                  file = file.path(tables_dir, "pqtl_significant_results.csv"),
                  row.names = FALSE)
      }
    } else {
      cat("No IVW results available.\n")
    }
  }
} else {
  cat("No analyses completed successfully.\n")
}

cat("\nOUTPUT DIRECTORY STRUCTURE:\n")
cat("1. Tables directory:", tables_dir, "\n")
cat("2. Figures directory:", figures_dir, "\n")
cat("3. Initial matching:", file.path(base_dir, "08.pQTL_ID_Mapping.csv"), "\n")
cat("4. Successful analyses:", file.path(base_dir, "09.pQTL_Successful_Genes.csv"), "\n")
cat("5. Analysis log:", log_file, "\n")
cat("6. Main results summary:", file.path(tables_dir, "pqtl_main_results_summary.csv"), "\n")

# ============================================================================
# 步骤5: 生成与其他分析格式一致的输出文件
# ============================================================================
cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 5: GENERATING STANDARD OUTPUT FILES\n")
cat(rep("=", 60), "\n")

# 初始化整合结果表格
integrated_results <- data.frame()
het_pleio_results <- data.frame()

# 处理所有分析结果
if (exists("all_analyses") && length(all_analyses) > 0) {
  cat("Processing", length(all_analyses), "analyses to generate summary files...\n")

  for (exp_id in names(all_analyses)) {
    analysis <- all_analyses[[exp_id]]

    # 检查是否有MR结果
    if (!is.null(analysis$mr_results) && nrow(analysis$mr_results) > 0) {

      # 获取SNP数量
      n_snps <- ifelse(!is.null(analysis$harmonised_data), nrow(analysis$harmonised_data), NA)

      # 提取各方法结果
      ivw_rows <- analysis$mr_results[grepl("Inverse variance weighted", analysis$mr_results$method, ignore.case = TRUE), ]
      egger_rows <- analysis$mr_results[grepl("MR Egger", analysis$mr_results$method, ignore.case = TRUE), ]
      wm_rows <- analysis$mr_results[grepl("Weighted median", analysis$mr_results$method, ignore.case = TRUE), ]
      wmode_rows <- analysis$mr_results[grepl("Weighted mode", analysis$mr_results$method, ignore.case = TRUE), ]
      smode_rows <- analysis$mr_results[grepl("Simple mode", analysis$mr_results$method, ignore.case = TRUE), ]

      # 为每个结局创建结果行
      for (outcome_id_val in outcome_id) {

        # IVW结果
        if (nrow(ivw_rows) > 0) {
          ivw_row <- ivw_rows[1, ]

          # 提取OR和置信区间
          ivw_beta <- ivw_row$b
          ivw_se <- ivw_row$se
          ivw_pval <- ifelse(!is.na(ivw_row$pval), format.pval(ivw_row$pval, digits = 3), "NA")
          ivw_or <- ifelse(!is.na(ivw_row$or), ivw_row$or, NA)
          ivw_or_lci <- ifelse(!is.na(ivw_row$or_lci95), ivw_row$or_lci95, NA)
          ivw_or_uci <- ifelse(!is.na(ivw_row$or_uci95), ivw_row$or_uci95, NA)
          ivw_or_95ci <- ifelse(!is.na(ivw_or_lci) && !is.na(ivw_or_uci),
                                paste0(round(ivw_or_lci, 3), "-", round(ivw_or_uci, 3)), "NA")

          # 判断显著性
          ivw_pval_num <- safe_parse_pval(ivw_pval)
          ivw_significant <- ifelse(!is.na(ivw_pval_num) && ivw_pval_num < 0.05, "YES", "NO")

          # 判断方向
          ivw_direction <- ifelse(!is.na(ivw_or) && ivw_or > 1, "Risk",
                                  ifelse(!is.na(ivw_or) && ivw_or < 1, "Protective", "NA"))

          # 提取异质性结果
          ivw_het_q <- NA
          ivw_het_pval <- NA
          if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
            ivw_het_rows <- analysis$heterogeneity[grepl("Inverse variance weighted", analysis$heterogeneity$method, ignore.case = TRUE), ]
            if (nrow(ivw_het_rows) > 0) {
              ivw_het_q <- round(ivw_het_rows$Q[1], 3)
              ivw_het_pval <- ifelse(!is.na(ivw_het_rows$Q_pval[1]),
                                     format.pval(ivw_het_rows$Q_pval[1], digits = 3, eps = 1e-99), "NA")
            }
          }

          # 提取MR-Egger截距检验
          egger_intercept <- NA
          egger_intercept_pval <- NA
          if (!is.null(analysis$pleiotropy) && nrow(analysis$pleiotropy) > 0) {
            egger_intercept <- round(analysis$pleiotropy$egger_intercept[1], 6)
            egger_intercept_pval <- ifelse(!is.na(analysis$pleiotropy$pval[1]),
                                           format.pval(analysis$pleiotropy$pval[1], digits = 6, eps = 1e-99), "NA")
          }

          # 提取Steiger结果并预计算派生变量
          steiger_correct <- NA
          steiger_pval <- NA
          steiger_pval_clean <- NA
          steiger_pval_num <- NA_real_
          if (!is.null(analysis$steiger) && nrow(analysis$steiger) > 0) {
            steiger_correct <- ifelse(!is.na(analysis$steiger$correct_causal_direction[1]),
                                     ifelse(analysis$steiger$correct_causal_direction[1], "TRUE", "FALSE"), "NA")
            steiger_pval <- ifelse(!is.na(analysis$steiger$steiger_pval[1]),
                                  format.pval(analysis$steiger$steiger_pval[1], digits = 3, eps = 1e-99), "NA")
            steiger_pval_num <- safe_parse_pval(steiger_pval)
          }

          # 预计算水平多效性派生变量
          egger_intercept_pval_num <- safe_parse_pval(egger_intercept_pval)

          # 提取各方法结果
          egger_beta <- NA
          egger_se <- NA
          egger_pval <- NA
          egger_or <- NA
          if (nrow(egger_rows) > 0) {
            egger_beta <- egger_rows[1, ]$b
            egger_se <- egger_rows[1, ]$se
            egger_pval <- ifelse(!is.na(egger_rows[1, ]$pval), format.pval(egger_rows[1, ]$pval, digits = 3), "NA")
            egger_or <- egger_rows[1, ]$or
          }

          wm_beta <- NA
          wm_se <- NA
          wm_pval <- NA
          wm_or <- NA
          if (nrow(wm_rows) > 0) {
            wm_beta <- wm_rows[1, ]$b
            wm_se <- wm_rows[1, ]$se
            wm_pval <- ifelse(!is.na(wm_rows[1, ]$pval), format.pval(wm_rows[1, ]$pval, digits = 3), "NA")
            wm_or <- wm_rows[1, ]$or
          }

          wmode_beta <- NA
          wmode_se <- NA
          wmode_pval <- NA
          wmode_or <- NA
          if (nrow(wmode_rows) > 0) {
            wmode_beta <- wmode_rows[1, ]$b
            wmode_se <- wmode_rows[1, ]$se
            wmode_pval <- ifelse(!is.na(wmode_rows[1, ]$pval), format.pval(wmode_rows[1, ]$pval, digits = 3), "NA")
            wmode_or <- wmode_rows[1, ]$or
          }

          smode_beta <- NA
          smode_se <- NA
          smode_pval <- NA
          smode_or <- NA
          if (nrow(smode_rows) > 0) {
            smode_beta <- smode_rows[1, ]$b
            smode_se <- smode_rows[1, ]$se
            smode_pval <- ifelse(!is.na(smode_rows[1, ]$pval), format.pval(smode_rows[1, ]$pval, digits = 3), "NA")
            smode_or <- smode_rows[1, ]$or
          }

          # 创建整合结果行
          result_row <- data.frame(
            gene_symbol = ifelse(exp_id %in% names(analysis_display_map),
                                 unname(analysis_display_map[exp_id]),
                                 exp_id),
            ensembl_id = exp_id,
            exposure_id = paste0("pqtl-", exp_id),
            outcome_id = outcome_id_val,
            n_snps = n_snps,
            method = "IVW",
            ivw_beta = round(ivw_beta, 4),
            ivw_se = round(ivw_se, 4),
            ivw_pval = ivw_pval,
            ivw_or = round(ivw_or, 3),
            ivw_or_lci = round(ivw_or_lci, 3),
            ivw_or_uci = round(ivw_or_uci, 3),
            ivw_or_95ci = ivw_or_95ci,
            egger_beta = round(egger_beta, 4),
            egger_se = round(egger_se, 4),
            egger_pval = egger_pval,
            egger_or = round(egger_or, 3),
            wm_beta = round(wm_beta, 4),
            wm_se = round(wm_se, 4),
            wm_pval = wm_pval,
            wm_or = round(wm_or, 3),
            wmode_beta = round(wmode_beta, 4),
            wmode_se = round(wmode_se, 4),
            wmode_pval = wmode_pval,
            wmode_or = round(wmode_or, 3),
            smode_beta = round(smode_beta, 4),
            smode_se = round(smode_se, 4),
            smode_pval = smode_pval,
            smode_or = round(smode_or, 3),
            heterogeneity_q = ivw_het_q,
            heterogeneity_pval = ivw_het_pval,
            egger_intercept = egger_intercept,
            egger_intercept_pval = egger_intercept_pval,
            steiger_correct_direction = steiger_correct,
            steiger_pval = steiger_pval,
            ivw_significant = ivw_significant,
            ivw_direction = ivw_direction,
            # 三重筛选: 方向性 + 水平多效性 (派生变量已预计算)
            directionality_significant = ifelse(
              !is.na(steiger_correct) && (steiger_correct == "TRUE" || steiger_correct == TRUE) &&
              !is.na(steiger_pval_num) && steiger_pval_num < 0.05,
              "YES", "NO"
            ),
            # 水平多效性筛选: MR-Egger intercept p > 0.05 (无显著水平多效性)
            pleiotropy_pass = ifelse(
              !is.na(egger_intercept_pval_num) && egger_intercept_pval_num > 0.05,
              "YES", "NO"
            ),
            stringsAsFactors = FALSE
          )

          integrated_results <- rbind(integrated_results, result_row)
        }

        # 添加异质性和多效性结果
        if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
          het_rows <- analysis$heterogeneity
          for (h_idx in 1:nrow(het_rows)) {
            het_row <- het_rows[h_idx, ]
            het_pleio_row <- data.frame(
              exposure = ifelse(exp_id %in% names(analysis_display_map),
                                unname(analysis_display_map[exp_id]),
                                exp_id),
              outcome = outcome_id_val,
              Method = het_row$method,
              Heterogeneity_Q = round(het_row$Q, 4),
              Q_df = het_row$Q_df,
              Q_P.value = ifelse(!is.na(het_row$Q_pval), format.pval(het_row$Q_pval, digits = 4), "NA"),
              Egger_intercept = NA,
              SE = NA,
              P.value = NA,
              stringsAsFactors = FALSE
            )
            het_pleio_results <- rbind(het_pleio_results, het_pleio_row)
          }
        }

        if (!is.null(analysis$pleiotropy) && nrow(analysis$pleiotropy) > 0) {
          pleio_row <- analysis$pleiotropy[1, ]
          het_pleio_row <- data.frame(
            exposure = ifelse(exp_id %in% names(analysis_display_map),
                              unname(analysis_display_map[exp_id]),
                              exp_id),
            outcome = outcome_id_val,
            Method = "MR Egger",
            Heterogeneity_Q = NA,
            Q_df = NA,
            Q_P.value = NA,
            Egger_intercept = round(pleio_row$egger_intercept, 6),
            SE = round(pleio_row$se, 6),
            P.value = ifelse(!is.na(pleio_row$pval), format.pval(pleio_row$pval, digits = 4), "NA"),
            stringsAsFactors = FALSE
          )
          het_pleio_results <- rbind(het_pleio_results, het_pleio_row)
        }
      }
    }
  }
}

# ===== 从 all_analyses 提取 MR-PRESSO 结果并添加到 integrated_results =====
if (nrow(integrated_results) > 0 && exists("all_analyses")) {
  presso_data <- lapply(integrated_results$ensembl_id, function(g) {
    p <- all_analyses[[g]]$mr_presso
    if (!is.null(p) && is.data.frame(p) && nrow(p) > 0) {
      data.frame(
        presso_global_p = ifelse(!is.null(p$global_test_p[1]), safe_parse_pval(p$global_test_p[1]), NA_real_),
        presso_n_outliers = ifelse(!is.null(p$n_outliers[1]), as.integer(p$n_outliers[1]), NA_integer_),
        presso_outlier_snps = ifelse(!is.null(p$outlier_snps[1]), as.character(p$outlier_snps[1]), NA_character_),
        presso_outlier_corrected_beta = ifelse(!is.null(p$outlier_corrected_estimate[1]), p$outlier_corrected_estimate[1], NA_real_),
        presso_outlier_corrected_se = ifelse(!is.null(p$outlier_corrected_se[1]), p$outlier_corrected_se[1], NA_real_),
        presso_outlier_corrected_p = ifelse(!is.null(p$outlier_corrected_p[1]), safe_parse_pval(p$outlier_corrected_p[1]), NA_real_)
      )
    } else {
      data.frame(
        presso_global_p = NA_real_,
        presso_n_outliers = NA_integer_,
        presso_outlier_snps = NA_character_,
        presso_outlier_corrected_beta = NA_real_,
        presso_outlier_corrected_se = NA_real_,
        presso_outlier_corrected_p = NA_real_
      )
    }
  })
  presso_df <- do.call(rbind, presso_data)
  integrated_results <- cbind(integrated_results, presso_df)
}

# 保存整合结果
if (nrow(integrated_results) > 0) {
  cat("Saving integrated results...\n")

  # 提取正确的outcome名称（从all_analyses的heterogeneity数据中获取完整GWAS名称）
  outcome_name <- "disease"
  first_gene <- names(all_analyses)[1]
  if (!is.null(all_analyses[[first_gene]]$heterogeneity) &&
      nrow(all_analyses[[first_gene]]$heterogeneity) > 0) {
    raw_outcome <- as.character(all_analyses[[first_gene]]$heterogeneity$outcome[1])
    if (!is.na(raw_outcome) && grepl("\\|\\|", raw_outcome)) {
      outcome_name <- gsub(" \\|\\|.*$", "", raw_outcome)
    } else if (!is.na(raw_outcome) && raw_outcome != "") {
      outcome_name <- raw_outcome
    }
  }
  cat("  Outcome name:", outcome_name, "\n")

  # 按IVW P值排序
  integrated_results$ivw_pval_numeric <- vapply(integrated_results$ivw_pval, safe_parse_pval, numeric(1))
  integrated_results <- integrated_results[order(integrated_results$ivw_pval_numeric, na.last = NA), ]
  integrated_results$ivw_pval_numeric <- NULL

  # ===== 保存合并格式文件（00-07）=====
  cat("\nSaving merged format files (00-07)...\n")

  # 00. 最完整初始结果（1行/基因，所有column）
  integrated_results$outcome <- outcome_name
  file_00 <- file.path(base_dir, "00.Complete_MR_Results.csv")
  write.csv(integrated_results, file_00, row.names = FALSE)
  cat("  ✓ 00.Complete_MR_Results.csv saved\n")

  # 01. Table 1格式（5行/蛋白：5种MR方法，包含所有column）
  mr_table1_list <- list()
  for (i in 1:nrow(integrated_results)) {
    row <- integrated_results[i, ]
    gene_name <- as.character(row$gene_symbol)
    methods <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode", "Simple mode")
    for (m in methods) {
      mr_table1_list[[length(mr_table1_list) + 1]] <- data.frame(
        Exposure = gene_name,
        Outcome = outcome_name,
        Method = m,
        n_SNPs = row$n_snps,
        gene_symbol = gene_name,
        ensembl_id = as.character(row$ensembl_id),
        exposure_id = as.character(row$exposure_id),
        outcome_id = as.character(row$outcome_id),
        ivw_beta = row$ivw_beta,
        ivw_se = row$ivw_se,
        ivw_pval = row$ivw_pval,
        ivw_or = row$ivw_or,
        ivw_or_lci = row$ivw_or_lci,
        ivw_or_uci = row$ivw_or_uci,
        ivw_or_95ci = row$ivw_or_95ci,
        egger_beta = row$egger_beta,
        egger_se = row$egger_se,
        egger_pval = row$egger_pval,
        egger_or = row$egger_or,
        wm_beta = row$wm_beta,
        wm_se = row$wm_se,
        wm_pval = row$wm_pval,
        wm_or = row$wm_or,
        wmode_beta = row$wmode_beta,
        wmode_se = row$wmode_se,
        wmode_pval = row$wmode_pval,
        wmode_or = row$wmode_or,
        smode_beta = row$smode_beta,
        smode_se = row$smode_se,
        smode_pval = row$smode_pval,
        smode_or = row$smode_or,
        heterogeneity_q = row$heterogeneity_q,
        heterogeneity_pval = row$heterogeneity_pval,
        egger_intercept = row$egger_intercept,
        egger_intercept_pval = row$egger_intercept_pval,
        steiger_correct_direction = as.character(row$steiger_correct_direction),
        steiger_pval = as.character(row$steiger_pval),
        ivw_significant = as.character(row$ivw_significant),
        ivw_direction = as.character(row$ivw_direction),
        directionality_significant = as.character(row$directionality_significant),
        pleiotropy_pass = as.character(row$pleiotropy_pass),
        presso_global_p = row$presso_global_p,
        stringsAsFactors = FALSE
      )
    }
  }
  mr_table1_df <- do.call(rbind, mr_table1_list)
  file_01 <- file.path(base_dir, "01.MR_Results_All_Genes.csv")
  write.csv(mr_table1_df, file_01, row.names = FALSE)
  cat("  ✓ 01.MR_Results_All_Genes.csv saved\n")

  # 02. 异质性检验汇总（IVW + MR Egger）+ MR-PRESSO global test
  het_list <- list()
  safe_pick <- function(df, col_name) {
    if (is.null(df) || nrow(df) == 0 || !(col_name %in% colnames(df))) return(NA)
    val <- df[[col_name]][1]
    if (length(val) == 0) NA else val
  }
  for (gene_name in names(all_analyses)) {
    analysis <- all_analyses[[gene_name]]
    display_gene_name <- if (!is.null(analysis_display_map) &&
                             length(analysis_display_map) > 0 &&
                             gene_name %in% names(analysis_display_map)) {
      mapped <- unname(analysis_display_map[gene_name])
      if (length(mapped) == 0 || is.na(mapped) || !nzchar(mapped)) gene_name else mapped
    } else {
      gene_name
    }
    presso_p <- NA
    egger_int_pval <- NA
    if (!is.null(analysis$mr_presso) && !is.null(analysis$mr_presso$global_test_p)) {
      presso_p <- analysis$mr_presso$global_test_p
    }
    if (!is.null(analysis$pleiotropy) && nrow(analysis$pleiotropy) > 0) {
      egger_int_pval <- analysis$pleiotropy$egger_intercept_pval[1]
    }
    if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
      het_res <- analysis$heterogeneity
      for (meth in c("Inverse variance weighted", "MR Egger")) {
        mrow <- het_res[grep(meth, het_res$method, ignore.case = TRUE), ]
        if (nrow(mrow) > 0) {
          mrow <- mrow[1, ]
          het_list[[length(het_list) + 1]] <- data.frame(
            Exposure = display_gene_name,
            Outcome = outcome_name,
            Method = ifelse(grepl("Inverse", meth, ignore.case = TRUE),
                          "Inverse variance weighted", "MR Egger"),
            Heterogeneity_Q = safe_pick(mrow, "Q"),
            Q_df = safe_pick(mrow, "Q_df"),
            Q_P_value = safe_pick(mrow, "Q_pval"),
            Egger_intercept_pval = egger_int_pval,
            MR_PRESSO_global_p = presso_p,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(het_list) > 0) {
    het_df <- do.call(rbind, het_list)
  } else {
    het_df <- data.frame(
      Exposure = character(0),
      Outcome = character(0),
      Method = character(0),
      Heterogeneity_Q = numeric(0),
      Q_df = numeric(0),
      Q_P_value = numeric(0),
      Egger_intercept_pval = numeric(0),
      MR_PRESSO_global_p = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  het_file <- file.path(base_dir, "02.Heterogeneity_All_Genes.csv")
  write.csv(het_df, het_file, row.names = FALSE)
  cat("  ✓ 02.Heterogeneity_All_Genes.csv saved\n")

  # 03. 水平多效性检验汇总
  pleio_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = outcome_name,
    Egger_intercept = integrated_results$egger_intercept,
    Egger_intercept_pval = integrated_results$egger_intercept_pval,
    Pleiotropy_Significant = ifelse(
      safe_parse_pval(integrated_results$egger_intercept_pval) < 0.05,
      "YES", "NO"
    ),
    stringsAsFactors = FALSE
  )
  pleio_file <- file.path(base_dir, "03.Pleiotropy_All_Genes.csv")
  write.csv(pleio_all, pleio_file, row.names = FALSE)
  cat("  ✓ 03.Pleiotropy_All_Genes.csv saved\n")

  # 04. MR-PRESSO汇总
  presso_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = outcome_name,
    MR_PRESSO_global_p = integrated_results$presso_global_p,
    n_outliers = integrated_results$presso_n_outliers,
    outlier_snps = integrated_results$presso_outlier_snps,
    outlier_corrected_beta = integrated_results$presso_outlier_corrected_beta,
    outlier_corrected_se = integrated_results$presso_outlier_corrected_se,
    outlier_corrected_p = integrated_results$presso_outlier_corrected_p,
    stringsAsFactors = FALSE
  )
  presso_file <- file.path(base_dir, "04.MR_PRESSO_Results.csv")
  write.csv(presso_all, presso_file, row.names = FALSE)
  cat("  ✓ 04.MR_PRESSO_Results.csv saved\n")

  # 05. Steiger方向性分析汇总
  steiger_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = outcome_name,
    Steiger_correct_direction = integrated_results$steiger_correct_direction,
    Steiger_pval = integrated_results$steiger_pval,
    Directionality_significant = integrated_results$directionality_significant,
    stringsAsFactors = FALSE
  )
  steiger_file <- file.path(base_dir, "05.Steiger_Direction_Results.csv")
  write.csv(steiger_all, steiger_file, row.names = FALSE)
  cat("  ✓ 05.Steiger_Direction_Results.csv saved\n")

  # 06. 三重筛选结果汇总
  three_filter_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = outcome_name,
    IVW_significant = integrated_results$ivw_significant,
    Directionality_significant = integrated_results$directionality_significant,
    Pleiotropy_pass = integrated_results$pleiotropy_pass,
    Three_filter_pass = ifelse(
      integrated_results$ivw_significant == "YES" &
      integrated_results$directionality_significant == "YES" &
      integrated_results$pleiotropy_pass == "YES",
      "YES", "NO"
    ),
    stringsAsFactors = FALSE
  )
  tf_file <- file.path(base_dir, "06.Three_Filter_Results.csv")
  write.csv(three_filter_all, tf_file, row.names = FALSE)
  cat("  ✓ 06.Three_Filter_Results.csv saved\n")

  # 07. 因果靶点筛选结果（三重筛选通过）
  causal_targets <- integrated_results[
    integrated_results$ivw_significant == "YES" &
    integrated_results$directionality_significant == "YES" &
    integrated_results$pleiotropy_pass == "YES",
  ]
  if (nrow(causal_targets) > 0) {
    causal_file <- file.path(base_dir, "07.MR_res_gene.csv")
    write.csv(causal_targets[, c("gene_symbol", "ensembl_id", "n_snps", "ivw_beta", "ivw_se", "ivw_pval", "ivw_or")],
              causal_file, row.names = FALSE)
    cat("  ✓ 07.MR_res_gene.csv saved (", nrow(causal_targets), " causal targets)\n")
  }

} else {
  cat("No integrated results to save.\n")
  # 创建空的00文件
  empty_00 <- integrated_results[0, ]
  file_00 <- file.path(base_dir, "00.Complete_MR_Results.csv")
  write.csv(empty_00, file_00, row.names = FALSE)
  cat("  ✓ 00.Complete_MR_Results.csv saved (empty)\n")
}

cat("\n✓ All merged format files (00-07) generated successfully!\n")

# 如果有失败分析，显示失败文件位置
if (exists("failed_file") && file.exists(failed_file)) {
  cat("7. Failed analyses:", failed_file, "\n")
}

cat("\nAnalysis completed at:", format(Sys.time()), "\n")

# 关闭日志文件
sink()

cat("\n✓ pQTL MR analysis completed!\n")
cat("Total proteins successfully analyzed:", sum(results_summary_pqtl$status == "Success"), "\n")
cat("Total proteins failed:", sum(results_summary_pqtl$status != "Success"), "\n")

cat("\n主要结果文件 (00-07 合并格式):\n")
cat("   00. 00.Complete_MR_Results.csv       - 最完整初始结果（1行/基因）\n")
cat("   01. 01.MR_Results_All_Genes.csv     - Table 1格式（5行/基因：5种MR方法）\n")
cat("   02. 02.Heterogeneity_All_Genes.csv  - 异质性检验（IVW + MR Egger + PRESSO）\n")
cat("   03. 03.Pleiotropy_All_Genes.csv     - 水平多效性检验\n")
cat("   04. 04.MR_PRESSO_Results.csv        - MR-PRESSO结果\n")
cat("   05. 05.Steiger_Direction_Results.csv - Steiger方向性分析\n")
cat("   06. 06.Three_Filter_Results.csv     - 三重筛选结果\n")
cat("   07. 07.MR_res_gene.csv             - 因果靶点筛选结果\n")
cat("   08. 08.pQTL_ID_Mapping.csv                  - 蛋白ID映射表\n")
cat("   09. 09.pQTL_Successful_Genes.csv       - 成功分析的蛋白列表\n")
cat("   10. MR_pQTL_Analysis_Log.txt       - 分析日志\n\n")

cat("\n", rep("=", 70), "\n", sep="")
cat("pQTL MR ANALYSIS PIPELINE FINISHED\n")
cat(rep("=", 70), "\n")

# Save session info + qs2 snapshot
run_tag <- if (!is.null(opt$timestamp) && nzchar(opt$timestamp)) {
  opt$timestamp
} else {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}
dataset_name <- basename(normalizePath(base_dir, mustWork = FALSE))
qs2_path <- file.path(base_dir, sprintf("%s_pqtl_node.qs2", dataset_name))
save_qs2_snapshot(list(
  integrated_results = integrated_results,
  results_summary = results_summary_pqtl,
  successful_genes = successful_genes,
  matched_info = matched_info
), qs2_path)
save_session_info_snapshot(file.path(base_dir, "logs"), dataset_name, run_tag)

# Auto-generate Word report
cat("\n", rep("=", 60), "\n", sep = "")
cat("STEP 11: AUTO-GENERATING MR REPORT\n")
cat(rep("=", 60), "\n\n")
cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- grep("^--file=", cmd_args, value = TRUE)
if (length(script_path) > 0) {
  script_dir <- dirname(gsub("^--file=", "", script_path[1]))
} else {
  script_dir <- "/media/desk16/share/secure/MR/scripts"
}
report_script <- file.path(script_dir, "generate_mr_report.R")
if (file.exists(report_script)) {
  report_output <- file.path(base_dir, paste0("MR_Report_", opt$timestamp, ".docx"))
  report_cmd <- paste("Rscript", report_script,
                      "-i", base_dir,
                      "-o", report_output,
                      "-t", "pqtl",
                      "-T", opt$timestamp)
  cat("Generating report...\n")
  system(report_cmd)
  cat("✓ Report saved to:", report_output, "\n")
} else {
  cat("Warning: generate_mr_report.R not found, skipping auto-report generation\n")
}