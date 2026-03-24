# 0. 辅助函数：解析可能带 "<" 前缀的 p 值
safe_parse_pval <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  x <- as.character(x)
  x <- gsub("^\\s*<\\s*", "", x)  # 去掉 "<" 前缀和空白
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
      cat("   Ensembl->Symbol mapping failed:", e$message, "\n")
      setNames(rep(NA_character_, length(ids)), ids)
    }
  )

  mapped_chr <- as.character(mapped)
  names(mapped_chr) <- names(mapped)
  mapped_chr
}

# 0c. 统一展示名解析（与pQTL/traditional对齐）
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

# 与pQTL/traditional保持一致的安全取值函数
safe_pick <- function(df, col_name) {
  if (is.null(df) || nrow(df) == 0 || is.null(col_name) || !(col_name %in% colnames(df))) return(NA)
  val <- df[[col_name]][1]
  if (is.null(val) || length(val) == 0) return(NA)
  val
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
font_add("Liberation Sans", regular = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
         bold = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
         italic = "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf")
library(configr)
library(clusterProfiler)
# Config reading function
read_mr_config <- function() {
  # Try to find config file
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

  if (!is.null(config_file)) {
    cat("[Config] Loading config from:", config_file, "\n")
    config <- read.config(config_file)
    return(config)
  } else {
    cat("[Config] No config file found, using defaults\n")
    return(NULL)
  }
}

# Timestamp function
timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

# =============================================================================
# Token Management Functions
# =============================================================================

#' Test if token is valid by making a lightweight API call
#' @return list with 'valid' (TRUE/FALSE) and 'message'
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
#' @param expr Expression to evaluate
#' @param max_retries Maximum number of retries for non-token errors
#' @param on_token_expired Callback function when token expires
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

    # Check if it's a token expiry error
    if (!is.null(last_error) && is_token_expired_error(last_error$message)) {
      cat("[RETRY] Token expired! No more retries for expired token.\n")
      if (!is.null(on_token_expired)) {
        on_token_expired()
      }
      stop(last_error)
    }

    # Non-token error, retry with exponential backoff
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
  # Priority 1: Command line token
  if (!is.null(cmd_token) && nchar(cmd_token) > 0) {
    cat("[Token] Using token from command line\n")
    Sys.setenv(OPENGWAS_JWT = cmd_token)
    return(cmd_token)
  }

  # Priority 2: Environment variable (already set)
  token <- Sys.getenv("OPENGWAS_JWT")
  if (nchar(token) > 0) {
    cat("[Token] Using token from environment variable\n")
    return(token)
  }

  # Priority 3: .Renviron file
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
  make_option(c("-t", "--token"), type = "character", default = NULL,
              help = "OpenGWAS JWT token (optional, overrides env/file token)"),
  make_option(c("--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run [default auto-generated]")
)

# Read config file if exists
mr_config <- read_mr_config()

opt <- parse_args(OptionParser(option_list = option_list))

# Load token with priority: command line > env > .Renviron
token <- load_token(opt$token)
token_loaded <- !is.null(token) && nchar(token) > 0

# Print header
cat("\n")
cat("================================================================================\n")
cat("  MR eQTL Analysis\n")
cat("  Started:", timestamp(), "\n")
cat("================================================================================\n\n")

# Set Liberation Sans as default font
par(family = "Liberation Sans")

# Test token validity at startup
cat("\n[Token Check] Validating OpenGWAS token...\n")
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
  paste0("type = eqtl"),
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

cat("\n", rep("=", 60), "\n", sep="") 
cat("STEP 0: PROCESSING CANDIDATE GENE LIST\n")
cat(rep("=", 60), "\n")


# 1. 智能识别基因列
cat("\n1. Loading and inspecting the gene list...\n")
if (is.null(opt$input)) {
  stop("Error: Input file is required. Use -i or --input to specify the gene list file.")
}
genes <- read.csv(opt$input, check.names = FALSE)

cat("   Original dimensions:", nrow(genes), "rows,", ncol(genes), "columns\n")
cat("   Original column names:", paste(colnames(genes), collapse = ", "), "\n")
cat("   First few rows:\n")
print(head(genes))

# 智能识别基因列：按优先级匹配列名，同时用内容验证
gene_col <- NULL
id_col   <- NULL
gene_col_names <- c("gene", "symbol", "hgnc_symbol", "gene_symbol", "geneid",
                    "entrez", "entrezid", "entrez_gene_id",
                    "ensembl", "ensembl_id", "ensembl_gene_id")
id_col_names <- c("id", "snp_id", "snpid", "rsid", "variant_id", "chr", "chromosome")

for (nm in gene_col_names) {
  hit <- which(tolower(colnames(genes)) == nm)
  if (length(hit) > 0) { gene_col <- colnames(genes)[hit[1]]; break }
}
for (nm in id_col_names) {
  hit <- which(tolower(colnames(genes)) == nm)
  if (length(hit) > 0) { id_col <- colnames(genes)[hit[1]]; break }
}

# 备用策略：如果列名未匹配到，用内容验证（排除全数字列和过短列）
if (is.null(gene_col)) {
  for (i in seq_len(ncol(genes))) {
    col_vals <- as.character(genes[[i]])
    col_vals <- na.omit(col_vals)
    if (length(col_vals) == 0) next
    # 排除全数字列（序号列）或只有1-2字符的列
    is_numeric <- all(grepl("^\\d+(\\.\\d+)?$", col_vals))
    is_too_short <- all(nchar(col_vals) <= 2)
    is_gene_like <- any(grepl("^[A-Z][A-Z0-9]{2,}$", toupper(col_vals)))
    if (!is_numeric && !is_too_short && is_gene_like) {
      gene_col <- colnames(genes)[i]
      cat("   ✓ Gene column inferred by content validation: '", gene_col, "'\n", sep = "")
      break
    }
  }
}

if (is.null(gene_col)) {
  gene_col <- colnames(genes)[1]
  if (ncol(genes) > 1) id_col <- colnames(genes)[2]
  cat("   [WARN] No recognized gene column found. Using first column: '", gene_col, "'\n", sep = "")
} else {
  cat("   ✓ Gene column detected: '", gene_col, "'\n", sep = "")
}
if (!is.null(id_col)) {
  cat("   ✓ ID column detected: '", id_col, "'\n", sep = "")
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
  cat("   ✓ Added auto-generated id column\n")
}

# 验证基因列包含有效标识符
raw_gene_list <- unique(na.omit(genes$gene))
cat("   Total unique identifiers:", length(raw_gene_list), "\n")
sample_sample <- head(raw_gene_list, min(5, length(raw_gene_list)))
cat("   Sample values:", paste(sample_sample, collapse = ", "), "\n")

# 2. 标识符清洗与转换
cat("\n2. Cleaning and converting gene identifiers to Ensembl ID...\n")

# 判断输入ID类型
sample_genes <- head(raw_gene_list, 10)
is_ensembl <- any(grepl("^ENSG[0-9]+", sample_genes))
is_entrez <- any(grepl("^\\d+$", sample_genes)) & !is_ensembl

if (is_ensembl) {
  cat("   Input appears to be Ensembl IDs. Proceeding directly.\n")
  ensembl_id_list <- raw_gene_list
  conversion_method <- "None (already Ensembl)"
  conversion_successful <- TRUE

  # 创建映射表 - 检查是否已有formatted_exposure_id
  if ("formatted_exposure_id" %in% colnames(genes)) {
    # 使用输入文件中已有的formatted_exposure_id
    mapping_df <- data.frame(
      original_id = genes$gene,
      ensembl_id = genes$gene,
      formatted_exposure_id = genes$formatted_exposure_id
    )
  } else {
    # 自动生成formatted_exposure_id
    mapping_df <- data.frame(
      original_id = genes$gene,
      ensembl_id = genes$gene,
      formatted_exposure_id = paste0("eqtl-a-", genes$gene)
    )
  }
} else {
  cat("   Input does not appear to be Ensembl IDs. Attempting conversion...\n")
  conversion_successful <- FALSE
  ensembl_id_list <- c()
  mapping_df <- NULL
  
  # 方案A: 使用 clusterProfiler 和 org.Hs.eg.db
  tryCatch({
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is not installed.")
    }
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is not installed.")
    }
    
    library(clusterProfiler)
    library(org.Hs.eg.db)

    # 尝试从基因符号转换
    cat("   Attempting conversion from gene symbol to Ensembl ID...\n")
    cat("   Loading clusterProfiler and org.Hs.eg.db...\n")
    mapped <- tryCatch({
      bitr(raw_gene_list,
             fromType = "SYMBOL",
             toType = "ENSEMBL",
             OrgDb = org.Hs.eg.db,
             drop = FALSE)
    }, error = function(e) {
      cat("   bitr error:", e$message, "\n")
      NULL
    })
    
    if (nrow(mapped) > 0 && sum(!is.na(mapped$ENSEMBL)) > 0) {
      # 处理一对多映射：每个原始基因只保留第一个Ensembl ID
      mapped_unique <- mapped[!duplicated(mapped$SYMBOL) | !is.na(mapped$ENSEMBL), ]
      
      # 创建完整的映射表，确保行数与原始genes数据框一致
      mapping_df <- data.frame(
        original_id = genes$gene  # 改为使用"gene"列
      )
      
      # 合并转换结果
      mapping_df <- merge(mapping_df, 
                          mapped_unique[, c("SYMBOL", "ENSEMBL")], 
                          by.x = "original_id", 
                          by.y = "SYMBOL", 
                          all.x = TRUE)
      
      # 移除版本号
      mapping_df$ensembl_id_clean <- sub("\\..*", "", mapping_df$ENSEMBL)
      ensembl_id_list <- unique(na.omit(mapping_df$ensembl_id_clean))
      
      conversion_method <- "Gene Symbol to Ensembl (via org.Hs.eg.db)"
      conversion_successful <- TRUE
      cat("   ✓ Successfully converted", length(ensembl_id_list), "unique Ensembl IDs.\n")
      cat("   Note:", nrow(mapped) - nrow(mapped_unique), "duplicate mappings were removed.\n")
    } else {
      cat("   ✗ Conversion from SYMBOL failed.\n")
    }
  }, error = function(e) {
    cat("   Conversion attempt failed with error:", e$message, "\n")
  })
  
  # 方案B: 如果A失败，尝试从Entrez ID转换
  if (!conversion_successful && is_entrez) {
    tryCatch({
      cat("   Attempting conversion from Entrez ID to Ensembl ID...\n")
      mapped <- bitr(as.character(raw_gene_list),
                     fromType = "ENTREZID",
                     toType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     drop = FALSE)
      
      if (nrow(mapped) > 0 && sum(!is.na(mapped$ENSEMBL)) > 0) {
        # 处理一对多映射
        mapped_unique <- mapped[!duplicated(mapped$ENTREZID) | !is.na(mapped$ENSEMBL), ]
        
        mapping_df <- data.frame(
          original_id = genes$gene  # 改为使用"gene"列
        )
        
        mapping_df <- merge(mapping_df, 
                            mapped_unique[, c("ENTREZID", "ENSEMBL")], 
                            by.x = "original_id", 
                            by.y = "ENTREZID", 
                            all.x = TRUE)
        
        mapping_df$ensembl_id_clean <- sub("\\..*", "", mapping_df$ENSEMBL)
        ensembl_id_list <- unique(na.omit(mapping_df$ensembl_id_clean))
        
        conversion_method <- "Entrez ID to Ensembl (via org.Hs.eg.db)"
        conversion_successful <- TRUE
        cat("   ✓ Successfully converted", length(ensembl_id_list), "unique Ensembl IDs.\n")
      }
    }, error = function(e) {
      cat("   Entrez ID conversion also failed:", e$message, "\n")
    })
  }
  
  # 方案C: 如果所有自动转换都失败
  if (!conversion_successful) {
    warning(paste(
      "Could not automatically convert gene identifiers to Ensembl format.\n",
      "Please ensure your input is either:\n",
      "  1. Ensembl IDs (starting with 'ENSG')\n",
      "  2. Official gene symbols (e.g., 'TP53')\n",
      "  3. Entrez IDs (numeric)\n",
      "Falling back to using raw input as exposure IDs."
    ))
    
    # 创建基本映射表
    mapping_df <- data.frame(
      original_id = genes$gene,  # 改为使用"gene"列
      ensembl_id_clean = genes$gene  # 改为使用"gene"列
    )
    
    ensembl_id_list <- unique(na.omit(genes$gene))  # 改为使用"gene"列
    conversion_method <- "None (conversion failed, using raw input)"
  }
  
  # 确保mapping_df有所有必要的列
  if (!"ensembl_id_clean" %in% colnames(mapping_df)) {
    mapping_df$ensembl_id_clean <- mapping_df$original_id
  }

  # 生成用于展示的 gene_symbol（优先 Ensembl 反查 symbol）
  symbol_map <- map_ensembl_to_symbol(mapping_df$ensembl_id_clean)
  symbol_vals <- unname(symbol_map[as.character(mapping_df$ensembl_id_clean)])
  symbol_vals <- ifelse(!is.na(symbol_vals) & nzchar(symbol_vals), symbol_vals, as.character(mapping_df$original_id))
  mapping_df$gene_symbol <- symbol_vals
  
  # 创建最终映射表
  mapping_df$formatted_exposure_id <- ifelse(
    !is.na(mapping_df$ensembl_id_clean) & mapping_df$ensembl_id_clean != "",
    paste0("eqtl-a-", mapping_df$ensembl_id_clean),
    NA
  )
}
# 3. 生成最终的eQTL暴露ID列表
cat("\n3. Finalizing eQTL exposure ID list...\n")

# 从mapping_df中提取唯一的、非空的formatted_exposure_id（优先使用）
if (exists("mapping_df") && !is.null(mapping_df) && "formatted_exposure_id" %in% colnames(mapping_df)) {
  # 优先使用输入文件中已有的formatted_exposure_id
  valid_ids <- unique(na.omit(mapping_df$formatted_exposure_id[mapping_df$formatted_exposure_id != ""]))
  if (length(valid_ids) > 0) {
    exposure_id_list <- valid_ids
  } else {
    # 回退：使用ensembl_id_clean
    valid_ensembl_ids <- unique(na.omit(mapping_df$ensembl_id_clean[mapping_df$ensembl_id_clean != ""]))
    exposure_id_list <- paste0("eqtl-a-", valid_ensembl_ids)
  }
} else if (exists("mapping_df") && !is.null(mapping_df)) {
  # 没有formatted_exposure_id列，使用ensembl_id_clean
  valid_ensembl_ids <- unique(na.omit(mapping_df$ensembl_id_clean[mapping_df$ensembl_id_clean != ""]))
  exposure_id_list <- paste0("eqtl-a-", valid_ensembl_ids)
} else {
  # 回退方案
  exposure_id_list <- paste0("eqtl-a-", ensembl_id_list)
}

# 去重
exposure_id_list <- unique(exposure_id_list)
cat("   Conversion method:", conversion_method, "\n")
cat("   Final unique eQTL exposure IDs:", length(exposure_id_list), "\n")

if (length(exposure_id_list) > 0) {
  cat("   Sample of first 5 formatted exposure IDs:\n")
  print(head(exposure_id_list, 5))
} else {
  cat("   WARNING: No exposure IDs generated!\n")
}

# 4. 保存映射关系
cat("\n4. Saving conversion mapping for verification...\n")
if (exists("mapping_df") && !is.null(mapping_df)) {
  # 简化输出，只保留关键列
  output_mapping <- data.frame(
    original_id = mapping_df$original_id,
    gene_symbol = mapping_df$gene_symbol,
    ensembl_id = if ("ensembl_id_clean" %in% colnames(mapping_df)) {
      mapping_df$ensembl_id_clean
    } else {
      mapping_df$ensembl_id
    },
    formatted_exposure_id = mapping_df$formatted_exposure_id
  )
  
  # 移除完全重复的行
  output_mapping <- output_mapping[!duplicated(output_mapping), ]
  
  write.csv(output_mapping,
            file = file.path(base_dir, "08.eQTL_ID_Mapping.csv"),
            row.names = FALSE, quote = FALSE)
  cat("   ✓ ID conversion mapping saved to '08.eQTL_ID_Mapping.csv'\n")
  cat("   Mapping file has", nrow(output_mapping), "rows.\n")
} else {
  cat("   ✗ Could not create mapping file.\n")
}

cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 0 COMPLETE. Proceeding to main MR analysis loop with", 
    length(exposure_id_list), "exposures.\n")
cat(rep("=", 60), "\n\n")

# 检查：如果exposure_id_list为空，则停止
if (length(exposure_id_list) == 0) {
  stop("ERROR: No exposure IDs generated. Cannot proceed with MR analysis.")
}



outcome_id <- opt$outcome 

# 2.2 工具变量筛选参数 (use command line argument)
clump_params <- list(p1 = as.numeric(opt$pval), r2 = as.numeric(opt$r2), kb = as.numeric(opt$kb))
f_stat_threshold <- as.numeric(opt$fstat)
eaf_threshold <- as.numeric(opt$eaf)

# 2.3 定义输出目录结构

tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")

# 创建子目录
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "01_scatter"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "02_forest"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "03_funnel"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "04_leaveoneout"), showWarnings = FALSE, recursive = TRUE)

# 3. 创建主日志文件
log_file <- file.path(base_dir, "MR_Analysis_Log.txt")
sink(log_file, split = TRUE)
cat("=== Exposures Mendelian Randomization Analysis Log ===\n")
cat("Start Time:", format(Sys.time()), "\n")
cat("Exposures: EQTL\n")
cat("Total exposures to process:", length(exposure_id_list), "\n")
cat("Outcome:", outcome_id, "\n")
cat("===========================================\n\n")

# 4. 初始化结果列表
results_summary_eqtl <- data.frame()
all_analyses <- list()

# 创建 exp_id -> gene_symbol 的映射（用于图和表文件名）
expid_to_genesymbol <- NULL
if (exists("mapping_df") && !is.null(mapping_df) &&
    "formatted_exposure_id" %in% colnames(mapping_df) &&
    "gene_symbol" %in% colnames(mapping_df)) {
  # 优先从 mapping_df 构建映射
  expid_to_genesymbol <- setNames(mapping_df$gene_symbol, mapping_df$formatted_exposure_id)
} else if (exists("output_mapping") && !is.null(output_mapping)) {
  display_col <- if ("gene_symbol" %in% colnames(output_mapping)) "gene_symbol" else "original_id"
  expid_to_genesymbol <- setNames(output_mapping[[display_col]], output_mapping$formatted_exposure_id)
}

# 5. 主分析循环
for (i in seq_along(exposure_id_list)) {
  exp_id <- exposure_id_list[i]
  exp_name_simple <- gsub("[-._]", "_", exp_id)
  # 基因符号用于图文件名
  gene_symbol_file <- if (!is.null(expid_to_genesymbol) && exp_id %in% names(expid_to_genesymbol)) {
    expid_to_genesymbol[exp_id]
  } else {
    exp_name_simple
  }
  
  cat("\n", rep("=", 60), "\n", sep="")
  cat("ANALYSIS [", i, "/", length(exposure_id_list), "]: ", exp_id, "\n", sep="")
  cat(rep("=", 60), "\n\n", sep="")
  
  # 5.1 提取工具变量
  cat("STEP 1: Extracting instrumental variables...\n")
  exposure_dat <- tryCatch({
    # Use retry wrapper for API call
    exp_dat <- retry_with_token_check(
      extract_instruments(
        outcomes = exp_id,
        p1 = clump_params$p1,
        clump = TRUE,
        r2 = clump_params$r2,
        kb = clump_params$kb,
        opengwas_jwt = ieugwasr::get_opengwas_jwt(),
        force_server = FALSE
      ),
      max_retries = 3,
      on_token_expired = function() {
        cat("  FATAL: Token expired during extraction!\n")
      }
    )

    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      stop("No significant SNPs found (p < 5e-8)")
    }
    cat("  Extracted", nrow(exp_dat), "IVs\n")

    # 质量控制
    exp_dat <- exp_dat %>%
      filter(eaf.exposure > eaf_threshold & eaf.exposure < (1 - eaf_threshold))
    cat("  After EAF filter (0.01-0.99):", nrow(exp_dat), "IVs\n")

    # 尝试从API获取真实样本量，如果不可用则使用NA
    if (!"samplesize.exposure" %in% colnames(exp_dat)) {
      # 尝试从OpenGWAS API获取样本量信息 (with retry)
      gwasi <- retry_with_token_check(
        ieugwasr::gwasinfo(exp_id),
        max_retries = 3
      )
      if (!is.null(gwasi) && "nsample" %in% names(gwasi)) {
        exp_dat$samplesize.exposure <- gwasi$nsample
        cat("  Retrieved sample size from API:", exp_dat$samplesize.exposure[1], "\n")
      } else {
        exp_dat$samplesize.exposure <- NA
        cat("  WARNING: Sample size not available from API\n")
      }
    }
    exp_dat$f_statistic <- (exp_dat$beta.exposure^2) / (exp_dat$se.exposure^2)
    exp_dat <- exp_dat %>% filter(f_statistic > f_stat_threshold)
    cat("  After F-statistic filter (>", f_stat_threshold, "):", nrow(exp_dat), "IVs\n")
    
    if (nrow(exp_dat) < 3) {
      cat("  WARNING: Fewer than 3 IVs remaining\n")
      cat("  SKIPPING: Insufficient SNPs after F-statistic filter\n")
      return(NULL)
    }

    exp_dat
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(exposure_dat)) {
    results_summary_eqtl <- rbind(results_summary_eqtl, 
                                  data.frame(exposure = exp_id, status = "Failed_IV_Extraction"))
    next
  }
  
  # 5.2 提取结局数据
  cat("\nSTEP 2: Extracting outcome data...\n")
  outcome_dat <- tryCatch({
    out_dat <- retry_with_token_check(
      extract_outcome_data(
        snps = exposure_dat$SNP,
        outcomes = outcome_id
      ),
      max_retries = 3,
      on_token_expired = function() {
        cat("  FATAL: Token expired during outcome extraction!\n")
      }
    )
    if (is.null(out_dat) || nrow(out_dat) == 0) {
      stop("No matching SNPs in outcome")
    }
    cat("  Found", nrow(out_dat), "SNPs in outcome\n")
    out_dat
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(outcome_dat)) {
    results_summary_eqtl <- rbind(results_summary_eqtl,
                                  data.frame(exposure = exp_id, status = "Failed_Outcome_Extraction"))
    next
  }
  
  # 5.3 协调数据
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
    results_summary_eqtl <- rbind(results_summary_eqtl,
                                  data.frame(exposure = exp_id, status = "Failed_Harmonisation"))
    next
  }

  # 5.3.1 剔除与结局显著相关的SNP（排除水平多效性，p < 1e-5）
  cat("  Checking exclusivity assumption (removing outcome-associated SNPs, p < 1e-5)...\n")
  if ("pval.outcome" %in% colnames(dat)) {
    n_before <- nrow(dat)
    dat <- dat[dat$pval.outcome >= 1e-5, ]
    n_after <- nrow(dat)
    if (n_before > n_after) {
      cat("    Removed", n_before - n_after, "SNPs associated with outcome (p < 1e-5)\n")
    } else {
      cat("    ✓ All SNPs pass exclusivity check (no outcome-associated SNPs)\n")
    }
  }

  # 检查过滤后是否仍有足够SNP
  if (nrow(dat) < 3) {
    cat("  WARNING: Fewer than 3 SNPs after outcome-associated filter\n")
    cat("  SKIPPING: Insufficient SNPs for MR analysis\n")
    return(NULL)
  }

  # 5.4 异质性检验 (先运行以决定使用哪种IVW模型)
  cat("\nSTEP 4: Testing heterogeneity to select MR method...\n")
  het_res <- tryCatch({
    het <- mr_heterogeneity(dat)
    if (!is.null(het)) {
      het$exposure <- exp_id
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
    # 查找包含"Inverse variance"的方法 (更精确的匹配)
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

  # 5.5 运行主MR分析
  cat("\nSTEP 5: Running MR analysis...\n")
  mr_res <- tryCatch({
    res <- mr(dat, method_list = method_list)
    
    if (nrow(res) == 0) {
      stop("MR analysis returned no results")
    }
    
    cat("  Methods used:", paste(unique(res$method), collapse = ", "), "\n")
    
    # 计算OR
    res_or <- generate_odds_ratios(res)
    res_or$exposure <- exp_id
    res_or$outcome <- outcome_id
    res_or$n_snps <- nrow(dat)
    
    res_or
  }, error = function(e) {
    cat("  ERROR in MR analysis:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(mr_res)) {
    results_summary_eqtl <- rbind(results_summary_eqtl,
                                  data.frame(exposure = exp_id, status = "Failed_MR_Analysis"))
    next
  }

  # 5.6 多效性检验
  cat("\nSTEP 6: Pleiotropy test...\n")
  pleio_res <- tryCatch({
    pleio <- mr_pleiotropy_test(dat)
    if (!is.null(pleio)) {
      pleio$exposure <- exp_id
    }
    pleio
  }, error = function(e) {
    cat("  Pleiotropy test failed:", e$message, "\n")
    return(NULL)
  })
  
  # 5.7 Steiger方向性分析 - 修复版
  cat("\nSTEP 7: Steiger directionality test (fixed)...\n")
  steiger_res <- tryCatch({
    # 检查并添加必要的列 - 使用真实样本量或NA
    if (!"samplesize.exposure" %in% colnames(dat)) {
      dat$samplesize.exposure <- NA
    }
    if (!"samplesize.outcome" %in% colnames(dat)) {
      dat$samplesize.outcome <- NA
    }
    
    # 获取真实样本量用于r值计算
    n_exp <- dat$samplesize.exposure[1]
    n_out <- dat$samplesize.outcome[1]
    
    # 检查是否有r.exposure和r.outcome列
    has_r_exposure <- "r.exposure" %in% colnames(dat)
    has_r_outcome <- "r.outcome" %in% colnames(dat)
    
    if (!has_r_exposure || !has_r_outcome) {
      cat("  Adding r values for Steiger test...\n")
      
      # 使用实际样本量，如果不可用则设为NA
      ncase_exp <- ifelse(!is.na(n_exp), floor(n_exp * 0.5), NA)
      ncontrol_exp <- ifelse(!is.na(n_exp), floor(n_exp * 0.5), NA)
      ncase_out <- ifelse(!is.na(n_out), floor(n_out * 0.5), NA)
      ncontrol_out <- ifelse(!is.na(n_out), floor(n_out * 0.5), NA)
      
      # 使用NA让get_r_from_lor使用默认行为
      prevalence <- 0.1
      
      # 为每个SNP计算r值
      for (j in 1:nrow(dat)) {
        # 计算r.exposure
        if (!has_r_exposure) {
          beta <- dat$beta.exposure[j]
          eaf <- ifelse(is.na(dat$eaf.exposure[j]), 0.5, dat$eaf.exposure[j])
          if (!is.na(beta) && !is.na(eaf)) {
            dat$r.exposure[j] <- get_r_from_lor(
              lor = beta,
              af = eaf,
              ncase = ncase_exp,
              ncontrol = ncontrol_exp,
              prevalence = prevalence
            )
          } else {
            dat$r.exposure[j] <- 0.1
          }
        }
        
        # 计算r.outcome
        if (!has_r_outcome) {
          beta <- dat$beta.outcome[j]
          eaf <- ifelse(is.na(dat$eaf.outcome[j]), 0.5, dat$eaf.outcome[j])
          if (!is.na(beta) && !is.na(eaf)) {
            dat$r.outcome[j] <- get_r_from_lor(
              lor = beta,
              af = eaf,
              ncase = ncase_out,
              ncontrol = ncontrol_out,
              prevalence = prevalence
            )
          } else {
            dat$r.outcome[j] <- 0.1
          }
        }
      }
    }
    
    # 运行Steiger测试
    steiger <- directionality_test(dat)
    if (!is.null(steiger) && nrow(steiger) > 0) {
      steiger$exposure <- exp_id
      cat("  Steiger test completed successfully\n")
    } else {
      cat("  WARNING: Steiger test returned empty result\n")
      steiger <- data.frame(
        exposure = exp_id,
        snp_r2.exposure = NA,
        snp_r2.outcome = NA,
        correct_causal_direction = FALSE,
        steiger_pval = NA
      )
    }
    steiger
  }, error = function(e) {
    cat("  Steiger test failed:", e$message, "\n")
    # 返回一个空的Steiger结果
    data.frame(
      exposure = exp_id,
      snp_r2.exposure = NA,
      snp_r2.outcome = NA,
      correct_causal_direction = FALSE,
      steiger_pval = NA
    )
  })
  
  # 5.8 MR-PRESSO分析
  cat("\nSTEP 8: MR-PRESSO analysis...\n")
  presso_res <- tryCatch({
    if (nrow(dat) >= 4) {  # MR-PRESSO需要至少4个SNP
      # 直接调用MRPRESSO::mr_presso，传入字符串列名
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
        exposure = exp_id,
        outcome = outcome_id,
        global_test_p = NA,
        n_outliers = 0,
        outlier_snps = NA,
        outlier_corrected_estimate = NA,
        outlier_corrected_se = NA,
        outlier_corrected_p = NA
      )

      if (!is.null(presso) && !is.null(presso$`MR-PRESSO results`)) {
        # 全局检验
        gt <- presso$`MR-PRESSO results`$`Global Test`
        if (!is.null(gt) && length(gt) > 0) {
          cat("  Global Test result:\n")
          print(gt)
          if (!is.null(gt$Pvalue)) {
            presso_summary$global_test_p <- safe_parse_pval(gt$Pvalue)
          }
        }

        # 异常值
        if (!is.null(presso$`MR-PRESSO results`$`Outlier Test`)) {
          outlier_test <- presso$`MR-PRESSO results`$`Outlier Test`
          if (!is.null(outlier_test) && nrow(outlier_test) > 0) {
            # outlier SNP 名称
            outlier_names <- outlier_test$Outliers
            if (!is.null(outlier_names) && length(outlier_names) > 0) {
              presso_summary$n_outliers <- nrow(outlier_test)
              presso_summary$outlier_snps <- paste(outlier_names, collapse = ", ")
            }
          }
        }

        # 校正后估计
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
        exposure = exp_id,
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

  # 5.9 可视化
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

  # 散点图
  tryCatch({
    p <- mr_scatter_plot(mr_res, dat)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12, legend.position = "top") +
        ggtitle("MR Scatter Plot")

      num_prefix <- sprintf("%02d", i)
      # 保存PNG
      showtext_begin()
      ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", gene_symbol_file, ".scatter.png")),
             plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
      showtext_end()
      # 保存PDF
      showtext_begin()
      ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", gene_symbol_file, ".scatter.pdf")),
             plot = p1, width = 10, height = 9, device = "pdf")
      showtext_end()
      cat("  ✓ Scatter plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Scatter plot failed:", e$message, "\n")
  })

  # 森林图
  tryCatch({
    res_single <- mr_singlesnp(dat)
    p <- mr_forest_plot(res_single)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12) +
        ggtitle("MR Forest Plot")

      num_prefix <- sprintf("%02d", i)
      # 保存PNG - 增加宽度确保X轴标签不超出
      showtext_begin()
      ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_symbol_file, ".forest.png")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.5), dpi = 300, bg = "white")
      showtext_end()
      # 保存PDF
      showtext_begin()
      ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_symbol_file, ".forest.pdf")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.5), device = "pdf")
      showtext_end()
      # 保存rsID.csv
      if (!is.null(res_single) && nrow(res_single) > 0 && "SNP" %in% colnames(res_single)) {
        rsid_df <- data.frame(SNP = res_single$SNP, stringsAsFactors = FALSE)
        write.csv(rsid_df, file.path(figures_dir, "02_forest", paste0(num_prefix, ".", gene_symbol_file, ".rsID.csv")),
                  row.names = FALSE)
      }
      cat("  ✓ Forest plot saved (PNG + PDF + rsID.csv)\n")
    }
  }, error = function(e) {
    cat("  Forest plot failed:", e$message, "\n")
  })

  # 漏斗图
  tryCatch({
    res_single <- mr_singlesnp(dat)
    p <- mr_funnel_plot(res_single)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12, legend.position = "top") +
        ggtitle("MR Funnel Plot")

      num_prefix <- sprintf("%02d", i)
      # 保存PNG
      showtext_begin()
      ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", gene_symbol_file, ".funnel.png")),
             plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
      showtext_end()
      # 保存PDF
      showtext_begin()
      ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", gene_symbol_file, ".funnel.pdf")),
             plot = p1, width = 10, height = 9, device = "pdf")
      showtext_end()
      cat("  ✓ Funnel plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Funnel plot failed:", e$message, "\n")
  })

  # 留一法图
  tryCatch({
    res_loo <- mr_leaveoneout(dat)
    p <- mr_leaveoneout_plot(res_loo)
    if (!is.null(p)) {
      p1 <- p[[1]] +
        scale_color_npg(palette = "nrc") +
        theme_mr(base_size = 12) +
        ggtitle("MR Leave-One-Out Plot")

      num_prefix <- sprintf("%02d", i)
      # 保存PNG - 增加宽度确保X轴标签不超出
      showtext_begin()
      ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", gene_symbol_file, ".leaveoneout.png")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.4), dpi = 300, bg = "white")
      showtext_end()
      # 保存PDF
      showtext_begin()
      ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", gene_symbol_file, ".leaveoneout.pdf")),
             plot = p1, width = 16, height = max(10, nrow(dat) * 0.4), device = "pdf")
      showtext_end()
      cat("  ✓ Leave-one-out plot saved (PNG + PDF)\n")
    }
  }, error = function(e) {
    cat("  Leave-one-out plot failed:", e$message, "\n")
  })
  
  # 5.10 保存结果
  cat("\nSTEP 10: Saving results...\n")
  
  # 保存协调数据
  write.csv(dat, file.path(tables_dir, paste0(gene_symbol_file, "_harmonised.csv")),
            row.names = FALSE)

  # 保存MR结果
  write.csv(mr_res, file.path(tables_dir, paste0(gene_symbol_file, "_mr_results.csv")),
            row.names = FALSE)

  # 保存异质性结果
  if (!is.null(het_res)) {
    write.csv(het_res, file.path(tables_dir, paste0(gene_symbol_file, "_heterogeneity.csv")),
              row.names = FALSE)
  }

  # 保存多效性结果
  if (!is.null(pleio_res)) {
    write.csv(pleio_res, file.path(tables_dir, paste0(gene_symbol_file, "_pleiotropy.csv")),
              row.names = FALSE)
  }

  # 保存Steiger结果
  if (!is.null(steiger_res)) {
    write.csv(steiger_res, file.path(tables_dir, paste0(gene_symbol_file, "_steiger.csv")),
              row.names = FALSE)
  }

  # 保存MR-PRESSO结果
  if (!is.null(presso_res)) {
    write.csv(presso_res, file.path(tables_dir, paste0(gene_symbol_file, "_mr_presso.csv")),
              row.names = FALSE)
  }
  
  # 存储到汇总列表
  all_analyses[[exp_id]] <- list(
    harmonised_data = dat,
    mr_results = mr_res,
    heterogeneity = het_res,
    pleiotropy = pleio_res,
    steiger = steiger_res,
    mr_presso = presso_res
  )
  
  # 记录成功
  results_summary_eqtl <- rbind(results_summary_eqtl,
                                data.frame(exposure = exp_id, status = "Success"))
  
  cat("\n✓ Analysis for", exp_id, "completed successfully!\n")
}

# 6. 生成综合报告 - 完全修复版
cat("\n", rep("=", 60), "\n", sep="")
cat("COMPREHENSIVE ANALYSIS REPORT - eQTL MR ANALYSIS\n")
cat(rep("=", 60), "\n\n", sep="")

cat("STUDY OVERVIEW:\n")
cat("• Exposures: eQTL (expression quantitative trait loci) genetic variants\n")
cat("• Total exposures processed:", length(exposure_id_list), "\n")
cat("• Successfully analyzed:", sum(results_summary_eqtl$status == "Success"), "\n")
cat("• Failed:", sum(results_summary_eqtl$status != "Success"), "\n")
cat("• Outcome:", outcome_id, "\n\n")

# 6.1 汇总主要结果
cat("MAIN MR RESULTS (IVW method):\n")
all_ivw_results <- data.frame()
for (exp_id in names(all_analyses)) {
  analysis <- all_analyses[[exp_id]]
  mr_res <- analysis$mr_results
  
  if (!is.null(mr_res) && nrow(mr_res) > 0) {
    # 查找IVW方法
    ivw_rows <- mr_res[grepl("Inverse variance weighted", mr_res$method, ignore.case = TRUE), ]
    
    if (nrow(ivw_rows) > 0) {
      ivw_row <- ivw_rows[1, ]
      
      all_ivw_results <- rbind(all_ivw_results, 
                               data.frame(
                                 exposure = exp_id,
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
  
  # 统计显著结果
  p_values_numeric <- suppressWarnings(as.numeric(all_ivw_results$p_value))
  n_sig <- sum(p_values_numeric < 0.05, na.rm = TRUE)
  cat("\nSignificant results (p < 0.05):", n_sig, "/", nrow(all_ivw_results), "\n")
  
  if (n_sig > 0) {
    cat("Significant exposures:\n")
    sig_indices <- which(p_values_numeric < 0.05)
    sig_exposures <- all_ivw_results[sig_indices, ]
    print(sig_exposures[, c("exposure", "OR", "CI_95", "p_value")])
  }
} else {
  cat("No IVW results available.\n")
}

# 6.2 异质性总结 - 修复版
cat("\nHETEROGENEITY SUMMARY (Cochran's Q test):\n")
heterogeneity_summary <- data.frame()

for (exp_id in names(all_analyses)) {
  analysis <- all_analyses[[exp_id]]
  het <- analysis$heterogeneity
  
  if (!is.null(het) && is.data.frame(het) && nrow(het) > 0) {
    # 查找IVW的异质性结果
    ivw_het_rows <- het[grepl("Inverse variance weighted", het$method, ignore.case = TRUE), ]
    
    if (nrow(ivw_het_rows) > 0) {
      ivw_het <- ivw_het_rows[1, ]
    } else {
      ivw_het <- het[1, ]
    }
    
    # 安全提取值
    extract_safe <- function(x) {
      if (length(x) > 0 && !is.na(x)) {
        return(x)
      } else {
        return(NA)
      }
    }
    
    Q_val <- extract_safe(ivw_het$Q)
    Q_pval_val <- extract_safe(ivw_het$Q_pval)
    I2_val <- extract_safe(ivw_het$I2)
    method_val <- ifelse("method" %in% colnames(ivw_het), as.character(ivw_het$method), "NA")
    
    # 构建行 - 确保每个值都是标量
    result_row <- data.frame(
      exposure = exp_id,
      Q = ifelse(is.numeric(Q_val), round(Q_val, 2), NA),
      Q_pval = ifelse(is.numeric(Q_pval_val), format.pval(Q_pval_val, digits = 3), "NA"),
      I2 = ifelse(is.numeric(I2_val), paste0(round(I2_val, 1), "%"), "NA"),
      method_used = method_val,
      stringsAsFactors = FALSE
    )
    
    heterogeneity_summary <- rbind(heterogeneity_summary, result_row)
  }
}

if (nrow(heterogeneity_summary) > 0) {
  print(heterogeneity_summary)
  
  # 判断是否使用随机效应
  cat("\nINTERPRETATION OF HETEROGENEITY:\n")
  for (i in 1:nrow(heterogeneity_summary)) {
    row <- heterogeneity_summary[i, ]
    
    if (row$Q_pval != "NA") {
      pval <- suppressWarnings(as.numeric(row$Q_pval))
      if (!is.na(pval)) {
        if (pval < 0.05) {
          cat("  ", row$exposure, ": Significant heterogeneity (p = ", row$Q_pval, 
              ") → Use RANDOM EFFECTS IVW\n", sep = "")
        } else {
          cat("  ", row$exposure, ": No significant heterogeneity (p = ", row$Q_pval, 
              ") → Use FIXED EFFECTS IVW\n", sep = "")
        }
      } else {
        cat("  ", row$exposure, ": Cannot determine heterogeneity (p-value: ", row$Q_pval, ")\n", sep = "")
      }
    } else {
      cat("  ", row$exposure, ": Heterogeneity p-value not available\n", sep = "")
    }
  }
} else {
  cat("No heterogeneity results available.\n")
}

# 6.3 多效性总结
cat("\nPLEIOTROPY SUMMARY (MR-Egger intercept test):\n")
pleiotropy_summary <- data.frame()
for (exp_id in names(all_analyses)) {
  analysis <- all_analyses[[exp_id]]
  pleio <- analysis$pleiotropy
  
  if (!is.null(pleio) && nrow(pleio) > 0) {
    pleiotropy_summary <- rbind(pleiotropy_summary,
                                data.frame(
                                  exposure = exp_id,
                                  egger_intercept = round(pleio$egger_intercept, 4),
                                  se = round(pleio$se, 4),
                                  p_value = ifelse(!is.na(pleio$pval), format.pval(pleio$pval, digits = 3), "NA"),
                                  stringsAsFactors = FALSE
                                )
    )
  }
}

if (nrow(pleiotropy_summary) > 0) {
  print(pleiotropy_summary)
  
  # 判断多效性
  cat("\nINTERPRETATION OF PLEIOTROPY:\n")
  for (i in 1:nrow(pleiotropy_summary)) {
    row <- pleiotropy_summary[i, ]
    pval <- suppressWarnings(as.numeric(row$p_value))
    
    if (!is.na(pval)) {
      if (pval < 0.05) {
        cat("  ", row$exposure, ": Significant horizontal pleiotropy (p = ", row$p_value, 
            ") → Results may be biased\n", sep = "")
      } else {
        cat("  ", row$exposure, ": No significant horizontal pleiotropy (p = ", row$p_value, 
            ") → Results are reliable\n", sep = "")
      }
    } else {
      cat("  ", row$exposure, ": Cannot determine pleiotropy (p-value: ", row$p_value, ")\n", sep = "")
    }
  }
} else {
  cat("No pleiotropy results available.\n")
}

# 6.4 Steiger方向性总结
cat("\nSTEIGER DIRECTIONALITY TEST SUMMARY:\n")
steiger_summary <- data.frame()
for (exp_id in names(all_analyses)) {
  analysis <- all_analyses[[exp_id]]
  steiger <- analysis$steiger
  
  if (!is.null(steiger) && nrow(steiger) > 0) {
    # 检查列是否存在
    has_snp_r2_exposure <- "snp_r2.exposure" %in% colnames(steiger)
    has_snp_r2_outcome <- "snp_r2.outcome" %in% colnames(steiger)
    has_correct_causal_direction <- "correct_causal_direction" %in% colnames(steiger)
    has_steiger_pval <- "steiger_pval" %in% colnames(steiger)
    
    steiger_summary <- rbind(steiger_summary,
                             data.frame(
                               exposure = exp_id,
                               snp_r2.exposure = ifelse(has_snp_r2_exposure, round(steiger$snp_r2.exposure, 4), NA),
                               snp_r2.outcome = ifelse(has_snp_r2_outcome, round(steiger$snp_r2.outcome, 4), NA),
                               correct_causal_direction = ifelse(has_correct_causal_direction, 
                                                                 ifelse(steiger$correct_causal_direction, "TRUE", "FALSE"), 
                                                                 "FALSE"),
                               steiger_pval = ifelse(has_steiger_pval, 
                                                     ifelse(!is.na(steiger$steiger_pval), 
                                                            format.pval(steiger$steiger_pval, digits = 3), "NA"), 
                                                     "NA"),
                               stringsAsFactors = FALSE
                             )
    )
  }
}

if (nrow(steiger_summary) > 0) {
  print(steiger_summary)
  
  # 判断方向正确性
  cat("\nINTERPRETATION OF DIRECTIONALITY:\n")
  for (i in 1:nrow(steiger_summary)) {
    row <- steiger_summary[i, ]
    
    if (row$correct_causal_direction == "TRUE") {
      cat("  ", row$exposure, ": Correct causal direction\n", sep = "")
    } else {
      cat("  ", row$exposure, ": WARNING - Possible reverse causality\n", sep = "")
    }
  }
} else {
  cat("No Steiger test results available.\n")
}


cat("\nOUTPUT FILES:\n")
cat("1. Tables directory:", tables_dir, "\n")
cat("   - *_harmonised.csv: Harmonised data for each exposure\n")
cat("   - *_mr_results.csv: MR results for each exposure\n")
cat("   - *_heterogeneity.csv: Heterogeneity test results\n")
cat("   - *_pleiotropy.csv: Pleiotropy test results\n")
cat("   - *_steiger.csv: Steiger directionality test results\n")
cat("   - *_mr_presso.csv: MR-PRESSO results\n")
cat("2. Figures directory:", figures_dir, "\n")
cat("   - *_scatter_plot.png: Exposure-outcome scatter plot\n")
cat("   - *_forest_plot.png: Forest plot of SNP effects\n")
cat("   - *_funnel_plot.png: Funnel plot for symmetry\n")
cat("   - *_leaveoneout_plot.png: Leave-one-out sensitivity analysis\n")

cat("\nAnalysis completed at:", format(Sys.time()), "\n")
cat("Output directory:", base_dir, "\n")

# 关闭日志文件
sink()

cat("\nComplete MR analysis for eQTL exposures finished successfully!\n")
cat("All results saved to:", base_dir, "\n")
cat("Total exposures successfully analyzed:", sum(results_summary_eqtl$status == "Success"), "\n")

cat("\nSaving successfully analyzed eQTL genes...\n")

# 提取成功分析的基因
successful_genes <- results_summary_eqtl %>% 
  filter(status == "Success") %>% 
  pull(exposure) %>% 
  unique()

# 从暴露ID中提取Ensembl ID
successful_ensembl_ids <- gsub("^eqtl-a-", "", successful_genes)

# 获取基因符号（尝试从mapped_unique映射）
# 优先从 output_mapping 取展示名，缺失时再按 Ensembl 反查
if (exists("output_mapping") && !is.null(output_mapping) &&
    "ensembl_id" %in% colnames(output_mapping)) {
  display_col <- if ("gene_symbol" %in% colnames(output_mapping)) "gene_symbol" else "original_id"
  map_from_output <- setNames(as.character(output_mapping[[display_col]]), as.character(output_mapping$ensembl_id))
  gene_symbols <- ifelse(successful_ensembl_ids %in% names(map_from_output),
                         map_from_output[successful_ensembl_ids],
                         successful_ensembl_ids)
} else {
  symbol_map <- map_ensembl_to_symbol(successful_ensembl_ids)
  gene_symbols <- ifelse(
    !is.na(symbol_map[successful_ensembl_ids]) & nzchar(symbol_map[successful_ensembl_ids]),
    unname(symbol_map[successful_ensembl_ids]),
    successful_ensembl_ids
  )
}

# 保存
result_df <- data.frame(exposure = gene_symbols)
write.csv(result_df,
          file.path(base_dir, "09.eQTL_Successful_Genes.csv"),
          row.names = FALSE)

cat("✓ Saved", nrow(result_df), "successfully analyzed eQTL genes to 09.eQTL_Successful_Genes.csv\n")


# 8. 创建详细的异质性和多效性汇总表
cat("\n", rep("=", 60), "\n", sep="")
cat("CREATING DETAILED HETEROGENEITY AND PLEIOTROPY SUMMARY TABLE\n")
cat(rep("=", 60), "\n\n")

# 初始化详细汇总表
detailed_summary <- data.frame()

for (exp_id in names(all_analyses)) {
  analysis <- all_analyses[[exp_id]]
  
  if (!is.null(analysis$mr_results) && nrow(analysis$mr_results) > 0) {
    # 提取基因的Ensembl ID（从exposure_id中提取）
    ensembl_id <- gsub("^eqtl-a-", "", exp_id)
    
    # 尝试获取基因符号
    gene_symbol <- ensembl_id
    if (exists("output_mapping")) {
      # 在映射表中查找对应的原始ID
      matching_row <- output_mapping[grepl(paste0("^eqtl-a-", ensembl_id, "$"), 
                                           output_mapping$formatted_exposure_id), ]
      if (nrow(matching_row) > 0) {
        if ("gene_symbol" %in% colnames(matching_row) &&
            !is.na(matching_row$gene_symbol[1]) &&
            nzchar(matching_row$gene_symbol[1])) {
          gene_symbol <- matching_row$gene_symbol[1]
        } else {
          gene_symbol <- matching_row$original_id[1]
        }
      }
    }
    
    # 提取IVW和MR-Egger结果
    ivw_results <- analysis$mr_results[grepl("Inverse variance weighted", 
                                             analysis$mr_results$method, ignore.case = TRUE), ]
    egger_results <- analysis$mr_results[grepl("MR Egger", 
                                               analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取IVW异质性结果
    ivw_het <- NA
    ivw_het_df <- NA
    ivw_het_pval <- NA
    if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
      ivw_het_rows <- analysis$heterogeneity[grepl("Inverse variance weighted", 
                                                   analysis$heterogeneity$method, ignore.case = TRUE), ]
      if (nrow(ivw_het_rows) > 0) {
        ivw_het <- round(ivw_het_rows$Q[1], 4)
        ivw_het_df <- ifelse("Q_df" %in% colnames(ivw_het_rows), 
                             ivw_het_rows$Q_df[1], 
                             ivw_het_rows$nsnp[1] - 1)  # 如果Q_df不存在，用nsnp-1估计
        ivw_het_pval <- ifelse(!is.na(ivw_het_rows$Q_pval[1]), 
                               format.pval(ivw_het_rows$Q_pval[1], digits = 4, eps = 1e-99), 
                               "NA")
      }
    }
    
    # 提取MR-Egger异质性结果
    egger_het <- NA
    egger_het_df <- NA
    egger_het_pval <- NA
    if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
      egger_het_rows <- analysis$heterogeneity[grepl("MR Egger", 
                                                     analysis$heterogeneity$method, ignore.case = TRUE), ]
      if (nrow(egger_het_rows) > 0) {
        egger_het <- round(egger_het_rows$Q[1], 4)
        egger_het_df <- ifelse("Q_df" %in% colnames(egger_het_rows), 
                               egger_het_rows$Q_df[1], 
                               egger_het_rows$nsnp[1] - 2)  # MR-Egger自由度：nsnp-2
        egger_het_pval <- ifelse(!is.na(egger_het_rows$Q_pval[1]), 
                                 format.pval(egger_het_rows$Q_pval[1], digits = 4, eps = 1e-99), 
                                 "NA")
      }
    }
    
    # 提取MR-Egger截距检验结果
    egger_intercept <- NA
    egger_intercept_se <- NA
    egger_intercept_pval <- NA
    if (!is.null(analysis$pleiotropy) && nrow(analysis$pleiotropy) > 0) {
      egger_intercept <- round(analysis$pleiotropy$egger_intercept[1], 6)
      egger_intercept_se <- round(analysis$pleiotropy$se[1], 6)
      egger_intercept_pval <- ifelse(!is.na(analysis$pleiotropy$pval[1]), 
                                     format.pval(analysis$pleiotropy$pval[1], digits = 6, eps = 1e-99), 
                                     "NA")
    }
    
    # 1. 添加MR-Egger行
    if (nrow(egger_results) > 0) {
      egger_row <- data.frame(
        exposure = gene_symbol,
        outcome = outcome_id,
        Method = "MR Egger",
        Heterogeneity_Q = ifelse(!is.na(egger_het), as.character(egger_het), ""),
        Q_df = ifelse(!is.na(egger_het_df), as.character(egger_het_df), ""),
        Q_P.value = ifelse(!is.na(egger_het_pval), egger_het_pval, ""),
        Egger_intercept = ifelse(!is.na(egger_intercept), as.character(egger_intercept), ""),
        SE = ifelse(!is.na(egger_intercept_se), as.character(egger_intercept_se), ""),
        P.value = ifelse(!is.na(egger_intercept_pval), egger_intercept_pval, ""),
        stringsAsFactors = FALSE
      )
      detailed_summary <- rbind(detailed_summary, egger_row)
    }
    
    # 2. 添加IVW行（Egger截距列为空）
    if (nrow(ivw_results) > 0) {
      ivw_row <- data.frame(
        exposure = "",  # 与示例一样，IVW行exposure为空
        outcome = outcome_id,
        Method = "IVW",
        Heterogeneity_Q = ifelse(!is.na(ivw_het), as.character(ivw_het), ""),
        Q_df = ifelse(!is.na(ivw_het_df), as.character(ivw_het_df), ""),
        Q_P.value = ifelse(!is.na(ivw_het_pval), ivw_het_pval, ""),
        Egger_intercept = "",  # IVW方法没有截距
        SE = "",  # IVW方法没有截距标准误
        P.value = "",  # IVW方法没有截距P值
        stringsAsFactors = FALSE
      )
      detailed_summary <- rbind(detailed_summary, ivw_row)
    }
  }
}


cat("\n", rep("=", 60), "\n", sep="")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 60), "\n\n")

# 9. 创建并输出"有效基因的MR结果整合表格"
cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 9: CREATING INTEGRATED MR RESULTS TABLE FOR SUCCESSFUL GENES\n")
cat(rep("=", 60), "\n\n")

# 初始化整合结果表格
integrated_results <- data.frame()

# 计数器
genes_with_results <- 0
total_genes <- length(names(all_analyses))

cat("整合", total_genes, "个基因的MR结果...\n")

for (i in seq_along(names(all_analyses))) {
  exp_id <- names(all_analyses)[i]
  analysis <- all_analyses[[exp_id]]
  
  # 提取基因的Ensembl ID
  ensembl_id <- gsub("^eqtl-a-", "", exp_id)
  
  # 尝试获取基因符号
  gene_symbol <- ensembl_id
  if (exists("output_mapping")) {
    # 在映射表中查找对应的原始ID
    matching_row <- output_mapping[grepl(paste0("^eqtl-a-", ensembl_id, "$"), 
                                         output_mapping$formatted_exposure_id), ]
    if (nrow(matching_row) > 0) {
      if ("gene_symbol" %in% colnames(matching_row) &&
          !is.na(matching_row$gene_symbol[1]) &&
          nzchar(matching_row$gene_symbol[1])) {
        gene_symbol <- matching_row$gene_symbol[1]
      } else {
        gene_symbol <- matching_row$original_id[1]
      }
    }
  }
  
  # 检查是否有MR结果
  if (!is.null(analysis$mr_results) && nrow(analysis$mr_results) > 0) {
    genes_with_results <- genes_with_results + 1
    
    # 提取IVW结果
    ivw_rows <- analysis$mr_results[grepl("Inverse variance weighted", 
                                          analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取MR-Egger结果
    egger_rows <- analysis$mr_results[grepl("MR Egger", 
                                            analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取Weighted median结果
    wm_rows <- analysis$mr_results[grepl("Weighted median", 
                                         analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取Weighted mode结果
    wmode_rows <- analysis$mr_results[grepl("Weighted mode", 
                                            analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取Simple mode结果
    smode_rows <- analysis$mr_results[grepl("Simple mode", 
                                            analysis$mr_results$method, ignore.case = TRUE), ]
    
    # 提取SNP数量
    n_snps <- ifelse(!is.null(analysis$harmonised_data), nrow(analysis$harmonised_data), NA)
    
    # 创建IVW行
    if (nrow(ivw_rows) > 0) {
      ivw_row <- ivw_rows[1, ]
      
      # 提取OR和置信区间
      ivw_or <- ifelse("or" %in% colnames(ivw_row), ivw_row$or, NA)
      ivw_or_lci <- ifelse("or_lci95" %in% colnames(ivw_row), ivw_row$or_lci95, NA)
      ivw_or_uci <- ifelse("or_uci95" %in% colnames(ivw_row), ivw_row$or_uci95, NA)
      
      # 提取异质性结果
      ivw_het_q <- NA
      ivw_het_pval <- NA
      if (!is.null(analysis$heterogeneity) && nrow(analysis$heterogeneity) > 0) {
        ivw_het_rows <- analysis$heterogeneity[grepl("Inverse variance weighted", 
                                                     analysis$heterogeneity$method, ignore.case = TRUE), ]
        if (nrow(ivw_het_rows) > 0) {
          ivw_het_q <- round(ivw_het_rows$Q[1], 3)
          ivw_het_pval <- ifelse(!is.na(ivw_het_rows$Q_pval[1]), 
                                 format.pval(ivw_het_rows$Q_pval[1], digits = 3, eps = 1e-99), 
                                 "NA")
        }
      }
      
      # 提取MR-Egger截距检验
      egger_intercept <- NA
      egger_intercept_pval <- NA
      if (!is.null(analysis$pleiotropy) && nrow(analysis$pleiotropy) > 0) {
        egger_intercept <- round(analysis$pleiotropy$egger_intercept[1], 6)
        egger_intercept_pval <- ifelse(!is.na(analysis$pleiotropy$pval[1]), 
                                       format.pval(analysis$pleiotropy$pval[1], digits = 6, eps = 1e-99), 
                                       "NA")
      }
      
        # 提取Steiger结果
      steiger_correct <- "NA"
      steiger_pval <- "NA"
      if (!is.null(analysis$steiger) && nrow(analysis$steiger) > 0) {
        steiger_correct <- tryCatch({
          ifelse(!is.na(analysis$steiger$correct_causal_direction[1]),
                                  ifelse(analysis$steiger$correct_causal_direction[1], "TRUE", "FALSE"),
                                  "NA")
        }, error = function(e) "NA")
        steiger_pval <- tryCatch({
          ifelse(!is.na(analysis$steiger$steiger_pval[1]),
                               format.pval(analysis$steiger$steiger_pval[1], digits = 3, eps = 1e-99),
                               "NA")
        }, error = function(e) "NA")
      }

      # 预处理筛选变量（在data.frame之外计算）
      steiger_pval_num_val <- safe_parse_pval(steiger_pval)
      directionality_significant_val <- ifelse(
        !is.na(steiger_correct) && (steiger_correct == "TRUE" || steiger_correct == TRUE) &&
        !is.na(steiger_pval_num_val) && steiger_pval_num_val < 0.05,
        "YES", "NO"
      )

      egger_intercept_pval_num_val <- safe_parse_pval(egger_intercept_pval)
      pleiotropy_pass_val <- ifelse(
        !is.na(egger_intercept_pval_num_val) && egger_intercept_pval_num_val > 0.05,
        "YES", "NO"
      )

      # 创建主结果行
      main_row <- data.frame(
        # 基因信息
        gene_symbol = gene_symbol,
        ensembl_id = ensembl_id,
        exposure_id = exp_id,
        outcome_id = outcome_id,
        
        # SNP信息
        n_snps = n_snps,
        
        # IVW方法结果
        method = "IVW",
        ivw_beta = round(ivw_row$b, 4),
        ivw_se = round(ivw_row$se, 4),
        ivw_pval = ifelse(!is.na(ivw_row$pval), format.pval(ivw_row$pval, digits = 4, eps = 1e-99), "NA"),
        ivw_or = round(ivw_or, 3),
        ivw_or_lci = round(ivw_or_lci, 3),
        ivw_or_uci = round(ivw_or_uci, 3),
        ivw_or_95ci = paste0(round(ivw_or_lci, 3), "-", round(ivw_or_uci, 3)),
        
        # MR-Egger方法结果
        egger_beta = ifelse(nrow(egger_rows) > 0, round(egger_rows$b[1], 4), NA),
        egger_se = ifelse(nrow(egger_rows) > 0, round(egger_rows$se[1], 4), NA),
        egger_pval = ifelse(nrow(egger_rows) > 0 && !is.na(egger_rows$pval[1]), 
                            format.pval(egger_rows$pval[1], digits = 4, eps = 1e-99), "NA"),
        egger_or = ifelse(nrow(egger_rows) > 0 && "or" %in% colnames(egger_rows), 
                          round(egger_rows$or[1], 3), NA),
        
        # Weighted median方法结果
        wm_beta = ifelse(nrow(wm_rows) > 0, round(wm_rows$b[1], 4), NA),
        wm_se = ifelse(nrow(wm_rows) > 0, round(wm_rows$se[1], 4), NA),
        wm_pval = ifelse(nrow(wm_rows) > 0 && !is.na(wm_rows$pval[1]), 
                         format.pval(wm_rows$pval[1], digits = 4, eps = 1e-99), "NA"),
        wm_or = ifelse(nrow(wm_rows) > 0 && "or" %in% colnames(wm_rows), 
                       round(wm_rows$or[1], 3), NA),
        
        # Weighted mode方法结果
        wmode_beta = ifelse(nrow(wmode_rows) > 0, round(wmode_rows$b[1], 4), NA),
        wmode_se = ifelse(nrow(wmode_rows) > 0, round(wmode_rows$se[1], 4), NA),
        wmode_pval = ifelse(nrow(wmode_rows) > 0 && !is.na(wmode_rows$pval[1]), 
                            format.pval(wmode_rows$pval[1], digits = 4, eps = 1e-99), "NA"),
        wmode_or = ifelse(nrow(wmode_rows) > 0 && "or" %in% colnames(wmode_rows), 
                          round(wmode_rows$or[1], 3), NA),
        
        # Simple mode方法结果
        smode_beta = ifelse(nrow(smode_rows) > 0, round(smode_rows$b[1], 4), NA),
        smode_se = ifelse(nrow(smode_rows) > 0, round(smode_rows$se[1], 4), NA),
        smode_pval = ifelse(nrow(smode_rows) > 0 && !is.na(smode_rows$pval[1]), 
                            format.pval(smode_rows$pval[1], digits = 4, eps = 1e-99), "NA"),
        smode_or = ifelse(nrow(smode_rows) > 0 && "or" %in% colnames(smode_rows), 
                          round(smode_rows$or[1], 3), NA),
        
        # 异质性检验
        heterogeneity_q = ivw_het_q,
        heterogeneity_pval = ivw_het_pval,
        
        # 多效性检验
        egger_intercept = egger_intercept,
        egger_intercept_pval = egger_intercept_pval,
        
        # Steiger方向性检验
        steiger_correct_direction = steiger_correct,
        steiger_pval = steiger_pval,

        # 结果解释
        ivw_significant = ifelse(!is.na(ivw_row$pval) && ivw_row$pval < 0.05, "YES", "NO"),
        ivw_direction = ifelse(!is.na(ivw_row$b) && ivw_row$b > 0, "Risk", "Protective"),

        # 新增筛选条件列 (使用预先计算的值)
        directionality_significant = directionality_significant_val,
        pleiotropy_pass = pleiotropy_pass_val,

        stringsAsFactors = FALSE
      )
      
      integrated_results <- rbind(integrated_results, main_row)
    }
  }
  
  # 进度显示
  if (i %% 10 == 0 || i == total_genes) {
    cat("  已处理", i, "/", total_genes, "个基因\n")
  }
}

cat("✓ 成功整合", genes_with_results, "个基因的MR结果\n")

# ===== 从 all_analyses 提取 MR-PRESSO 结果并添加到 integrated_results =====
if (nrow(integrated_results) > 0 && exists("all_analyses")) {
  # 构建 gene_symbol -> exp_id 的反向映射
  gene_to_expid <- sapply(names(all_analyses), function(eid) {
    ensembl <- gsub("^eqtl-a-", "", eid)
    gs <- ensembl
    if (exists("output_mapping")) {
      mr <- output_mapping[grepl(paste0("^eqtl-a-", ensembl, "$"), output_mapping$formatted_exposure_id), ]
      if (nrow(mr) > 0) {
        if ("gene_symbol" %in% colnames(mr) &&
            !is.na(mr$gene_symbol[1]) &&
            nzchar(mr$gene_symbol[1])) {
          gs <- mr$gene_symbol[1]
        } else {
          gs <- mr$original_id[1]
        }
      }
    }
    return(gs)
  })
  names(gene_to_expid) <- names(all_analyses)

  # 提取PRESSO全局检验和异常值信息
  presso_data <- lapply(integrated_results$gene_symbol, function(gs) {
    exp_id <- names(gene_to_expid)[gene_to_expid == gs][1]
    if (is.na(exp_id) || is.null(exp_id)) return(data.frame(
      presso_global_p = NA_real_,
      presso_n_outliers = NA_integer_,
      presso_outlier_snps = NA_character_,
      presso_outlier_corrected_beta = NA_real_,
      presso_outlier_corrected_se = NA_real_,
      presso_outlier_corrected_p = NA_real_
    ))
    p <- all_analyses[[exp_id]]$mr_presso
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

if (nrow(integrated_results) > 0) {
  # 提取正确的outcome名称（从all_analyses的heterogeneity数据中获取完整GWAS名称）
  outcome_name <- "disease"
  first_exp_id <- names(all_analyses)[1]
  if (!is.null(all_analyses[[first_exp_id]]$heterogeneity) &&
      nrow(all_analyses[[first_exp_id]]$heterogeneity) > 0) {
    raw_outcome <- as.character(all_analyses[[first_exp_id]]$heterogeneity$outcome[1])
    # 格式如 "Nonalcoholic fatty liver disease || id:ebi-a-GCST90091033"
    if (!is.na(raw_outcome) && grepl("\\|\\|", raw_outcome)) {
      outcome_name <- gsub(" \\|\\|.*$", "", raw_outcome)
    } else if (!is.na(raw_outcome) && raw_outcome != "") {
      outcome_name <- raw_outcome
    }
  }
  cat("Outcome name:", outcome_name, "\n")

  # 按IVW P值排序
  integrated_results$ivw_pval_numeric <- vapply(integrated_results$ivw_pval, safe_parse_pval, numeric(1))
  integrated_results <- integrated_results[order(integrated_results$ivw_pval_numeric, na.last = NA), ]
  integrated_results$ivw_pval_numeric <- NULL

  # ===== 保存合并格式文件（00-07）=====
  # 00. 最完整初始结果（1行/基因，所有column）
  integrated_results$outcome <- outcome_name
  file_00 <- file.path(base_dir, "00.Complete_MR_Results.csv")
  write.csv(integrated_results, file_00, row.names = FALSE)
  cat("✓ 保存完整MR结果到:", file_00, "\n")

  # 01. Table 1格式（5行/基因：5种MR方法，包含所有column）
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
  cat("✓ 保存Table 1格式（5种方法）到:", file_01, "\n")

  # 02. 异质性检验汇总（2行/基因：IVW + MR Egger）+ MR-PRESSO global test
  het_list <- list()
  for (i in 1:nrow(integrated_results)) {
    row <- integrated_results[i, ]
    exp_id <- row$exposure_id
    gene_name <- row$gene_symbol

    # 获取MR-PRESSO global p
    presso_p <- NA
    if (!is.null(all_analyses[[exp_id]]$mr_presso) &&
        !is.null(all_analyses[[exp_id]]$mr_presso$global_test_p)) {
      presso_p <- all_analyses[[exp_id]]$mr_presso$global_test_p
    }

    # 获取Egger intercept pval
    egger_int_pval <- row$egger_intercept_pval

    # 从 all_analyses 读取完整的异质性结果（包含IVW和MR Egger两行）
    if (!is.null(all_analyses[[exp_id]]$heterogeneity) &&
        nrow(all_analyses[[exp_id]]$heterogeneity) > 0) {
      het_res <- all_analyses[[exp_id]]$heterogeneity
      # 按method遍历
      for (meth in c("Inverse variance weighted", "MR Egger")) {
        mrow <- het_res[grep(meth, het_res$method, ignore.case = TRUE), ]
        if (nrow(mrow) > 0) {
          mrow <- mrow[1, ]
          het_list[[length(het_list) + 1]] <- data.frame(
            Exposure = gene_name,
            Outcome = outcome_name,
            Method = ifelse(grepl("Inverse", meth, ignore.case = TRUE),
                          "Inverse variance weighted", "MR Egger"),
            Heterogeneity_Q = mrow$Q,
            Q_df = mrow$Q_df,
            Q_P_value = mrow$Q_pval,
            Egger_intercept_pval = egger_int_pval,
            MR_PRESSO_global_p = presso_p,
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      # fallback：从integrated_results取IVW结果
      het_list[[length(het_list) + 1]] <- data.frame(
        Exposure = gene_name,
        Outcome = outcome_name,
        Method = "Inverse variance weighted",
        Heterogeneity_Q = row$heterogeneity_q,
        Q_df = NA,
        Q_P_value = row$heterogeneity_pval,
        Egger_intercept_pval = egger_int_pval,
        MR_PRESSO_global_p = presso_p,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(het_list) > 0) {
    het_all <- do.call(rbind, het_list)
  } else {
    het_all <- data.frame(
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
  write.csv(het_all, het_file, row.names = FALSE)
  cat("✓ 保存异质性检验结果到:", het_file, "\n")

  # 03. 水平多效性检验汇总
  pleio_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = outcome_name,
    Egger_intercept = integrated_results$egger_intercept,
    Egger_intercept_pval = integrated_results$egger_intercept_pval,
    Pleiotropy_Significant = ifelse(
      vapply(integrated_results$egger_intercept_pval, safe_parse_pval, numeric(1)) < 0.05,
      "YES", "NO"
    ),
    stringsAsFactors = FALSE
  )
  pleio_file <- file.path(base_dir, "03.Pleiotropy_All_Genes.csv")
  write.csv(pleio_all, pleio_file, row.names = FALSE)
  cat("✓ 保存水平多效性检验结果到:", pleio_file, "\n")

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
  cat("✓ 保存MR-PRESSO结果到:", presso_file, "\n")

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
  cat("✓ 保存Steiger方向性分析结果到:", steiger_file, "\n")

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
  cat("✓ 保存三重筛选结果到:", tf_file, "\n")

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
    cat("✓ 保存因果靶点筛选结果到:", causal_file, "（共", nrow(causal_targets), "个）\n")
  }

  # ===== 保存合并格式文件（参考文本格式）=====

  # 11. 显示整合结果摘要
  cat("\n", rep("-", 120), "\n", sep="")
  cat("INTEGRATED MR RESULTS SUMMARY\n")
  cat("(Total: ", nrow(integrated_results), " genes)\n", sep="")
  cat(rep("-", 120), "\n", sep="")
  
  # 统计摘要
  cat("\nSTATISTICAL SUMMARY:\n")
  cat("• Total genes analyzed:", nrow(integrated_results), "\n")
  cat("• Genes with significant IVW (p < 0.05):", 
      sum(integrated_results$ivw_significant == "YES", na.rm = TRUE), "\n")
  cat("• Genes with protective effect (OR < 1):", 
      sum(integrated_results$ivw_or < 1, na.rm = TRUE), "\n")
  cat("• Genes with risk effect (OR > 1):", 
      sum(integrated_results$ivw_or > 1, na.rm = TRUE), "\n")
  cat("• Genes with heterogeneity (p < 0.05):", 
      sum(integrated_results$heterogeneity_pval != "NA" & 
            as.numeric(integrated_results$heterogeneity_pval) < 0.05, na.rm = TRUE), "\n")
  cat("• Genes with pleiotropy (p < 0.05):", 
      sum(integrated_results$egger_intercept_pval != "NA" & 
            as.numeric(integrated_results$egger_intercept_pval) < 0.05, na.rm = TRUE), "\n")
  cat("• Genes with correct direction (Steiger):", 
      sum(integrated_results$steiger_correct_direction == "TRUE", na.rm = TRUE), "\n")
  
  # 显示前20个基因
  cat("\nTOP 20 GENES (sorted by IVW P-value):\n")
  top_20 <- head(integrated_results, 20)
  
  display_table <- data.frame(
    Rank = 1:nrow(top_20),
    Gene = top_20$gene_symbol,
    Ensembl_ID = top_20$ensembl_id,
    n_SNPs = top_20$n_snps,
    IVW_OR = sprintf("%.3f", top_20$ivw_or),
    IVW_95CI = top_20$ivw_or_95ci,
    IVW_P = ifelse(top_20$ivw_pval == "NA", "NA", 
                   ifelse(as.numeric(top_20$ivw_pval) < 0.001, 
                          format(as.numeric(top_20$ivw_pval), scientific = TRUE, digits = 3),
                          sprintf("%.4f", as.numeric(top_20$ivw_pval)))),
    Significant = top_20$ivw_significant,
    Direction = top_20$ivw_direction
  )
  
  print(display_table, row.names = FALSE)
  
  # 显示最显著的前5个基因
  if (nrow(top_20) > 0) {
    cat("\nMOST SIGNIFICANT GENES:\n")
    for (i in 1:min(5, nrow(top_20))) {
      gene <- top_20[i, ]
      cat(sprintf("  %d. %s: OR = %.3f (95%% CI: %s), p = %s\n", 
                  i, 
                  gene$gene_symbol,
                  gene$ivw_or,
                  gene$ivw_or_95ci,
                  gene$ivw_pval))
    }
  }
  

} else {
  cat("⚠️ 警告: 没有MR结果可供整合\n")

  # 创建空的00文件
  empty_table <- data.frame(
    gene_symbol = character(0),
    ensembl_id = character(0),
    exposure_id = character(0),
    outcome_id = character(0),
    n_snps = numeric(0),
    method = character(0),
    stringsAsFactors = FALSE
  )
  file_00 <- file.path(base_dir, "00.Complete_MR_Results.csv")
  write.csv(empty_table, file_00, row.names = FALSE)
  cat("  ✓ 创建空的00文件:", file_00, "\n")
}

# 13. 最终输出文件列表
cat("\n", rep("=", 60), "\n", sep="")
cat("FINAL OUTPUT FILES - INTEGRATED RESULTS\n")
cat(rep("=", 60), "\n\n")

cat("📁 输出目录: ", base_dir, "\n\n")

cat("📄 主要结果文件 (00-07 合并格式):\n")
cat("   00. 00.Complete_MR_Results.csv       - 最完整初始结果（1行/基因）\n")
cat("   01. 01.MR_Results_All_Genes.csv     - Table 1格式（5行/基因：5种MR方法）\n")
cat("   02. 02.Heterogeneity_All_Genes.csv  - 异质性检验（IVW + MR Egger + PRESSO）\n")
cat("   03. 03.Pleiotropy_All_Genes.csv     - 水平多效性检验\n")
cat("   04. 04.MR_PRESSO_Results.csv        - MR-PRESSO结果\n")
cat("   05. 05.Steiger_Direction_Results.csv - Steiger方向性分析\n")
cat("   06. 06.Three_Filter_Results.csv     - 三重筛选结果\n")
cat("   07. 07.MR_res_gene.csv             - 因果靶点筛选结果\n")
cat("   08. 08.eQTL_ID_Mapping.csv                  - 基因ID映射表\n")
cat("   09. 09.eQTL_Successful_Genes.csv       - 成功分析的基因列表\n")
cat("   10. MR_Analysis_Log.txt             - 分析日志\n\n")

cat("📂 子目录:\n")
cat("   • tables/ - 每个基因的详细结果文件\n")
cat("   • figures/ - 所有可视化图表\n\n")

cat("📊 分析统计:\n")
cat("   • 总基因数: ", length(exposure_id_list), "\n")
cat("   • 成功分析: ", genes_with_results, "\n")
cat("   • 失败分析: ", length(exposure_id_list) - genes_with_results, "\n")
if (exists("integrated_results") && nrow(integrated_results) > 0) {
  cat("   • IVW显著基因 (p<0.05): ", sum(integrated_results$ivw_significant == "YES", na.rm = TRUE), "\n")
  cat("   • 方向性正确基因 (Steiger p<0.05 & direction=TRUE): ", sum(integrated_results$directionality_significant == "YES", na.rm = TRUE), "\n")
  cat("   • 无水平多效性基因 (Egger intercept p>0.05): ", sum(integrated_results$pleiotropy_pass == "YES", na.rm = TRUE), "\n")
  causal_count <- sum(integrated_results$ivw_significant == "YES" &
                      integrated_results$directionality_significant == "YES" &
                      integrated_results$pleiotropy_pass == "YES", na.rm = TRUE)
  cat("   • 三重筛选通过基因: ", causal_count, "\n")
  cat("   • 保护性基因: ", sum(as.numeric(integrated_results$ivw_or) < 1, na.rm = TRUE), "\n")
  cat("   • 风险性基因: ", sum(as.numeric(integrated_results$ivw_or) > 1, na.rm = TRUE), "\n")
}

cat("\n⏰ 分析完成时间: ", format(Sys.time()), "\n")
cat("📋 日志文件: ", log_file, "\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("ANALYSIS COMPLETE - ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep="")
cat(rep("=", 60), "\n")

# 10. Save session info + qs2 snapshot
run_tag <- if (!is.null(opt$timestamp) && nzchar(opt$timestamp)) {
  opt$timestamp
} else {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}
dataset_name <- basename(normalizePath(base_dir, mustWork = FALSE))
qs2_path <- file.path(base_dir, sprintf("%s_eqtl_node.qs2", dataset_name))
save_qs2_snapshot(list(
  integrated_results = integrated_results,
  results_summary = results_summary_eqtl,
  successful_genes = if (exists("result_df")) result_df else NULL
), qs2_path)
save_session_info_snapshot(file.path(base_dir, "logs"), dataset_name, run_tag)

# 11. Auto-generate Word report
cat("\n", rep("=", 60), "\n", sep="")
cat("STEP 10: AUTO-GENERATING MR REPORT\n")
cat(rep("=", 60), "\n\n")

# Determine script directory
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
                      "-t", "eqtl",
                      "-T", opt$timestamp)
  cat("Generating report...\n")
  system(report_cmd)
  cat("✓ Report saved to:", report_output, "\n")
} else {
  cat("Warning: generate_mr_report.R not found, skipping auto-report generation\n")
}