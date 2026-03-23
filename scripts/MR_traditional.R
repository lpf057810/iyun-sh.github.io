#!/usr/bin/env Rscript

# =============================================================================
# MR Analysis - Traditional MR
# =============================================================================
# Supports exposure type selection: exposure1, exposure2, or both
# =============================================================================

# 0. 辅助函数：解析可能带 "<" 前缀的 p 值（sync with eQTL/pQTL）
safe_parse_pval <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  x <- as.character(x)
  x <- gsub("^\\s*<\\s*", "", x)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0 || is.na(v)) return(NA_real_)
  v[1]
}

# Ensembl -> Gene Symbol（用于结果展示名）
map_ensembl_to_symbol <- function(ensembl_ids) {
  ids <- unique(sub("\\..*", "", as.character(ensembl_ids)))
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0) return(setNames(character(0), character(0)))

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

safe_pick <- function(df, col_name) {
  if (is.null(df) || nrow(df) == 0 || is.null(col_name) || !(col_name %in% colnames(df))) return(NA)
  val <- df[[col_name]][1]
  if (is.null(val) || length(val) == 0) return(NA)
  val
}

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

# Command line argument parsing
library(optparse)

option_list <- list(
  make_option(c("-t", "--type"), type = "character", default = "both",
              help = "Exposure type: exposure1, exposure2, or both [default %default]"),
  make_option(c("-o", "--output"), type = "character", default = "./output",
              help = "Output directory [default %default]"),
  make_option(c("-p", "--pval"), type = "character", default = "5e-8",
              help = "P-value threshold for SNP selection [default %default]"),
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
  make_option(c("-m", "--minsnps"), type = "character", default = "3",
              help = "Minimum number of SNPs [default %default]"),
  make_option(c("-T", "--token"), type = "character", default = NULL,
              help = "OpenGWAS JWT token (optional, overrides env/file token)"),
  make_option(c("--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run [default auto-generated]")
)

mr_config <- read_mr_config()
opt <- parse_args(OptionParser(option_list = option_list))
run_timestamp <- if (!is.null(opt$timestamp) && nzchar(opt$timestamp)) {
  opt$timestamp
} else {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

library(TwoSampleMR)
library(ieugwasr)
library(data.table)
library(dplyr)
library(ggplot2)
library(MRPRESSO)
library(ggsci)
library(showtext)
showtext_auto()
showtext_opts(dpi = 96)
font_add("Liberation Sans", regular = "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", bold = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf", italic = "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf", bolditalic = "/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf")

# Timestamp function
timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
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

# =============================================================================
# Token Management Functions (must be defined before use)
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

# Load token with priority: command line > env > .Renviron
token <- load_token(opt$token)
token_loaded <- !is.null(token) && nchar(token) > 0

cat("\n")
cat("================================================================================\n")
cat("  MR Traditional Analysis\n")
cat("  Started:", timestamp(), "\n")
cat("  Exposure type:", opt$type, "\n")
cat("  Output:", opt$output, "\n")
cat("================================================================================\n\n")

# Test token validity at startup
cat("[Token Check] Validating OpenGWAS token...\n")
token_test <- test_token_validity(token)
if (token_test$valid) {
  cat("[Token Check] OK - Token is valid\n")
  cat("[", timestamp(), "] JWT token loaded, starting with:", substr(token, 1, 30), "...\n")
} else {
  cat("[Token Check] FAILED -", token_test$message, "\n")
  cat("[Token Check] To update token, run: ./update_opengwas_token.sh <new_token>\n")
}

# Set Arial as default font
par(family = "Liberation Sans")

# Create output directory
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE)

# Save analysis config to result directory
config_file <- file.path(output_dir, "analysis_config.ini")
config_lines <- c(
  "[Analysis]",
  paste0("type = traditional"),
  paste0("exposure_type = ", opt$type),
  paste0("timestamp = ", opt$timestamp),
  "",
  "[Parameters]",
  paste0("pval_threshold = ", opt$pval),
  paste0("clump_kb = ", opt$kb),
  paste0("clump_r2 = ", opt$r2),
  paste0("fstat_threshold = ", opt$fstat),
  paste0("eaf_threshold = ", opt$eaf),
  paste0("outcome_gwas = ", opt$outcome),
  paste0("min_snps = ", opt$minsnps)
)
writeLines(config_lines, config_file)
cat("[Config] Parameters saved to:", config_file, "\n")

# Exposure configurations
exposure_configs <- list(
  exposure1 = list(
    id = c(
      "ukb-d-20480"
    ),
    name = "Self-harm"
  ),
  exposure2 = list(
    id = c(
      "ukb-e-20523_AFR",
      "ukb-e-20523_CSA"
    ),
    name = "Violence"
  )
)

# Determine which exposures to run
if (tolower(opt$type) == "exposure1") {
  exposures_to_run <- c("exposure1")
} else if (tolower(opt$type) == "exposure2") {
  exposures_to_run <- c("exposure2")
} else {
  exposures_to_run <- c("exposure1", "exposure2")
}

outcome_id <- "finn-b-TRAUMBRAIN_NONCONCUS"

# MR parameters
params <- list(
  p1 = as.numeric(opt$pval),
  r2 = as.numeric(opt$r2),
  kb = as.numeric(opt$kb),
  min_f = as.numeric(opt$fstat),
  min_snps = as.numeric(opt$minsnps)
)

cat("[", timestamp(), "] Outcome:", outcome_id, "\n")
cat("[", timestamp(), "] P-value threshold:", params$p1, "\n")

# Analysis summary data frame
analysis_summary <- data.frame(
  gene_symbol = character(),
  exposure = character(),
  outcome = character(),
  exposure_type = character(),
  status = character(),
  error_stage = character(),
  error_message = character(),
  n_snps = numeric(),
  f_stat_mean = numeric(),
  f_stat_min = numeric(),
  f_stat_max = numeric(),
  ivw_beta = numeric(),
  ivw_se = numeric(),
  ivw_pval = numeric(),
  ivw_or = numeric(),
  ivw_or_lci95 = numeric(),
  ivw_or_uci95 = numeric(),
  ivw_direction = character(),
  heterogeneity_significant = logical(),
  pleiotropy_significant = logical(),
  egger_intercept_pval = numeric(),
  steiger_pval = numeric(),
  steiger_correct_direction = character(),
  ivw_significant = character(),
  directionality_significant = character(),
  pleiotropy_pass = character(),
  three_filter_pass = character(),
  presso_global_p = numeric(),
  presso_n_outliers = numeric()
)

# Define theme with Arial font
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

# Loop through each exposure type
for (exp_type in exposures_to_run) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Processing:", exposure_configs[[exp_type]]$name, "\n")
  cat(rep("=", 70), "\n\n")

  exposure_id_list <- exposure_configs[[exp_type]]$id
  exp_idx <- 0  # Per-type index counter

  for (exp_id in exposure_id_list) {
    exp_idx <- exp_idx + 1
    cat("\n[", timestamp(), "] ---", exp_id, "---\n")

    # 1. Extract instrumental variables
    exposure_ivs <- tryCatch({
      retry_with_token_check(
        extract_instruments(outcomes = exp_id, p1 = params$p1, r2 = params$r2, kb = params$kb),
        max_retries = 3,
        on_token_expired = function() {
          cat("[", timestamp(), "] FATAL: Token expired during IV extraction!\n")
        }
      )
    }, error = function(e) {
      cat("[", timestamp(), "] IV extraction failed:", e$message, "\n")
      NULL
    })

    if (is.null(exposure_ivs) || nrow(exposure_ivs) == 0) {
      cat("[", timestamp(), "] No SNPs extracted\n")
      analysis_summary <- rbind(analysis_summary, data.frame(
        gene_symbol = exp_id, exposure = exp_id, outcome = outcome_id, exposure_type = exp_type, status = "Failed",
        error_stage = "IV_Extraction", error_message = "No SNPs extracted",
        n_snps = 0, f_stat_mean = NA, f_stat_min = NA, f_stat_max = NA,
        ivw_beta = NA, ivw_se = NA, ivw_pval = NA, ivw_or = NA,
        ivw_or_lci95 = NA, ivw_or_uci95 = NA, ivw_direction = NA,
        heterogeneity_significant = NA, pleiotropy_significant = NA, egger_intercept_pval = NA,
        steiger_pval = NA, steiger_correct_direction = NA,
        ivw_significant = NA, directionality_significant = NA, pleiotropy_pass = NA, three_filter_pass = NA,
        presso_global_p = NA, presso_n_outliers = NA
      ))
      next
    }

    cat("[", timestamp(), "] SNPs extracted:", nrow(exposure_ivs), "\n")

    # 2. Extract outcome data
    snps <- exposure_ivs$SNP
    outcome_data <- tryCatch({
      retry_with_token_check(
        extract_outcome_data(snps = snps, outcomes = outcome_id),
        max_retries = 3,
        on_token_expired = function() {
          cat("[", timestamp(), "] FATAL: Token expired during outcome extraction!\n")
        }
      )
    }, error = function(e) {
      cat("[", timestamp(), "] Outcome extraction failed:", e$message, "\n")
      NULL
    })

    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      cat("[", timestamp(), "] No outcome data\n")
      analysis_summary <- rbind(analysis_summary, data.frame(
        gene_symbol = exp_id, exposure = exp_id, outcome = outcome_id, exposure_type = exp_type, status = "Failed",
        error_stage = "Outcome_Extraction", error_message = "No outcome data",
        n_snps = nrow(exposure_ivs), f_stat_mean = NA, f_stat_min = NA, f_stat_max = NA,
        ivw_beta = NA, ivw_se = NA, ivw_pval = NA, ivw_or = NA,
        ivw_or_lci95 = NA, ivw_or_uci95 = NA, ivw_direction = NA,
        heterogeneity_significant = NA, pleiotropy_significant = NA, egger_intercept_pval = NA,
        steiger_pval = NA, steiger_correct_direction = NA,
        ivw_significant = NA, directionality_significant = NA, pleiotropy_pass = NA, three_filter_pass = NA,
        presso_global_p = NA, presso_n_outliers = NA
      ))
      next
    }

    # 3. Harmonise data
    harmonised <- harmonise_data(exposure_ivs, outcome_data, action = 2)

    # 3.1 Remove SNPs associated with outcome (exclusivity assumption, p < 1e-5)
    if ("pval.outcome" %in% colnames(harmonised)) {
      n_before <- nrow(harmonised)
      harmonised <- harmonised[harmonised$pval.outcome >= 1e-5, ]
      n_after <- nrow(harmonised)
      if (n_before > n_after) {
        cat("[", timestamp(), "] Removed", n_before - n_after, "SNPs associated with outcome (p < 1e-5)\n")
      }
      # Check if enough SNPs remain
      if (nrow(harmonised) < params$min_snps) {
        cat("[", timestamp(), "] Not enough SNPs after outcome-associated filter:", nrow(harmonised), "\n")
        analysis_summary <- rbind(analysis_summary, data.frame(
          exposure = exp_id, outcome = outcome_id, exposure_type = exp_type, status = "Failed",
          error_stage = "OutcomeFilter", error_message = paste0("Not enough SNPs after removing outcome-associated (n=", nrow(harmonised), ")"),
          n_snps = nrow(harmonised), f_stat_mean = NA, f_stat_min = NA, f_stat_max = NA,
          ivw_beta = NA, ivw_se = NA, ivw_pval = NA, ivw_or = NA,
          ivw_or_lci95 = NA, ivw_or_uci95 = NA, ivw_direction = NA,
          heterogeneity_significant = NA, pleiotropy_significant = NA, egger_intercept_pval = NA,
          steiger_pval = NA, steiger_correct_direction = NA,
          ivw_significant = NA, directionality_significant = NA, pleiotropy_pass = NA, three_filter_pass = NA,
          presso_global_p = NA, presso_n_outliers = NA
        ))
        next
      }
    }

    # 4. Calculate F-statistics
    if ("beta.exposure" %in% colnames(harmonised) && "se.exposure" %in% colnames(harmonised)) {
      harmonised$F.exposure <- (harmonised$beta.exposure / harmonised$se.exposure)^2
    }

    f_stats <- harmonised$F.exposure
    f_stat_mean <- mean(f_stats, na.rm = TRUE)
    f_stat_min <- min(f_stats, na.rm = TRUE)
    f_stat_max <- max(f_stats, na.rm = TRUE)

    cat("[", timestamp(), "] F-statistics: mean=", round(f_stat_mean, 2),
        ", min=", round(f_stat_min, 2), ", max=", round(f_stat_max, 2), "\n")

    # 5. Check for weak instruments
    if (any(f_stats < params$min_f, na.rm = TRUE)) {
      n_weak <- sum(f_stats < params$min_f, na.rm = TRUE)
      cat("[", timestamp(), "] Warning:", n_weak, " weak instruments (F < 10)\n")
    }

    # 6. Check SNP count
    if (nrow(harmonised) < params$min_snps) {
      cat("[", timestamp(), "] Not enough SNPs after harmonisation:", nrow(harmonised), "\n")
      analysis_summary <- rbind(analysis_summary, data.frame(
        exposure = exp_id, outcome = outcome_id, exposure_type = exp_type, status = "Failed",
        error_stage = "Harmonisation", error_message = paste0("Not enough SNPs (n=", nrow(harmonised), ")"),
        n_snps = nrow(harmonised), f_stat_mean = f_stat_mean, f_stat_min = f_stat_min, f_stat_max = f_stat_max,
        ivw_beta = NA, ivw_se = NA, ivw_pval = NA, ivw_or = NA,
        ivw_or_lci95 = NA, ivw_or_uci95 = NA, ivw_direction = NA,
        heterogeneity_significant = NA, pleiotropy_significant = NA, egger_intercept_pval = NA,
        steiger_pval = NA, steiger_correct_direction = NA,
        ivw_significant = NA, directionality_significant = NA, pleiotropy_pass = NA, three_filter_pass = NA,
        presso_global_p = NA, presso_n_outliers = NA
      ))
      next
    }

    # 7. Heterogeneity test
    het <- tryCatch({
      het_result <- mr_heterogeneity(harmonised)
      if (!is.null(het_result)) {
        cat("[", timestamp(), "] Heterogeneity methods:", paste(unique(het_result$method), collapse = ", "), "\n")
      }
      het_result
    }, error = function(e) {
      cat("[", timestamp(), "] Heterogeneity test failed:", e$message, "\n")
      return(NULL)
    })

    # Select MR method based on heterogeneity
    het_qvalue <- NA
    if (!is.null(het) && nrow(het) > 0) {
      print(het[, c("method", "Q", "Q_pval")])
      ivw_idx <- grep("Inverse variance", het$method, ignore.case = TRUE)
      if (length(ivw_idx) > 0) {
        het_qvalue <- het$Q_pval[ivw_idx[1]]
      }
    }

    if (!is.na(het_qvalue) && het_qvalue < 0.05) {
      cat("[", timestamp(), "] Heterogeneity detected → Using RANDOM EFFECTS IVW\n")
      method_list <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median",
                       "mr_weighted_mode", "mr_simple_mode")
    } else {
      cat("[", timestamp(), "] No significant heterogeneity → Using FIXED EFFECTS IVW\n")
      method_list <- c("mr_ivw_fe", "mr_egger_regression", "mr_weighted_median",
                       "mr_weighted_mode", "mr_simple_mode")
    }

    # 8. Run MR analysis
    mr_results <- mr(harmonised, method_list = method_list)

    # 9. Calculate OR and CI
    mr_results <- generate_odds_ratios(mr_results)
    mr_results$outcome <- outcome_id

    # 10. Pleiotropy test
    pleio <- mr_pleiotropy_test(harmonised)

    # 11. MR-PRESSO analysis
    presso <- tryCatch({
      run_mr_presso(harmonised, NbDistribution = 1000, SignifThreshold = 0.05)
    }, error = function(e) {
      cat("[", timestamp(), "] MR-PRESSO failed:", e$message, "\n")
      NULL
    })

    # 12. Steiger test
    steiger <- tryCatch({
      mr_steiger(harmonised)
    }, error = function(e) {
      cat("[", timestamp(), "] Steiger test failed:", e$message, "\n")
      NULL
    })

    # 13. Save results
    safe_name <- gsub("-", "_", exp_id)
    tables_dir <- file.path(output_dir, "tables", exp_type)
    figures_dir <- file.path(output_dir, "figures", exp_type)
    dir.create(tables_dir, recursive = TRUE)
    dir.create(figures_dir, recursive = TRUE)
    dir.create(file.path(figures_dir, "01_scatter"), recursive = TRUE)
    dir.create(file.path(figures_dir, "02_forest"), recursive = TRUE)
    dir.create(file.path(figures_dir, "03_funnel"), recursive = TRUE)
    dir.create(file.path(figures_dir, "04_leaveoneout"), recursive = TRUE)

    write.csv(mr_results, file.path(tables_dir, paste0(safe_name, "_mr_results.csv")), row.names = FALSE)
    write.csv(harmonised, file.path(tables_dir, paste0(safe_name, "_harmonised.csv")), row.names = FALSE)
    write.csv(het, file.path(tables_dir, paste0(safe_name, "_heterogeneity.csv")), row.names = FALSE)
    write.csv(pleio, file.path(tables_dir, paste0(safe_name, "_pleiotropy.csv")), row.names = FALSE)
    if (!is.null(steiger)) {
      write.csv(steiger, file.path(tables_dir, paste0(safe_name, "_steiger.csv")), row.names = FALSE)
    }

    # Save MR-PRESSO results
    if (!is.null(presso) && length(presso) > 0) {
      tryCatch({
        presso_res <- presso[[1]]
        if (!is.null(presso_res) && is.data.frame(presso_res)) {
          write.csv(presso_res, file.path(tables_dir, paste0(safe_name, "_mr_presso.csv")), row.names = FALSE)
          cat("[", timestamp(), "] MR-PRESSO results saved\n")
        }
      }, error = function(e) {
        cat("[", timestamp(), "] MR-PRESSO save failed:", e$message, "\n")
      })
    }

    # 14. Check significance
    ivw <- mr_results[mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw) > 0) {
      pval <- ivw$pval
      or <- ivw$or
      or_lci <- ivw$or_lci95
      or_uci <- ivw$or_uci95
      cat("[", timestamp(), "] IVW: b=", round(ivw$b, 4), ", OR=", round(or, 3),
          ", 95%CI=[", round(or_lci, 3), ",", round(or_uci, 3), "], p=", round(pval, 4), "\n")

      status <- ifelse(!is.na(pval) && pval < 0.05, "Success (p<0.05)", "Success")
    } else {
      pval <- NA
      or <- NA
      or_lci <- NA
      or_uci <- NA
      status <- "Success (no IVW)"
    }

    # Check heterogeneity and pleiotropy
    het_sig <- if (nrow(het) > 0) {
      any(vapply(het$pval, safe_parse_pval, numeric(1)) < 0.05, na.rm = TRUE)
    } else FALSE
    pleio_sig <- if (nrow(pleio) > 0) {
      any(vapply(pleio$pval, safe_parse_pval, numeric(1)) < 0.05, na.rm = TRUE)
    } else FALSE

    # Get Egger intercept p-value for pleiotropy assessment
    egger_pval <- NA
    if (!is.null(pleio) && nrow(pleio) > 0) {
      egger_row <- pleio[pleio$method == "MR Egger", ]
      if (nrow(egger_row) > 0) {
        egger_pval <- egger_row$pval[1]
      }
    }

    # Get Steiger results
    steiger_pval <- NA
    steiger_correct <- NA
    if (!is.null(steiger) && nrow(steiger) > 0) {
      steiger_pval <- steiger$steiger_pval[1]
      # correct_causal_direction is a logical (TRUE/FALSE), convert to string "TRUE"/"FALSE"
      steiger_correct <- ifelse(!is.na(steiger$correct_causal_direction[1]),
                                ifelse(steiger$correct_causal_direction[1], "TRUE", "FALSE"),
                                "NA")
    }

    # Calculate direction
    ivw_direction <- NA
    if (!is.na(pval) && !is.na(ivw$b)) {
      ivw_direction <- ifelse(ivw$b > 0, "positive", "negative")
    }

    # Three-filter logic (三重筛选)
    # Filter 1: IVW significant (p < 0.05)
    ivw_significant <- ifelse(!is.na(pval) && pval < 0.05, "YES", "NO")

    # Filter 2: Directionality significant (Steiger p < 0.05 AND correct direction)
    # 处理 "<1e-99" 格式的p值: 先去掉"<"符号再转数值
    steiger_pval_num <- safe_parse_pval(steiger_pval)
    # 兼容字符串 "TRUE" 和逻辑值 TRUE
    directionality_significant <- ifelse(
      !is.na(steiger_correct) && (steiger_correct == "TRUE" || steiger_correct == TRUE) &&
      !is.na(steiger_pval_num) && steiger_pval_num < 0.05,
      "YES", "NO"
    )

    # Filter 3: Pleiotropy pass (Egger intercept p > 0.05, i.e., no significant pleiotropy)
    egger_pval_num <- safe_parse_pval(egger_pval)
    pleiotropy_pass <- ifelse(
      !is.na(egger_pval_num) && egger_pval_num > 0.05,
      "YES", "NO"
    )

    # Combined three-filter: pass all three
    three_filter_pass <- ifelse(
      ivw_significant == "YES" & directionality_significant == "YES" & pleiotropy_pass == "YES",
      "YES", "NO"
    )

    cat("[", timestamp(), "] Three-filter: IVW=", ivw_significant, ", Direction=", directionality_significant, ", Pleiotropy=", pleiotropy_pass, ", Combined=", three_filter_pass, "\n")

    # Check MR-PRESSO results
    presso_pval <- NA
    presso_outliers <- NA
    if (!is.null(presso) && length(presso) > 0) {
      tryCatch({
        presso_res <- presso[[1]]
        if (!is.null(presso_res) && is.data.frame(presso_res) && nrow(presso_res) > 0) {
          if ("MR.PRESSO.results.Global.Test.Pvalue" %in% colnames(presso_res)) {
            presso_pval <- safe_parse_pval(presso_res$MR.PRESSO.results.Global.Test.Pvalue[1])
          } else if ("global_test_p" %in% colnames(presso_res)) {
            presso_pval <- safe_parse_pval(presso_res$global_test_p[1])
          }
          if ("n_outliers" %in% colnames(presso_res)) {
            presso_outliers <- suppressWarnings(as.numeric(presso_res$n_outliers[1]))
          }
        }
      }, error = function(e) {})
    }

    analysis_summary <- rbind(analysis_summary, data.frame(
      gene_symbol = exp_id,
      exposure = exp_id,
      outcome = outcome_id,
      exposure_type = exp_type,
      status = status,
      error_stage = "None",
      error_message = NA,
      n_snps = nrow(harmonised),
      f_stat_mean = f_stat_mean,
      f_stat_min = f_stat_min,
      f_stat_max = f_stat_max,
      ivw_beta = if (nrow(ivw) > 0) ivw$b else NA,
      ivw_se = if (nrow(ivw) > 0) ivw$se else NA,
      ivw_pval = pval,
      ivw_or = or,
      ivw_or_lci95 = or_lci,
      ivw_or_uci95 = or_uci,
      ivw_direction = ivw_direction,
      heterogeneity_significant = het_sig,
      pleiotropy_significant = pleio_sig,
      egger_intercept_pval = egger_pval,
      steiger_pval = steiger_pval,
      steiger_correct_direction = steiger_correct,
      ivw_significant = ivw_significant,
      directionality_significant = directionality_significant,
      pleiotropy_pass = pleiotropy_pass,
      three_filter_pass = three_filter_pass,
      presso_global_p = presso_pval,
      presso_n_outliers = presso_outliers
    ))

    # 15. Generate plots
    exp_name_simple <- gsub("-", "_", exp_id)

    # Scatter plot
    tryCatch({
      p <- mr_scatter_plot(mr_results, harmonised)
      if (!is.null(p) && length(p) > 0) {
        p1 <- p[[1]] +
          scale_color_npg(palette = "nrc") +
          theme_mr(base_size = 14, legend.position = "top") +
          ggtitle("MR Scatter Plot") +
          theme(axis.title.y = element_text(size = 10, hjust = 0.5, family = "Liberation Sans")) +
          labs(y = "Effect on TBI")
        num_prefix <- sprintf("%02d", exp_idx)
        showtext_begin()
        ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", exp_name_simple, ".scatter.png")),
               plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
        showtext_end()
        showtext_begin()
        ggsave(file.path(figures_dir, "01_scatter", paste0(num_prefix, ".", exp_name_simple, ".scatter.pdf")),
               plot = p1, width = 10, height = 9, device = "pdf")
        showtext_end()
        cat("[", timestamp(), "] Scatter plot saved\n")
      }
    }, error = function(e) cat("[", timestamp(), "] Scatter plot failed:", e$message, "\n"))

    # Forest plot
    tryCatch({
      res_single <- mr_singlesnp(harmonised)
      p <- mr_forest_plot(res_single)
      if (!is.null(p) && length(p) > 0) {
        p1 <- p[[1]] +
          scale_color_npg(palette = "nrc") +
          theme_mr(base_size = 12) +
          ggtitle("MR Forest Plot") +
          theme(axis.title.y = element_text(size = 10, hjust = 0.5, family = "Liberation Sans")) +
          labs(y = "Non-concussion TBI")
        num_prefix <- sprintf("%02d", exp_idx)
        showtext_begin()
        ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", exp_name_simple, ".forest.png")),
               plot = p1, width = 16, height = max(10, nrow(harmonised) * 0.5), dpi = 300, bg = "white")
        showtext_end()
        showtext_begin()
        ggsave(file.path(figures_dir, "02_forest", paste0(num_prefix, ".", exp_name_simple, ".forest.pdf")),
               plot = p1, width = 16, height = max(10, nrow(harmonised) * 0.5), device = "pdf")
        showtext_end()
        # 保存rsID.csv
        if (!is.null(res_single) && nrow(res_single) > 0 && "SNP" %in% colnames(res_single)) {
          rsid_df <- data.frame(SNP = res_single$SNP, stringsAsFactors = FALSE)
          write.csv(rsid_df, file.path(figures_dir, "02_forest", paste0(num_prefix, ".", exp_name_simple, ".rsID.csv")),
                    row.names = FALSE)
        }
        cat("[", timestamp(), "] Forest plot saved\n")
      }
    }, error = function(e) cat("[", timestamp(), "] Forest plot failed:", e$message, "\n"))

    # Funnel plot
    tryCatch({
      res_single <- mr_singlesnp(harmonised)
      p <- mr_funnel_plot(res_single)
      if (!is.null(p) && length(p) > 0) {
        p1 <- p[[1]] +
          scale_color_npg(palette = "nrc") +
          theme_mr(base_size = 14, legend.position = "top") +
          ggtitle("MR Funnel Plot")
        num_prefix <- sprintf("%02d", exp_idx)
        showtext_begin()
        ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", exp_name_simple, ".funnel.png")),
               plot = p1, width = 10, height = 9, dpi = 300, bg = "white")
        showtext_end()
        showtext_begin()
        ggsave(file.path(figures_dir, "03_funnel", paste0(num_prefix, ".", exp_name_simple, ".funnel.pdf")),
               plot = p1, width = 10, height = 9, device = "pdf")
        showtext_end()
        cat("[", timestamp(), "] Funnel plot saved\n")
      }
    }, error = function(e) cat("[", timestamp(), "] Funnel plot failed:", e$message, "\n"))

    # Leave-one-out plot
    tryCatch({
      res_loo <- mr_leaveoneout(harmonised)
      p <- mr_leaveoneout_plot(res_loo)
      if (!is.null(p) && length(p) > 0) {
        p1 <- p[[1]] +
          scale_color_npg(palette = "nrc") +
          theme_mr(base_size = 12) +
          ggtitle("MR Leave-One-Out Plot") +
          theme(axis.title.y = element_text(size = 10, hjust = 0.5, family = "Liberation Sans")) +
          labs(y = "Non-concussion TBI")
        num_prefix <- sprintf("%02d", exp_idx)
        showtext_begin()
        ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", exp_name_simple, ".leaveoneout.png")),
               plot = p1, width = 16, height = max(10, nrow(harmonised) * 0.4), dpi = 300, bg = "white")
        showtext_end()
        showtext_begin()
        ggsave(file.path(figures_dir, "04_leaveoneout", paste0(num_prefix, ".", exp_name_simple, ".leaveoneout.pdf")),
               plot = p1, width = 16, height = max(10, nrow(harmonised) * 0.4), device = "pdf")
        showtext_end()
        cat("[", timestamp(), "] Leave-one-out plot saved\n")
      }
    }, error = function(e) cat("[", timestamp(), "] Leave-one-out plot failed:", e$message, "\n"))
  }
}

# Save summary
write.csv(analysis_summary, file.path(output_dir, "tables", "analysis_summary.csv"), row.names = FALSE)

# Rename columns to match eQTL/pQTL format
integrated_results <- analysis_summary
if (nrow(integrated_results) > 0) {
  # Add 95% CI for OR if not already present
  if (!"ivw_or_95ci" %in% colnames(integrated_results)) {
    integrated_results$ivw_or_95ci <- paste0("(", round(integrated_results$ivw_or_lci95, 3), "-", round(integrated_results$ivw_or_uci95, 3), ")")
  }
  integrated_results$ivw_pval_numeric <- vapply(integrated_results$ivw_pval, safe_parse_pval, numeric(1))
  integrated_results <- integrated_results[order(integrated_results$ivw_pval_numeric, na.last = NA), ]
  integrated_results$ivw_pval_numeric <- NULL
}

cat("\n")

# ===== 保存合并格式文件（00-07）=====
if (nrow(integrated_results) > 0) {
  cat("[", timestamp(), "] Saving merged format files (00-07)...\n")

  # 00. 最完整初始结果（1行/暴露，所有column）
  file_00 <- file.path(output_dir, "00.Complete_MR_Results.csv")
  write.csv(integrated_results, file_00, row.names = FALSE)
  cat("[", timestamp(), "] 00.Complete_MR_Results.csv saved\n")

  # 01. Table 1格式（5行/暴露：5种MR方法，从per-exposure文件读取）
  mr_table1_list <- list()
  mr_files <- list.files(file.path(output_dir, "tables"), pattern = "_mr_results\\.csv$", recursive = TRUE, full.names = TRUE)
  for (f in mr_files) {
    mr_dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(mr_dt) || nrow(mr_dt) == 0) next
    safe_name <- gsub("_mr_results\\.csv$", "", basename(f))
    exp_name <- gsub("_", " ", safe_name)
    outcome_raw <- as.character(mr_dt$outcome[1])
    if (is.na(outcome_raw) || outcome_raw == "") outcome_raw <- as.character(mr_dt$id.outcome[1])
    if (grepl("\\|\\|", outcome_raw)) outcome_name <- gsub(" \\|\\|.*$", "", outcome_raw)
    else outcome_name <- outcome_raw
    for (i in 1:nrow(mr_dt)) {
      mr_table1_list[[length(mr_table1_list) + 1]] <- data.frame(
        Exposure = exp_name,
        Outcome = outcome_name,
        Method = as.character(mr_dt$method[i]),
        n_SNPs = mr_dt$nsnp[i],
        beta = mr_dt$b[i],
        se = mr_dt$se[i],
        pval = mr_dt$pval[i],
        or = mr_dt$or[i],
        or_lci = mr_dt$or_lci95[i],
        or_uci = mr_dt$or_uci95[i]
      )
    }
  }
  if (length(mr_table1_list) > 0) {
    mr_table1_df <- do.call(rbind, mr_table1_list)
    all_mr_file <- file.path(output_dir, "01.MR_Results_All_Genes.csv")
    write.csv(mr_table1_df, all_mr_file, row.names = FALSE)
    cat("[", timestamp(), "] 01.MR_Results_All_Genes.csv saved (5 methods per exposure)\n")
  }

  # 02. 异质性检验汇总（与eQTL/pQTL对齐：IVW + MR Egger）+ MR-PRESSO
  # 构建 exposure -> presso_global_p / egger_intercept_pval 映射
  presso_map <- setNames(integrated_results$presso_global_p, integrated_results$gene_symbol)
  egger_map <- setNames(integrated_results$egger_intercept_pval, integrated_results$gene_symbol)
  het_list <- list()
  het_files <- list.files(file.path(output_dir, "tables"), pattern = "_heterogeneity\\.csv$", recursive = TRUE, full.names = TRUE)
  for (f in het_files) {
    het_dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(het_dt) || nrow(het_dt) == 0) next
    safe_name <- gsub("_heterogeneity\\.csv$", "", basename(f))
    exp_candidates <- c(gsub("_", " ", safe_name), gsub("_", "-", safe_name), safe_name)
    exp_key <- exp_candidates[exp_candidates %in% names(presso_map)][1]
    if (is.na(exp_key) || is.null(exp_key) || !nzchar(exp_key)) exp_key <- exp_candidates[1]
    outcome_raw <- as.character(het_dt$outcome[1])
    if (grepl("\\|\\|", outcome_raw)) outcome_name <- gsub(" \\|\\|.*$", "", outcome_raw)
    else outcome_name <- outcome_raw
    presso_p <- if (exp_key %in% names(presso_map)) safe_parse_pval(presso_map[[exp_key]]) else NA_real_
    egger_int_pval <- if (exp_key %in% names(egger_map)) safe_parse_pval(egger_map[[exp_key]]) else NA_real_
    for (meth in c("Inverse variance weighted", "MR Egger")) {
      mrow <- het_dt[grep(meth, het_dt$method, ignore.case = TRUE), ]
      if (nrow(mrow) == 0) next
      mrow <- mrow[1, ]
      het_list[[length(het_list) + 1]] <- data.frame(
        Exposure = exp_key,
        Outcome = outcome_name,
        Method = ifelse(grepl("Inverse", meth, ignore.case = TRUE), "Inverse variance weighted", "MR Egger"),
        Heterogeneity_Q = safe_pick(mrow, "Q"),
        Q_df = safe_pick(mrow, "Q_df"),
        Q_P_value = safe_pick(mrow, "Q_pval"),
        Egger_intercept_pval = egger_int_pval,
        MR_PRESSO_global_p = presso_p,
        stringsAsFactors = FALSE
      )
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
  het_file <- file.path(output_dir, "02.Heterogeneity_All_Genes.csv")
  write.csv(het_df, het_file, row.names = FALSE)
  cat("[", timestamp(), "] 02.Heterogeneity_All_Genes.csv saved\n")

  # 03. 水平多效性检验汇总（与eQTL/pQTL对齐）
  pleio_df <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = integrated_results$outcome,
    Egger_intercept = if ("egger_intercept" %in% colnames(integrated_results)) integrated_results$egger_intercept else NA_real_,
    Egger_intercept_pval = integrated_results$egger_intercept_pval,
    Pleiotropy_Significant = ifelse(
      vapply(integrated_results$egger_intercept_pval, safe_parse_pval, numeric(1)) < 0.05,
      "YES", "NO"
    ),
    stringsAsFactors = FALSE
  )
  pleio_file <- file.path(output_dir, "03.Pleiotropy_All_Genes.csv")
  write.csv(pleio_df, pleio_file, row.names = FALSE)
  cat("[", timestamp(), "] 03.Pleiotropy_All_Genes.csv saved\n")

  # 04. MR-PRESSO汇总
  presso_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = integrated_results$outcome,
    MR_PRESSO_global_p = integrated_results$presso_global_p,
    n_outliers = integrated_results$presso_n_outliers,
    outlier_snps = if ("presso_outlier_snps" %in% colnames(integrated_results)) integrated_results$presso_outlier_snps else NA_character_,
    outlier_corrected_beta = if ("presso_outlier_corrected_beta" %in% colnames(integrated_results)) integrated_results$presso_outlier_corrected_beta else NA_real_,
    outlier_corrected_se = if ("presso_outlier_corrected_se" %in% colnames(integrated_results)) integrated_results$presso_outlier_corrected_se else NA_real_,
    outlier_corrected_p = if ("presso_outlier_corrected_p" %in% colnames(integrated_results)) integrated_results$presso_outlier_corrected_p else NA_real_,
    stringsAsFactors = FALSE
  )
  presso_file <- file.path(output_dir, "04.MR_PRESSO_Results.csv")
  write.csv(presso_all, presso_file, row.names = FALSE)
  cat("[", timestamp(), "] 04.MR_PRESSO_Results.csv saved\n")

  # 05. Steiger方向性分析汇总
  steiger_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = integrated_results$outcome,
    Steiger_correct_direction = integrated_results$steiger_correct_direction,
    Steiger_pval = integrated_results$steiger_pval,
    Directionality_significant = integrated_results$directionality_significant,
    stringsAsFactors = FALSE
  )
  steiger_file <- file.path(output_dir, "05.Steiger_Direction_Results.csv")
  write.csv(steiger_all, steiger_file, row.names = FALSE)
  cat("[", timestamp(), "] 05.Steiger_Direction_Results.csv saved\n")

  # 06. 三重筛选结果汇总
  three_filter_all <- data.frame(
    Exposure = integrated_results$gene_symbol,
    Outcome = integrated_results$outcome,
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
  tf_file <- file.path(output_dir, "06.Three_Filter_Results.csv")
  write.csv(three_filter_all, tf_file, row.names = FALSE)
  cat("[", timestamp(), "] 06.Three_Filter_Results.csv saved\n")

  # 07. 因果靶点筛选结果（三重筛选通过）
  causal_targets <- integrated_results[integrated_results$three_filter_pass == "YES", ]
  if (nrow(causal_targets) > 0) {
    causal_file <- file.path(output_dir, "07.MR_res_gene.csv")
    write.csv(causal_targets[, c("gene_symbol", "outcome", "n_snps", "ivw_beta", "ivw_se", "ivw_pval", "ivw_or")],
              causal_file, row.names = FALSE)
    cat("[", timestamp(), "] 07.MR_res_gene.csv saved (", nrow(causal_targets), " causal targets)\n")
  }

  # 08. ID映射与展示名（与eQTL/pQTL输出格式对齐）
  id_map <- data.frame(
    raw_id = as.character(integrated_results$gene_symbol),
    display_gene = resolve_display_ids(integrated_results$gene_symbol),
    stringsAsFactors = FALSE
  )
  id_map <- unique(id_map)
  file_08 <- file.path(output_dir, "08.Traditional_ID_Mapping.csv")
  write.csv(id_map, file_08, row.names = FALSE)
  cat("[", timestamp(), "] 08.Traditional_ID_Mapping.csv saved\n")

  # 09. 成功分析基因列表（展示名）
  success_genes <- unique(as.character(id_map$display_gene))
  success_genes <- success_genes[!is.na(success_genes) & nzchar(success_genes)]
  file_09 <- file.path(output_dir, "09.Traditional_Successful_Genes.csv")
  write.csv(data.frame(gene = success_genes, stringsAsFactors = FALSE), file_09, row.names = FALSE)
  cat("[", timestamp(), "] 09.Traditional_Successful_Genes.csv saved (", length(success_genes), " genes)\n")
}

cat("================================================================================\n")
cat("  Traditional MR Analysis Complete!\n")
cat("  Finished:", timestamp(), "\n")
cat("  Success:", sum(grepl("Success", analysis_summary$status)), "/", nrow(analysis_summary), "\n")
cat("  Results saved to:", file.path(output_dir, "tables"), "\n")
cat("================================================================================\n")

# Save session info + qs2 snapshot
dataset_name <- basename(normalizePath(output_dir, mustWork = FALSE))
qs2_path <- file.path(output_dir, sprintf("%s_traditional_node.qs2", dataset_name))
save_qs2_snapshot(list(
  integrated_results = integrated_results,
  analysis_summary = analysis_summary,
  exposures_run = exposures_to_run
), qs2_path)
save_session_info_snapshot(file.path(output_dir, "logs"), dataset_name, run_timestamp)

# Auto-generate Word report
cat("\n================================================================================\n")
cat("  STEP: AUTO-GENERATING MR REPORT\n")
cat("================================================================================\n")
cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- grep("^--file=", cmd_args, value = TRUE)
if (length(script_path) > 0) {
  script_dir <- dirname(gsub("^--file=", "", script_path[1]))
} else {
  script_dir <- "/media/desk16/share/secure/MR/scripts"
}
report_script <- file.path(script_dir, "generate_mr_report.R")
if (file.exists(report_script)) {
  report_output <- file.path(output_dir, paste0("MR_Report_", run_timestamp, ".docx"))
  report_cmd <- paste("Rscript", report_script,
                      "-i", output_dir,
                      "-o", report_output,
                      "-t", "traditional",
                      "-T", run_timestamp)
  cat("Generating report...\n")
  system(report_cmd)
  cat("✓ Report saved to:", report_output, "\n")
} else {
  cat("Warning: generate_mr_report.R not found, skipping auto-report generation\n")
}
