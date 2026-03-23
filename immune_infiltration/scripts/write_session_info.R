#!/usr/bin/env Rscript
# 将 sessionInfo() 写入文本文件，供复现与论文方法部分引用

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Output text file path [required]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$output) || !nzchar(opt$output)) {
  stop("Must specify -o / --output path")
}

out_dir <- dirname(opt$output)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

sink(opt$output, split = FALSE)
on.exit(sink(), add = TRUE)

cat("# R sessionInfo\n")
cat("# Generated: ", format(Sys.time(), usetz = TRUE), "\n\n", sep = "")
print(sessionInfo())
cat("\n")
