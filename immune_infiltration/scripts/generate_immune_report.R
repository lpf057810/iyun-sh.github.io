#!/usr/bin/env Rscript

# =============================================================================
# Immune Infiltration Analysis Report Generator
# 版式要求：
# 1) 中文宋体、英文 Times New Roman
# 2) 正文首行缩进（用两个中文全角空格）
# 3) 图表标题居中，分析正文左对齐不居中
# 4) 表格内容使用英文
# 5) 图片分析文本由统计表数值和文本拼接生成
# =============================================================================

library(optparse)
library(officer)
library(flextable)
library(magrittr)

# -----------------------------
# 参数解析
# -----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL),
  make_option(c("-o", "--output"), type = "character", default = NULL),
  make_option(c("-c", "--config"), type = "character", default = NULL),
  make_option(c("-t", "--timestamp"), type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$timestamp)) opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (is.null(opt$output)) opt$output <- paste0("Immune_Report_", opt$timestamp, ".docx")
if (is.null(opt$input) || !dir.exists(opt$input)) {
  stop("Input result directory is missing or does not exist.")
}
result_dir <- opt$input

cat("报告生成中...\n")

# -----------------------------
# 全局映射与辅助函数
# -----------------------------
method_names <- c(
  "cibersort" = "CIBERSORT",
  "epic" = "EPIC",
  "estimate" = "ESTIMATE",
  "ips" = "IPS",
  "mcpcounter" = "MCPcounter",
  "ssgsea" = "ssGSEA",
  "timer" = "TIMER",
  "xcell" = "xCell"
)

method_docs <- list(
  cibersort = "CIBERSORT基于线性支持向量机（SVM）算法，利用已知免疫细胞特征矩阵反卷积bulk表达数据，估算免疫细胞相对丰度。",
  epic = "EPIC基于约束最小二乘回归模型，估算肿瘤微环境中免疫细胞和肿瘤细胞的绝对丰度。",
  estimate = "ESTIMATE通过基因表达特征计算Stromal Score、Immune Score和ESTIMATE Score。",
  ips = "IPS综合免疫细胞浸润、抗原呈递和免疫检查点等特征，计算免疫表型评分。",
  mcpcounter = "MCPcounter基于细胞特异性标记基因估算多类免疫细胞及基质细胞丰度。",
  ssgsea = "ssGSEA通过单样本基因集富集分数评估各免疫细胞相关基因集活性。",
  timer = "TIMER综合多种算法评估多类核心免疫细胞在样本中的浸润水平。",
  xcell = "xCell通过细胞类型特异性基因集推断细胞富集评分。"
)

get_r_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("Not installed")
  as.character(utils::packageVersion(pkg))
}

get_python_version <- function() {
  out <- tryCatch(system2("python3", "--version", stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  if (!length(out)) return("Not found")
  trimws(gsub("^Python\\s+", "", out[1]))
}

get_python_pkg_version <- function(pkg) {
  out <- tryCatch(system2("python3", c("-m", "pip", "show", pkg), stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  if (!length(out)) return("Not found")
  vline <- grep("^Version:", out, value = TRUE)
  if (!length(vline)) return("Not installed")
  trimws(sub("^Version:\\s*", "", vline[1]))
}

build_software_catalog <- function() {
  r_meta <- data.frame(
    Category = rep("R package", 10),
    Name = c("optparse", "IOBR", "GSVA", "data.table", "tidyverse", "ggplot2", "rstatix", "officer", "flextable", "qs"),
    Version = c(
      get_r_pkg_version("optparse"),
      get_r_pkg_version("IOBR"),
      get_r_pkg_version("GSVA"),
      get_r_pkg_version("data.table"),
      get_r_pkg_version("tidyverse"),
      get_r_pkg_version("ggplot2"),
      get_r_pkg_version("rstatix"),
      get_r_pkg_version("officer"),
      get_r_pkg_version("flextable"),
      get_r_pkg_version("qs")
    ),
    Official_URL = c(
      "https://cran.r-project.org/package=optparse",
      "https://github.com/IOBR/IOBR",
      "https://bioconductor.org/packages/GSVA",
      "https://cran.r-project.org/package=data.table",
      "https://www.tidyverse.org/",
      "https://cran.r-project.org/package=ggplot2",
      "https://cran.r-project.org/package=rstatix",
      "https://cran.r-project.org/package=officer",
      "https://cran.r-project.org/package=flextable",
      "https://cran.r-project.org/package=qs"
    ),
    stringsAsFactors = FALSE
  )

  py_meta <- data.frame(
    Category = rep("Python package", 6),
    Name = c("pandas", "numpy", "matplotlib", "seaborn", "scipy", "python-docx"),
    Version = c(
      get_python_pkg_version("pandas"),
      get_python_pkg_version("numpy"),
      get_python_pkg_version("matplotlib"),
      get_python_pkg_version("seaborn"),
      get_python_pkg_version("scipy"),
      get_python_pkg_version("python-docx")
    ),
    Official_URL = c(
      "https://pandas.pydata.org/",
      "https://numpy.org/",
      "https://matplotlib.org/",
      "https://seaborn.pydata.org/",
      "https://scipy.org/",
      "https://python-docx.readthedocs.io/"
    ),
    stringsAsFactors = FALSE
  )

  env_meta <- data.frame(
    Category = c("Software", "Software"),
    Name = c("R", "Python"),
    Version = c(as.character(getRversion()), get_python_version()),
    Official_URL = c("https://www.r-project.org/", "https://www.python.org/"),
    stringsAsFactors = FALSE
  )

  rbind(env_meta, r_meta, py_meta)
}

software_info <- build_software_catalog()

find_methods <- function(dir) {
  dirs <- list.dirs(dir, full.names = FALSE, recursive = FALSE)
  dirs[grepl("_output$", dirs)]
}

get_method_key <- function(x) gsub("_output$", "", x)

load_res <- function(dir, name) {
  res_f <- list.files(dir, "^01\\..*_res\\.csv$", full.names = TRUE)
  stat_f <- list.files(dir, "^02\\.stat\\..*\\.csv$", full.names = TRUE)
  list(
    res = if (length(res_f)) tryCatch(read.csv(res_f[1], check.names = FALSE), error = function(e) NULL) else NULL,
    stat = if (length(stat_f)) tryCatch(read.csv(stat_f[1], check.names = FALSE), error = function(e) NULL) else NULL,
    name = name
  )
}

get_config_value <- function(lines, key) {
  hit <- grep(paste0("^\\s*", key, "\\s*="), lines, value = TRUE)
  if (!length(hit)) return(NULL)
  trimws(sub("^[^=]*=", "", hit[1]))
}

fmt_num <- function(x, digits = 4) {
  if (is.null(x) || is.na(x)) return("NA")
  if (abs(x) < 1e-4 && x != 0) return(format(x, scientific = TRUE, digits = 3))
  format(round(x, digits), nsmall = min(digits, 4), scientific = FALSE, trim = TRUE)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 3))
  format(round(p, 4), nsmall = 4, trim = TRUE)
}

clean_cell <- function(x) {
  if (is.null(x) || is.na(x)) return("Unknown")
  gsub("_", " ", as.character(x), fixed = TRUE)
}

# -----------------------------
# 字体与段落风格（中文宋体 + 英文Times New Roman）
# -----------------------------
make_text_prop <- function(size = 12, bold = FALSE) {
  fp_text(
    font.size = size,
    bold = bold,
    font.family = "Times New Roman",
    hansi.family = "Times New Roman",
    cs.family = "Times New Roman",
    eastasia.family = "SimSun"
  )
}

add_text <- function(doc, text, align = "left", indent = FALSE, bold = FALSE, size = 12) {
  prefix <- if (indent && nzchar(text)) "　　" else ""
  p <- fpar(
    ftext(paste0(prefix, text), prop = make_text_prop(size = size, bold = bold)),
    fp_p = fp_par(text.align = align)
  )
  body_add_fpar(doc, p)
}

add_blank <- function(doc) {
  body_add_par(doc, "", style = "Normal")
}

style_table <- function(ft, first_col_left = FALSE) {
  ft %>%
    flextable::bold(part = "header") %>%
    flextable::fontsize(size = 10, part = "all") %>%
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::align(align = "center", part = "body") %>%
    { if (first_col_left) flextable::align(., j = 1, align = "left", part = "body") else . } %>%
    flextable::border_remove() %>%
    flextable::border(part = "all", border.top = fp_border(width = 1.5, color = "black")) %>%
    flextable::border(part = "header", border.bottom = fp_border(width = 0.5, color = "black")) %>%
    flextable::border(part = "body", border.bottom = fp_border(width = 1.5, color = "black"), i = nrow(ft$body$dataset)) %>%
    flextable::set_table_properties(layout = "autofit")
}

compute_stat_enriched <- function(stat, res, group_info) {
  if (is.null(stat) || !nrow(stat)) return(stat)

  stat$mean_group1 <- NA_real_
  stat$mean_group2 <- NA_real_
  stat$delta_g2_minus_g1 <- NA_real_
  stat$direction <- "undetermined"

  needed <- c("ImmuneCell", "group1", "group2")
  if (!all(needed %in% colnames(stat))) return(stat)
  if (is.null(res) || is.null(group_info)) return(stat)
  if (!all(c("ID", "sample", "group") %in% c(colnames(res), colnames(group_info)))) return(stat)

  for (i in seq_len(nrow(stat))) {
    cell <- as.character(stat$ImmuneCell[i])
    g1 <- as.character(stat$group1[i])
    g2 <- as.character(stat$group2[i])
    if (!(cell %in% colnames(res))) next

    s1 <- group_info$sample[group_info$group == g1]
    s2 <- group_info$sample[group_info$group == g2]
    v1 <- suppressWarnings(as.numeric(res[res$ID %in% s1, cell]))
    v2 <- suppressWarnings(as.numeric(res[res$ID %in% s2, cell]))
    m1 <- mean(v1, na.rm = TRUE)
    m2 <- mean(v2, na.rm = TRUE)

    if (is.finite(m1) && is.finite(m2)) {
      stat$mean_group1[i] <- m1
      stat$mean_group2[i] <- m2
      stat$delta_g2_minus_g1[i] <- m2 - m1
      stat$direction[i] <- ifelse(
        m2 > m1,
        paste0(g2, " > ", g1),
        ifelse(m2 < m1, paste0(g2, " < ", g1), paste0(g2, " = ", g1))
      )
    }
  }
  stat
}

build_data_driven_summary <- function(stat_enriched, method_name) {
  if (is.null(stat_enriched) || !nrow(stat_enriched) || !("p" %in% colnames(stat_enriched))) {
    return(paste0(method_name, "统计结果文件不完整，无法进行差异统计解读。"))
  }

  pvals <- suppressWarnings(as.numeric(stat_enriched$p))
  sig <- stat_enriched[!is.na(pvals) & pvals < 0.05, , drop = FALSE]
  if (!nrow(sig)) {
    return(paste0("基于统计表，", method_name, "在两组间未检出显著差异条目（p < 0.05）。"))
  }

  ord_col <- if ("p.adj" %in% colnames(sig)) suppressWarnings(as.numeric(sig[["p.adj"]])) else pvals[!is.na(pvals) & pvals < 0.05]
  ord <- order(ord_col, na.last = TRUE)
  sig <- sig[ord, , drop = FALSE]
  top_n <- min(3, nrow(sig))
  top_sig <- sig[seq_len(top_n), , drop = FALSE]

  details <- character(0)
  for (i in seq_len(nrow(top_sig))) {
    row <- top_sig[i, , drop = FALSE]
    cell_txt <- clean_cell(row$ImmuneCell)
    g1 <- as.character(row$group1)
    g2 <- as.character(row$group2)
    p_txt <- fmt_p(suppressWarnings(as.numeric(row$p)))
    m1_txt <- fmt_num(suppressWarnings(as.numeric(row$mean_group1)))
    m2_txt <- fmt_num(suppressWarnings(as.numeric(row$mean_group2)))
    dir_txt <- as.character(row$direction)
    details <- c(details, paste0(cell_txt, " (", dir_txt, ", ", g1, "=", m1_txt, ", ", g2, "=", m2_txt, ", p=", p_txt, ")"))
  }

  paste0(
    "基于统计表（p < 0.05），", method_name, "共识别到", nrow(sig),
    "个显著差异条目。主要差异包括：", paste(details, collapse = "; "), "。"
  )
}

# -----------------------------
# 读取结果
# -----------------------------
available <- find_methods(result_dir)
cat("发现方法:", paste(available, collapse = ", "), "\n")

all_res <- list()
for (m in available) {
  d <- file.path(result_dir, m)
  k <- get_method_key(m)
  method_name <- method_names[k]
  if (is.na(method_name)) method_name <- m
  all_res[[m]] <- load_res(d, method_name)
}

# -----------------------------
# 读取配置与分组信息
# -----------------------------
config_params <- list()
group_info <- NULL

if (!is.null(opt$config) && file.exists(opt$config)) {
  cfg_lines <- readLines(opt$config, warn = FALSE)
  config_params$methods <- get_config_value(cfg_lines, "methods")
  config_params$expression <- get_config_value(cfg_lines, "expression")
  config_params$group <- get_config_value(cfg_lines, "group")
  config_params$signature <- get_config_value(cfg_lines, "signature")
  config_params$base_dir <- get_config_value(cfg_lines, "base_dir")
  config_params$create_subdirs <- get_config_value(cfg_lines, "create_subdirs")
  config_params$parallel_enabled <- get_config_value(cfg_lines, "parallel_enabled")
  config_params$max_workers <- get_config_value(cfg_lines, "max_workers")
  config_params$cache_enabled <- get_config_value(cfg_lines, "enabled")
  config_params$cache_dir <- get_config_value(cfg_lines, "cache_dir")
  config_params$force <- get_config_value(cfg_lines, "force")
  config_params$perm <- get_config_value(cfg_lines, "perm")
  config_params$tissue <- get_config_value(cfg_lines, "tissue")

  group_path <- get_config_value(cfg_lines, "group")
  if (!is.null(group_path) && nzchar(group_path)) {
    if (!grepl("^/", group_path)) {
      group_path <- file.path(dirname(opt$config), group_path)
    }
    if (file.exists(group_path)) {
      group_info <- tryCatch(read.csv(group_path, check.names = FALSE), error = function(e) NULL)
    }
  }
}

groups <- if (!is.null(group_info) && "group" %in% colnames(group_info)) unique(as.character(group_info$group)) else NULL

# 与 Python 主程序一致：cache_dir 为空时使用 <base_dir>/qs2
effective_cache_dir <- if (!is.null(config_params$cache_dir) && nzchar(config_params$cache_dir)) {
  config_params$cache_dir
} else if (!is.null(config_params$base_dir) && nzchar(config_params$base_dir)) {
  file.path(config_params$base_dir, "qs2")
} else {
  NA_character_
}

# -----------------------------
# 文档生成
# -----------------------------
doc <- read_docx()

# 标题
doc <- add_text(doc, "免疫浸润分析报告", align = "center", bold = TRUE, size = 16)
doc <- add_text(doc, paste("生成时间:", opt$timestamp), align = "center", size = 11)
doc <- add_blank(doc)

# 1. 方法
doc <- add_text(doc, "1. 方法", bold = TRUE, size = 14)
doc <- add_blank(doc)

doc <- add_text(doc, "1.1 软件、工具包与版本", bold = TRUE, size = 12)
doc <- add_text(
  doc,
  "本研究在Linux环境下运行免疫浸润分析流程，采用R与Python协同实现数据预处理、免疫浸润估计、统计检验、可视化及报告生成。表1系统列出本研究使用的软件与工具包名称、版本号及官方链接。",
  indent = TRUE,
  size = 12
)
doc <- add_blank(doc)

ft_sw <- flextable(software_info)
ft_sw <- style_table(ft_sw, first_col_left = TRUE)
doc <- body_add_flextable(doc, ft_sw)
doc <- add_blank(doc)

doc <- add_text(doc, "1.2 分析流程", bold = TRUE, size = 12)
flow_text <- paste0(
  "分析流程包括：",
  "（1）读取表达矩阵与分组文件并进行格式校验；",
  "（2）按配置调用免疫浸润算法（CIBERSORT、EPIC、ESTIMATE、IPS、MCPcounter、ssGSEA、TIMER、xCell）；",
  "（3）输出各方法浸润结果、组间比较统计表与图像（堆叠图、箱线图、差异图）；",
  "（4）基于统计表自动拼接结果描述并生成Word报告。"
)
doc <- add_text(doc, flow_text, indent = TRUE, size = 12)
doc <- add_blank(doc)

doc <- add_text(doc, "1.3 免疫浸润算法说明", bold = TRUE, size = 12)
doc <- add_blank(doc)

meth_count <- 0
for (m in available) {
  k <- get_method_key(m)
  info <- method_docs[[k]]
  method_name <- method_names[k]
  if (is.na(method_name)) method_name <- m
  if (!is.null(info)) {
    meth_count <- meth_count + 1
    doc <- add_text(doc, paste0("1.3.", meth_count, " ", method_name), bold = TRUE, size = 12)
    doc <- add_text(doc, info, indent = TRUE, size = 12)

    if (k == "cibersort" && !is.null(config_params$perm) && nzchar(config_params$perm)) {
      doc <- add_text(doc, paste0("参数设置：perm = ", config_params$perm), indent = TRUE, size = 12)
    }
    if (k == "timer" && !is.null(config_params$tissue) && nzchar(config_params$tissue)) {
      doc <- add_text(doc, paste0("参数设置：tissue = ", config_params$tissue), indent = TRUE, size = 12)
    }
    doc <- add_blank(doc)
  }
}

doc <- add_text(doc, "1.4 参数设置", bold = TRUE, size = 12)
param_text <- paste0(
  "本次分析参数为：methods = ", ifelse(is.null(config_params$methods) || !nzchar(config_params$methods), "NA", config_params$methods),
  "；expression = ", ifelse(is.null(config_params$expression) || !nzchar(config_params$expression), "NA", config_params$expression),
  "；group = ", ifelse(is.null(config_params$group) || !nzchar(config_params$group), "NA", config_params$group),
  "；signature = ", ifelse(is.null(config_params$signature) || !nzchar(config_params$signature), "NA", config_params$signature),
  "；base_dir = ", ifelse(is.null(config_params$base_dir) || !nzchar(config_params$base_dir), "NA", config_params$base_dir),
  "；create_subdirs = ", ifelse(is.null(config_params$create_subdirs) || !nzchar(config_params$create_subdirs), "NA", config_params$create_subdirs),
  "；parallel_enabled = ", ifelse(is.null(config_params$parallel_enabled) || !nzchar(config_params$parallel_enabled), "NA", config_params$parallel_enabled),
  "；max_workers = ", ifelse(is.null(config_params$max_workers) || !nzchar(config_params$max_workers), "NA", config_params$max_workers),
  "；cache_enabled = ", ifelse(is.null(config_params$cache_enabled) || !nzchar(config_params$cache_enabled), "NA", config_params$cache_enabled),
  "；cache_dir = ", ifelse(is.na(effective_cache_dir), "NA", effective_cache_dir),
  "；force = ", ifelse(is.null(config_params$force) || !nzchar(config_params$force), "NA", config_params$force),
  "；CIBERSORT perm = ", ifelse(is.null(config_params$perm) || !nzchar(config_params$perm), "NA", config_params$perm),
  "；TIMER tissue = ", ifelse(is.null(config_params$tissue) || !nzchar(config_params$tissue), "NA", config_params$tissue), "。"
)
doc <- add_text(doc, param_text, indent = TRUE, size = 12)
doc <- add_blank(doc)

doc <- add_text(doc, "1.5 统计分析", bold = TRUE, size = 12)
doc <- add_text(
  doc,
  "组间比较采用Mann-Whitney U检验，多重比较校正采用Benjamini-Hochberg（BH）方法，显著性阈值设定为p < 0.05。",
  indent = TRUE,
  size = 12
)
doc <- add_blank(doc)

# 2. 结果
doc <- add_text(doc, "2. 结果", bold = TRUE, size = 14)
doc <- add_blank(doc)

n_total <- if (!is.null(group_info) && "sample" %in% colnames(group_info)) nrow(group_info) else 0
if (n_total > 0) {
  sample_summary <- data.frame(
    Metric = c("Total samples", if (!is.null(groups)) paste0(groups, " samples") else character(0)),
    Value = c(as.character(n_total), if (!is.null(groups)) sapply(groups, function(g) as.character(sum(group_info$group == g))) else character(0)),
    stringsAsFactors = FALSE
  )
  doc <- add_text(doc, "2.1 Sample Summary", bold = TRUE, size = 12)
  doc <- add_blank(doc)
  ft_samp <- style_table(flextable(sample_summary), first_col_left = TRUE)
  doc <- body_add_flextable(doc, ft_samp)
  doc <- add_blank(doc)
}

method_summary <- data.frame(
  No = seq_along(available),
  Method = sapply(available, function(m) {
    k <- get_method_key(m)
    n <- method_names[k]
    if (is.na(n)) n <- m
    n
  }),
  stringsAsFactors = FALSE
)

doc <- add_text(doc, "2.2 Methods Included", bold = TRUE, size = 12)
doc <- add_blank(doc)
ft_meth <- style_table(flextable(method_summary), first_col_left = FALSE)
doc <- body_add_flextable(doc, ft_meth)
doc <- add_blank(doc)

# 各方法结果
fig_idx <- 0
table_idx <- 0
auto_desc <- list()

for (m in available) {
  k <- get_method_key(m)
  method_name <- method_names[k]
  if (is.na(method_name)) method_name <- m
  res_obj <- all_res[[m]]
  if (is.null(res_obj)) next

  stat_enriched <- compute_stat_enriched(res_obj$stat, res_obj$res, group_info)
  method_summary_text <- build_data_driven_summary(stat_enriched, method_name)
  auto_desc[[method_name]] <- method_summary_text

  sec_title <- paste0("2.", length(auto_desc) + 2, " ", method_name, "分析结果")
  doc <- add_text(doc, sec_title, bold = TRUE, size = 12)
  doc <- add_blank(doc)

  # Figure 1: stacked bar
  sb <- list.files(file.path(result_dir, m), "^01\\..*stacked_bar\\.png$", full.names = TRUE)
  if (length(sb)) {
    fig_idx <- fig_idx + 1
    doc <- add_text(doc, paste0("Figure ", fig_idx, ". ", method_name, " stacked composition plot"), align = "center", bold = TRUE, size = 11)
    doc <- body_add_img(doc, src = sb[1], width = 6, height = 4)
    doc <- add_text(doc, paste0("图像解读：", method_summary_text), indent = TRUE, size = 12)
    doc <- add_blank(doc)
  }

  # Figure 2: boxplot
  bx <- list.files(file.path(result_dir, m), "^02\\..*boxplot\\.png$", full.names = TRUE)
  if (length(bx)) {
    fig_idx <- fig_idx + 1
    doc <- add_text(doc, paste0("Figure ", fig_idx, ". ", method_name, " group comparison boxplot"), align = "center", bold = TRUE, size = 11)
    doc <- body_add_img(doc, src = bx[1], width = 6, height = 4)
    doc <- add_text(doc, paste0("图像解读：该图展示组间分布差异。", method_summary_text), indent = TRUE, size = 12)
    doc <- add_blank(doc)
  }

  # Figure 3: DE plot
  de <- list.files(file.path(result_dir, m), "^03\\.DE\\..*\\.png$", full.names = TRUE)
  if (length(de)) {
    fig_idx <- fig_idx + 1
    doc <- add_text(doc, paste0("Figure ", fig_idx, ". ", method_name, " differential analysis plot"), align = "center", bold = TRUE, size = 11)
    doc <- body_add_img(doc, src = de[1], width = 6, height = 4)
    doc <- add_text(doc, paste0("图像解读：该图反映了显著差异条目的变化方向与显著性水平。", method_summary_text), indent = TRUE, size = 12)
    doc <- add_blank(doc)
  }

  # 显著性统计表（英文）
  if (!is.null(stat_enriched) && nrow(stat_enriched) > 0 && "p" %in% colnames(stat_enriched)) {
    pvals <- suppressWarnings(as.numeric(stat_enriched$p))
    sig_dt <- stat_enriched[!is.na(pvals) & pvals < 0.05, , drop = FALSE]
    if (nrow(sig_dt) > 0) {
      keep_cols <- c("ImmuneCell", ".y.", "group1", "group2", "n1", "n2", "mean_group1", "mean_group2", "delta_g2_minus_g1", "statistic", "p", "p.adj", "p.signif", "direction")
      keep_cols <- keep_cols[keep_cols %in% colnames(sig_dt)]
      sig_dt <- sig_dt[, keep_cols, drop = FALSE]

      col_map <- c(
        "ImmuneCell" = "Feature",
        ".y." = "ScoreType",
        "group1" = "Group1",
        "group2" = "Group2",
        "n1" = "N1",
        "n2" = "N2",
        "mean_group1" = "Mean_Group1",
        "mean_group2" = "Mean_Group2",
        "delta_g2_minus_g1" = "Delta_G2_minus_G1",
        "statistic" = "Statistic",
        "p" = "P",
        "p.adj" = "P_adj",
        "p.signif" = "Signif",
        "direction" = "Direction"
      )
      colnames(sig_dt) <- unname(col_map[colnames(sig_dt)])

      table_idx <- table_idx + 1
      doc <- add_text(doc, paste0("Table ", table_idx, ". Significant features (p < 0.05) - ", method_name), align = "center", bold = TRUE, size = 11)
      doc <- add_blank(doc)
      ft <- style_table(flextable(sig_dt), first_col_left = TRUE)
      doc <- body_add_flextable(doc, ft)
      doc <- add_blank(doc)
    }
  }
}

# 统一的结果文字汇总
doc <- add_text(doc, "2. Result Description", bold = TRUE, size = 12)
doc <- add_blank(doc)
for (nm in names(auto_desc)) {
  doc <- add_text(doc, paste0(nm, "：", auto_desc[[nm]]), indent = TRUE, size = 12)
}
doc <- add_blank(doc)

# 保存
print(doc, target = opt$output)
cat("报告已生成:", opt$output, "\n")
