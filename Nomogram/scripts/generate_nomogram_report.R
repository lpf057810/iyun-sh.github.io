#!/usr/bin/env Rscript
# 列线图 Word 报告生成器
# - 首行缩进（Python XML 后处理）
# - 混合字体：中文 SimSun，英文 Times New Roman
# - 图片居中嵌入，图注居中
library(optparse)
library(officer)
library(flextable)
library(magrittr)

# ─────────────────────────────────────────
# 命令行参数
# ─────────────────────────────────────────
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "结果目录（默认 results/Nomogram）"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "输出 docx 路径"),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL,
              help = "运行时间戳")
)
opt <- parse_args(OptionParser(option_list = option_list))

timestamp <- if (!is.null(opt$timestamp)) opt$timestamp else format(Sys.time(), "%Y%m%d_%H%M%S")
result_dir <- if (!is.null(opt$input)) opt$input else "results/Nomogram"
output_docx <- if (!is.null(opt$output)) opt$output else paste0("nomogram_report_", timestamp, ".docx")

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- args[grep(file_arg, args)]
  if (length(hit) == 0) return(getwd())
  dirname(normalizePath(sub(file_arg, "", hit[[1]]), mustWork = FALSE))
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
report_template <- file.path(project_root, "templates", "nomogram_report_template.docx")

# ─────────────────────────────────────────
# 加载分析结果
# ─────────────────────────────────────────
params_file <- file.path(result_dir, "_report_params.RData")
if (file.exists(params_file)) {
  load(params_file)
} else {
  report_params <- list(C_INDEX = NA, R2 = NA, HL_PVAL = NA, AUC_VAL = NA,
                        CALIB_MAE = NA, CALIBRATE_B = 1000, LRM_MAXIT = 1000,
                        DATASET_NAME = "Nomogram", GENES = "Gene1",
                        COMPARISON = c("Control", "Case"))
  genes <- "Gene1"
}

if (!is.null(genes) && length(genes) == 1 && grepl(",", genes[1])) {
  genes <- strsplit(genes[1], ",\\s*")[[1]]
}
dataset_name <- if (!is.null(report_params$DATASET_NAME)) report_params$DATASET_NAME else basename(result_dir)
n_genes <- if (!is.null(genes)) length(genes) else 1
comparison <- report_params$COMPARISON

# 读取模型系数
model_coef <- NULL
coef_path <- file.path(result_dir, "05.model_coefficients.csv")
if (file.exists(coef_path)) {
  model_coef <- read.csv(coef_path, stringsAsFactors = FALSE)
}
nomogram_stats <- NULL
nom_stats_path <- file.path(result_dir, "01.nomogram_stats.csv")
if (file.exists(nom_stats_path)) {
  nomogram_stats <- read.csv(nom_stats_path, stringsAsFactors = FALSE)
}

# 读取 ROC 结果
roc_data <- NULL
roc_path <- file.path(result_dir, "04.roc_results.csv")
if (file.exists(roc_path)) {
  df <- read.csv(roc_path, stringsAsFactors = FALSE)
  roc_data <- setNames(as.character(df$Value), df$Metric)
}
calib_stats <- NULL
calib_stats_path <- file.path(result_dir, "02.calibrate_stats.csv")
if (file.exists(calib_stats_path)) {
  calib_stats <- read.csv(calib_stats_path, stringsAsFactors = FALSE)
}
dca_stats <- NULL
dca_stats_path <- file.path(result_dir, "03.dca_stats.csv")
if (file.exists(dca_stats_path)) {
  dca_stats <- read.csv(dca_stats_path, stringsAsFactors = FALSE)
}

auc_val  <- suppressWarnings(as.numeric(report_params$AUC_VAL))
c_index  <- suppressWarnings(as.numeric(report_params$C_INDEX))
hl_pval  <- report_params$HL_PVAL
hl_num   <- suppressWarnings(as.numeric(hl_pval))
calib_mae <- if (!is.null(report_params$CALIB_MAE) && !is.na(report_params$CALIB_MAE)) {
  suppressWarnings(as.numeric(report_params$CALIB_MAE))
} else NA
if (!is.null(calib_stats) && nrow(calib_stats) > 0 && all(c("metric", "value") %in% colnames(calib_stats))) {
  get_metric <- function(x) {
    hit <- calib_stats$value[calib_stats$metric == x]
    if (length(hit) == 0) return(NA)
    suppressWarnings(as.numeric(hit[1]))
  }
  mae_file <- get_metric("calibration_mae")
  hl_file <- get_metric("hosmer_lemeshow_p")
  if (!is.na(mae_file)) calib_mae <- mae_file
  if (!is.na(hl_file)) {
    hl_num <- hl_file
    hl_pval <- sprintf("%.4g", hl_file)
  }
}

# ─────────────────────────────────────────
# 字体配置
# ─────────────────────────────────────────
CN_FONT  <- "SimSun"          # 宋体（中文）
EN_FONT  <- "Times New Roman" # 新罗马（英文）
BODY_SIZE <- 12               # 正文 12pt
HEAD_SIZE <- 14               # 标题 14pt

# ─────────────────────────────────────────
# 辅助函数：构建混合字体段落
# ─────────────────────────────────────────
build_fpar <- function(text,
                       align      = "left",
                       font_size  = BODY_SIZE,
                       bold       = FALSE,
                       p_top      = 0,
                       p_bot      = 0,
                       first_line = 0,
                       word_style = NULL) {
  if (first_line > 0 && nzchar(trimws(text))) {
    text <- paste0("\u3000\u3000", text)
  }
  para_args <- list(
    text.align     = align,
    padding.top    = p_top,
    padding.bottom = p_bot
  )
  if (!is.null(word_style) && is.character(word_style) && nzchar(word_style)) {
    para_args$word_style <- word_style
  }
  para_props <- do.call(fp_par, para_args)

  chars   <- strsplit(text, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(fpar("", fp_p = para_props))

  is_cn   <- grepl("[一-龥]", chars)
  runs    <- list()
  start   <- 1L

  for (i in seq_along(chars)) {
    is_break <- i == length(chars) || is_cn[[i + 1L]] != is_cn[[i]]
    if (is_break) {
      seg <- paste(chars[start:i], collapse = "")
      ff  <- if (is_cn[[start]]) CN_FONT else EN_FONT
      runs[[length(runs) + 1L]] <- ftext(
        seg,
        prop = fp_text(font.family = ff, font.size = font_size, bold = bold)
      )
      start <- i + 1L
    }
  }

  do.call(fpar, c(runs, list(fp_p = para_props)))
}

# ─────────────────────────────────────────
# 辅助函数：空段落
# ─────────────────────────────────────────
mk_blank <- function(doc) {
  body_add_par(doc, "", style = "Normal")
}

# ─────────────────────────────────────────
# 辅助函数：图片 + 图注（与旧版一致）
# ─────────────────────────────────────────
mk_figure <- function(doc, img_path, width, height, figure_label,
                       caption_text, caption_note = NULL) {
  if (!is.null(img_path) && file.exists(img_path)) {
    ext <- tolower(tools::file_ext(img_path))
    img_block <- fpar(
      external_img(src = img_path, width = width, height = height),
      fp_p = fp_par(text.align = "center")
    )
    doc <- body_add_fpar(doc, value = img_block)
    doc <- mk_blank(doc)
  }

  if (!is.null(figure_label)) {
    doc <- body_add_fpar(doc,
      value = build_fpar(figure_label, align = "center", font_size = BODY_SIZE),
      style = "Normal"
    )
  }

  if (!is.null(caption_text)) {
    caption_line <- if (startsWith(caption_text, "图注：")) caption_text else paste0("图注：", caption_text)
    doc <- body_add_fpar(doc,
      value = build_fpar(caption_line, align = "center",
                         font_size = BODY_SIZE - 1),
      style = "Normal"
    )
  }

  doc <- mk_blank(doc)
  doc
}

# ─────────────────────────────────────────
# 构建文档
# ─────────────────────────────────────────
doc <- if (file.exists(report_template)) {
  cat("Using report template:", report_template, "\n")
  read_docx(path = report_template)
} else {
  read_docx()
}

# 标题
doc <- body_add_fpar(doc,
  value = build_fpar("列线图预测模型分析报告",
                     align = "center", font_size = HEAD_SIZE + 4, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 1. 方法 ──────────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("1. 方法", font_size = HEAD_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

pkg_rms  <- tryCatch(as.character(utils::packageVersion("rms")),        error = function(e) "unknown")
pkg_rmda <- tryCatch(as.character(utils::packageVersion("rmda")),        error = function(e) "unknown")
pkg_pROC <- tryCatch(as.character(utils::packageVersion("pROC")),        error = function(e) "unknown")
r_ver    <- as.character(getRversion())

method_text <- paste0(
  "采用R(v", r_ver, ", https://www.r-project.org)构建列线图；",
  "rms(v", pkg_rms, ")、pROC(v", pkg_pROC, ")、rmda(v", pkg_rmda, ")（https://cran.r-project.org）",
  "分别用于Logistic建模与Bootstrap校准(B=", report_params$CALIBRATE_B, ", maxit=", report_params$LRM_MAXIT, ")、AUC计算和DCA评估。"
)
doc <- body_add_fpar(doc,
  value = build_fpar(method_text, font_size = BODY_SIZE, first_line = 0.21),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 2. 结果 ──────────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2. 结果", font_size = HEAD_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 2.1 模型系数 ─────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2.1 模型系数", font_size = BODY_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

if (!is.null(model_coef) && nrow(model_coef) > 0) {
  has_or    <- "OR"         %in% colnames(model_coef)
  has_range <- "Range"      %in% colnames(model_coef)
  has_pts   <- "Points_pct" %in% colnames(model_coef)

  if (has_or) {
    model_coef$Coefficient <- round(as.numeric(model_coef$Coefficient), 4)
    model_coef$OR          <- round(as.numeric(model_coef$OR), 4)
    if (has_range && has_pts) {
      display <- data.frame(
        "Variable"          = model_coef$Variable,
        "Coefficient"      = model_coef$Coefficient,
        "OR"               = model_coef$OR,
        "95%CI"            = ifelse(is.na(model_coef$OR_lower), "-",
                          paste0("(", round(model_coef$OR_lower,4), "-", round(model_coef$OR_upper,4), ")")),
        "Expression Range"  = model_coef$Range,
        "Points Range"     = ifelse(is.na(model_coef$Points) | model_coef$Points == "-", "-", as.character(model_coef$Points)),
        "Contribution"      = ifelse(is.na(model_coef$Points_pct), "-", paste0(model_coef$Points_pct, "%")),
        stringsAsFactors   = FALSE,
        check.names        = FALSE
      )
    } else {
      display <- data.frame(
        "Variable"    = model_coef$Variable,
        "Coefficient" = model_coef$Coefficient,
        "OR"          = model_coef$OR,
        "95%CI"       = ifelse(is.na(model_coef$OR_lower), "-",
                      paste0("(", round(model_coef$OR_lower,4), "-", round(model_coef$OR_upper,4), ")")),
        stringsAsFactors = FALSE,
        check.names   = FALSE
      )
    }
  } else {
    model_coef$Coefficient <- round(as.numeric(model_coef$Coefficient), 4)
    display <- data.frame(
      "Variable"    = model_coef$Variable,
      "Coefficient" = model_coef$Coefficient,
      stringsAsFactors = FALSE
    )
  }

  ft <- flextable(display)
  ft <- ft %>% bold(part = "header")
  ft <- ft %>% fontsize(size = 10, part = "all")
  ft <- ft %>% font(fontname = EN_FONT, part = "all")
  ft <- ft %>% align(align = "center", part = "all")
  ft <- ft %>% border(part = "header", border.top    = officer::fp_border(width = 1.5))
  ft <- ft %>% border(part = "header", border.bottom = officer::fp_border(width = 0.5))
  ft <- ft %>% border(part = "body",   border.bottom = officer::fp_border(width = 1.5),
                       i = nrow(display))
  ft <- ft %>% set_table_properties(layout = "autofit", width = 1)
  doc <- body_add_fpar(doc,
    value = build_fpar("表1：Nomogram模型系数表", align = "center", font_size = BODY_SIZE),
    style = "Normal"
  )
  doc <- doc %>% body_add_flextable(ft, align = "center")
  doc <- mk_blank(doc)

  # 动态解读
  if (has_or) {
    gene_rows    <- model_coef[model_coef$Variable != "Intercept", ]
    risk_genes   <- gene_rows[as.numeric(gene_rows$OR) > 1, "Variable"]
    prot_genes   <- gene_rows[as.numeric(gene_rows$OR) < 1, "Variable"]

    coef_txt <- "表1结果显示，"
    if (length(risk_genes) > 0) {
      coef_txt <- paste0(coef_txt, paste(risk_genes, collapse = "、"),
                         "为危险因素（OR>1，系数为正），")
    }
    if (length(prot_genes) > 0) {
      coef_txt <- paste0(coef_txt, paste(prot_genes, collapse = "、"),
                         "为保护因素（OR<1，系数为负）。")
    }
    if (has_pts) {
      top_gene  <- gene_rows$Variable[which.max(as.numeric(gene_rows$Points_pct))]
      top_range <- gene_rows$Points[gene_rows$Variable == top_gene]
      top_pct   <- gene_rows$Points_pct[gene_rows$Variable == top_gene]
      coef_txt  <- paste0(coef_txt, sprintf(
        "其中%s对模型得分贡献最大（表达量范围对应的得分范围为%s，相对贡献为%s%%）。",
        top_gene, top_range, top_pct))
    }
    coef_txt <- paste0(coef_txt, "模型共纳入", n_genes, "个预测因子。相关结果表：`05.model_coefficients.csv`。")
    doc <- body_add_fpar(doc,
      value = build_fpar(coef_txt, font_size = BODY_SIZE, first_line = 0.21),
      style = "Normal"
    )
    doc <- mk_blank(doc)
  }
}

# ── 2.2 列线图 ───────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2.2 列线图模型", font_size = BODY_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

nom_img <- file.path(result_dir, "01.nomogram.png")
if (file.exists(nom_img)) {
  doc <- mk_figure(doc, nom_img, 6, 4,
    figure_label = "图1：列线图",
    caption_text = paste0(
      "图中最上方为Points刻度轴，中间为各预测因子取值轴，底部为Total Points与预测概率轴。",
      "各因子取值先映射为分值，分值累加得到总分，再映射为事件发生概率。",
      "模型共纳入", n_genes, "个预测因子（", paste(genes, collapse = "、"), "）。"
    )
  )
}
nom_txt <- paste0(
  "图1显示，列线图将各预测因子的表达水平映射为相应分数，各分数累加得到总分后映射为预测概率。",
  "列线图详细参数见表1，共纳入", n_genes, "个预测因子。",
  "相关结果表：`01.nomogram_stats.csv`。"
)
doc <- body_add_fpar(doc,
  value = build_fpar(nom_txt, font_size = BODY_SIZE, first_line = 0.21),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 2.3 ROC曲线 ──────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2.3 ROC曲线分析", font_size = BODY_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

roc_img <- file.path(result_dir, "04.roc.png")
auc_lvl <- if (is.na(auc_val)) "区分能力未知" else
  if (auc_val >= 0.9) "曲线紧贴左上角，区分能力极高" else
  if (auc_val >= 0.8) "曲线明显向左上角凸出，区分能力较高" else
  if (auc_val >= 0.7) "曲线向左上角偏移，区分能力中等" else
    "曲线接近对角线，区分能力较弱"

if (!is.null(roc_data) && !any(is.na(as.numeric(roc_data)))) {
  sens <- as.numeric(roc_data["Sensitivity"])
  spec <- as.numeric(roc_data["Specificity"])
  ppv  <- as.numeric(roc_data["PPV"])
  npv  <- as.numeric(roc_data["NPV"])
  thr  <- as.numeric(roc_data["Optimal Threshold"])
  roc_txt <- paste0(
    "图2显示，ROC曲线呈现\"", auc_lvl, "\"（AUC = ", round(auc_val, 4),
    "），一致性指数（C-index）= ", round(c_index, 4), "。",
    "在最佳阈值", round(thr, 4), "处，敏感度为", round(sens * 100, 1),
    "%，特异度为", round(spec * 100, 1), "%，",
    "阳性预测值（PPV）为", round(ppv * 100, 1),
    "%，阴性预测值（NPV）为", round(npv * 100, 1), "%。"
  )
} else {
  roc_txt <- paste0(
    "图2显示，ROC曲线呈现\"", auc_lvl, "\"（AUC = ",
    ifelse(is.na(auc_val), "NA", round(auc_val, 4)),
    "），一致性指数（C-index）= ",
    ifelse(is.na(c_index), "NA", round(c_index, 4)), "。"
  )
}

if (file.exists(roc_img)) {
  if (!is.null(roc_data)) {
    sens <- as.numeric(roc_data["Sensitivity"])
    spec <- as.numeric(roc_data["Specificity"])
  }
  doc <- mk_figure(doc, roc_img, 5, 4,
    figure_label = "图2：ROC曲线",
    caption_text = sprintf(
      "X轴为1-特异度（False Positive Rate），Y轴为敏感度（True Positive Rate）。AUC = %.4f；C-index = %.4f；敏感度 = %.1f%%；特异度 = %.1f%%。",
      auc_val, c_index, sens * 100, spec * 100
    )
  )
}
roc_full_txt <- roc_txt
if (!is.na(auc_val)) {
  roc_full_txt <- paste0(roc_full_txt, sprintf("客观统计：AUC=%.4f，C-index=%.4f。", auc_val, c_index))
}
roc_full_txt <- paste0(roc_full_txt, "相关结果表：`04.roc_results.csv`。")
doc <- body_add_fpar(doc,
  value = build_fpar(roc_full_txt, font_size = BODY_SIZE, first_line = 0.21),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 2.4 校准曲线 ─────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2.4 校准曲线分析", font_size = BODY_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

cal_img <- file.path(result_dir, "02.calibrate.png")

if (!is.na(hl_num)) {
  lvl <- if (hl_num > 0.05) "实线紧贴理想对角线，偏差较小" else
    if (hl_num > 0.01) "实线与对角线存在一定偏差，但整体趋势一致" else
      "实线与对角线偏离明显，校准欠佳"
  calib_txt <- paste0("图3显示，校准曲线呈现\"", lvl, "\"")
  if (!is.na(calib_mae)) calib_txt <- paste0(calib_txt, "，平均绝对误差（MAE）= ", round(calib_mae, 3))
  calib_txt <- paste0(calib_txt, "。Hosmer-Lemeshow检验结果为p = ", hl_pval,
    if (hl_num > 0.05) "，表明模型预测概率与实际观测概率具有良好的一致性。" else
      "，表明模型校准度存在一定偏差。"
  )
} else {
  calib_txt <- "图3显示，由于样本量限制，校准检验无法可靠执行。"
}

if (file.exists(cal_img)) {
  doc <- mk_figure(doc, cal_img, 5, 4,
    figure_label = "图3：校准曲线",
    caption_text = paste0(
      "X轴为模型预测概率，Y轴为实际观测概率；斜线为理想一致性线，实线为Bootstrap校准曲线。",
      "MAE = ", round(calib_mae, 3), "；H-L p = ", hl_pval, "。"
    )
  )
}
calib_full_txt <- paste0(
  calib_txt,
  "关键参数：校准评估采用Bootstrap重采样（B=", report_params$CALIBRATE_B, "）。",
  "相关结果表：`02.calibrate_stats.csv`。"
)
doc <- body_add_fpar(doc,
  value = build_fpar(calib_full_txt, font_size = BODY_SIZE, first_line = 0.21),
  style = "Normal"
)
doc <- mk_blank(doc)

# ── 2.5 决策曲线 ─────────────────────────
doc <- body_add_fpar(doc,
  value = build_fpar("2.5 决策曲线分析", font_size = BODY_SIZE, bold = TRUE),
  style = "Normal"
)
doc <- mk_blank(doc)

dca_img <- file.path(result_dir, "03.dca.png")
if (file.exists(dca_img)) {
  doc <- mk_figure(doc, dca_img, 6, 4,
    figure_label = "图4：决策曲线（DCA）",
    caption_text = "X轴为阈值概率（threshold probability），Y轴为净获益（net benefit）。红色曲线为模型策略，黑/灰线分别表示treat-all与treat-none策略。"
  )
}
dca_txt <- "图4显示，决策曲线分析显示，列线图模型在阈值区间内整体高于极端策略，提示具有潜在临床净获益。"
if (!is.null(dca_stats) && nrow(dca_stats) > 0 && all(c("metric", "value") %in% colnames(dca_stats))) {
  v_thr_min <- suppressWarnings(as.numeric(dca_stats$value[dca_stats$metric == "threshold_min"][1]))
  v_thr_max <- suppressWarnings(as.numeric(dca_stats$value[dca_stats$metric == "threshold_max"][1]))
  v_nb_max <- suppressWarnings(as.numeric(dca_stats$value[dca_stats$metric == "max_net_benefit"][1]))
  v_thr_at <- suppressWarnings(as.numeric(dca_stats$value[dca_stats$metric == "threshold_at_max_nb"][1]))
  if (!any(is.na(c(v_thr_min, v_thr_max, v_nb_max, v_thr_at)))) {
    dca_txt <- sprintf("图4显示，模型在阈值概率%.3f-%.3f范围内具有净获益，最大净获益为%.4f（阈值=%.3f）。", v_thr_min, v_thr_max, v_nb_max, v_thr_at)
  }
}
dca_full_txt <- paste0(dca_txt, "相关结果表：`03.dca_stats.csv`。")
doc <- body_add_fpar(doc,
  value = build_fpar(dca_full_txt, font_size = BODY_SIZE, first_line = 0.21),
  style = "Normal"
)
doc <- mk_blank(doc)

# ─────────────────────────────────────────
# 保存 DOCX
# ─────────────────────────────────────────
tmp_docx <- output_docx
print(doc, target = tmp_docx)

cat("\n报告已生成：", output_docx, "\n")
