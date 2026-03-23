# ==============================================================================
# 免疫浸润分析专用函数库
# ==============================================================================
# 说明: 本文件包含免疫浸润分析中常用的 R 函数
# 维护: 持续迭代优化，记录所有变更于 CHANGELOG.md
# 作者: [项目团队]
# 创建: 2026-02-15
# ==============================================================================

# 设置镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# ==============================================================================
# 字体设置
# ==============================================================================

# 设置字体以支持中文和英文
Sys.setenv(PATH = paste0(Sys.getenv('PATH'), ':/usr/share/fonts'))

# 设置基础图形字体
par(family = "Arial")

# 设置 ggplot2 主题和字体
theme_set(theme_bw() + theme(
  text = element_text(family = "Arial", size = 11),
  axis.text = element_text(family = "Arial", size = 10),
  axis.title = element_text(family = "Arial", size = 12),
  legend.text = element_text(family = "Arial", size = 10),
  legend.title = element_text(family = "Arial", size = 11)
))

# ==============================================================================
# 颜色定义
# ==============================================================================

#' 分组颜色
#'
#' 用于分组比较图的颜色方案
#' @export
group_color <- c("firebrick2", "#386CB0", "orange", "seagreen", "#BC80BD", "#17BECFFF")

#' 风险分组颜色
#' @export
risk_color <- c("#dc143c", "#6495ed")

#' 火山图颜色
#' @export
volcano_color <- c("seagreen", "darkgray", "firebrick2")

#' 渐变色
#' @export
ramp_color <- colorRampPalette(c("seagreen", "white", "firebrick2"))(100)

#' 调色板颜色
#' @export
palettes_color <- c(
  "#5f75ae", "#64a841", "#e5486e", "#de8e06", "#b5aa0f",
  "#7ba39d", "#b15928", "#6a3d9a", "#cab2d6", "#374E55FF",
  "#80796BFF", "#e31a1c", "#fb9a99", "#1f78b4", "#a6cee3",
  "#008280FF", "#8F7700FF", "#A20056FF", "#fdbf6f", "#E78AC3",
  "#b2df8a", "#CD534CFF", "#008B45FF", "#67001F", "#00A087FF",
  "#A73030FF", "#386CB0", "#F0027F", "#666666", "#EFC000FF",
  "#003C67FF", "#7AA6DCFF", "#8F7700FF", "#33a02c", "#66C2A5",
  "#A6D854", "#E5C494", "#6A3D9A", "#374E55FF", "#DF8F44FF",
  "#8DA0CB", "#80796BFF", "#FFFF99", "#E78AC3", "#7FC97F",
  "#3B3B3BFF", "#B24745FF", "#3B4992FF", "#631879FF", "#7AA6DCFF"
)

# ==============================================================================
# 保存图形函数
# ==============================================================================

#' 保存图形到文件
#'
#' 支持 PDF 和 PNG 格式，可同时保存两种格式
#'
#' @param filename 文件名（可带或不带扩展名）
#' @param plot 图形对象（ggplot、grob、绘图函数等）
#' @param outdir 输出目录，默认为当前目录
#' @param width 宽度，默认为 7
#' @param height 高度，默认为 7
#' @param both 逻辑值，是否同时保存 PDF 和 PNG，默认为 FALSE
#' @param bg 背景颜色，默认为 "white"
#' @param family 字体族，默认为 "Arial"
#' @param draw_func 绘图执行模式（1=自动识别，2=直接打印，3=表达式，4=原始模式）
#' @param units 单位，默认为 "in"
#' @param res 分辨率（DPI），默认为 600
#'
#' @return 无返回值，直接保存文件
#'
#' @export
save_plot <- function(filename, plot, outdir = ".",
                      width = 7, height = 7, both = FALSE,
                      bg = 'white', family = "DejaVu Sans",
                      draw_func = 1, units = "in", res = 300) {

  # 自动识别是否是 ggplot 对象
  is_gg <- inherits(plot, "ggplot")

  # 确保输出目录存在
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # 构造完整的文件路径（只用文件名，去掉目录和扩展名）
  base_name <- tools::file_path_sans_ext(basename(filename))
  file_path_pdf <- file.path(outdir, paste0(base_name, ".pdf"))
  file_path_png <- file.path(outdir, paste0(base_name, ".png"))

  # 绘图执行器
  draw_plot <- function() {
    if (is_gg) {
      print(plot)
    } else {
      if (is.function(plot)) {
        plot()
      } else {
        print(plot)
      }
    }
  }

  # 根据 draw_func 控制是否使用函数包装
  execute_plot <- function() {
    if (draw_func == 1) {
      draw_plot()
    } else if (draw_func == 2) {
      print(plot)
    } else if (draw_func == 3) {
      plot
    }
  }

  # 设置字体
  if (family != "sans") {
    par(family = family)
  }

  if (both) {
    # 同时保存 PNG 和 PDF
    # 保存 PNG (300dpi)
    png(filename = file_path_png, width = width, height = height, units = units, res = res)
    execute_plot()
    dev.off()
    cat("[OK] Saved PNG:", file_path_png, "\n")

    # 保存 PDF (矢量图，不指定字体)
    pdf(file = file_path_pdf, width = width, height = height)
    execute_plot()
    dev.off()
    cat("[OK] Saved PDF:", file_path_pdf, "\n")
  } else {
    # 默认只保存 PNG
    ext <- tolower(tools::file_ext(filename))
    if (ext == "png") {
      png(filename = file_path_png, width = width, height = height, units = units, res = res)
      execute_plot()
      dev.off()
      cat("[OK] File saved to:", file_path_png, "\n")
    } else if (ext == "pdf") {
      pdf(file = file_path_pdf, width = width, height = height, family = family)
      execute_plot()
      dev.off()
      cat("[OK] File saved to:", file_path_pdf, "\n")
    } else {
      # 其他格式也改为 PNG 保存
      png(filename = file_path_png, width = width, height = height, units = units, res = res)
      execute_plot()
      dev.off()
      cat("[OK] File saved to:", file_path_png, "\n")
    }
  }

  return(invisible(NULL))
}

# ==============================================================================
# 时间格式化函数
# ==============================================================================

#' 格式化当前时间
#'
#' @return 无返回值，打印格式化的时间
#' @export
formatted_time <- function() {
  cat(sprintf("[%s]\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

# ==============================================================================
# 相关性分析函数
# ==============================================================================

#' 计算两个数据集之间的相关性
#'
#' @param data1 数据矩阵1（行为特征，列为样本）
#' @param data2 数据矩阵2（行为特征，列为样本）
#' @param method 相关性方法，默认 "pearson"
#' @param cor_val 相关系数阈值，默认 0.6
#' @param onlycor 是否仅返回显著相关，默认 FALSE
#' @param drop 是否删除不显著的结果，默认 FALSE
#' @param plot 是否绑制热图，默认 TRUE
#'
#' @return 列表，包含 res（相关性结果）和 p（图形对象）
#'
#' @export
get_corr <- function(data1 = NULL, data2 = NULL, method = "pearson",
                     cor_val = 0.6, onlycor = FALSE, drop = FALSE, plot = TRUE) {

  # 加载必需的包
  if (!requireNamespace("ggcor", quietly = TRUE)) {
    stop("需要安装 ggcor 包")
  }
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("需要安装 WGCNA 包")
  }

  # 计算相关性和 p 值
  corr <- WGCNA::corAndPvalue(data1, data2, method = method)

  # 整理相关性结果
  cor_r2 <- corr$cor %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene1") %>%
    gather(key = "gene2", value = "Correlation", -gene1)

  cor_p2 <- corr$p %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene1") %>%
    gather(key = "gene2", value = "Pvalue", -gene1)

  cor_dat <- cbind(cor_r2, cor_p2)[, c("gene1", "gene2", "Correlation", "Pvalue")]
  cor_dat <- cor_dat[cor_dat$gene1 != cor_dat$gene2, ]

  # 添加显著性标记
  corr_plot <- cor_dat %>%
    dplyr::mutate(
      Pvalue.signif = dplyr::case_when(
        Pvalue < 0.001 ~ "***",
        Pvalue < 0.01 ~ "**",
        Pvalue < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  corr_plot$label <- paste0(round(corr_plot$Correlation, 2), corr_plot$Pvalue.signif)

  # 筛选显著相关
  if (onlycor) {
    corr_plot1 <- corr_plot %>%
      dplyr::filter(Pvalue.signif != "ns", abs(Correlation) >= cor_val)
    corr_plot <- corr_plot %>%
      dplyr::filter(gene2 %in% corr_plot1$gene2)
  }

  # 绑制热图
  p1 <- NULL
  if (plot && nrow(corr_plot) > 0) {
    p1 <- ggplot(corr_plot, aes(x = gene1, y = gene2)) +
      geom_tile(aes(fill = Correlation)) +
      geom_text(aes(label = label), color = "black", size = 4) +
      scale_fill_gradientn(
        colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
        limits = c(-1, 1),
        breaks = seq(-1, 1, by = 0.5)
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
  }

  return(list(res = corr_plot, p = p1))
}

# ==============================================================================
# 工具函数
# ==============================================================================

#' 创建目录（如果不存在）
#'
#' @param x 目录路径
#' @param recursive 是否递归创建，默认为 TRUE
#'
#' @return 无返回值
#' @export
dir_create <- function(x, recursive = TRUE) {
  if (!dir.exists(x)) {
    dir.create(x, recursive = recursive)
    cat(sprintf("[INFO] 创建目录: %s\n", x))
  }
}

#' 查看数据框头部
#'
#' @param dat 数据框
#' @export
head_view <- function(dat) {
  View(head(dat))
}

#' 格式化为科学计数法
#'
#' @param x 数值
#' @param digits 小数位数，默认为 1
#' @return 字符，科学计数法表示
#' @export
sci_number <- function(x, digits = 1) {
  sprintf(sprintf("%%.%de", digits), x)
}
