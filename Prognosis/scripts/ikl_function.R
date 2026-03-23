set.seed(123)

# if (!requireNamespace("pacman", quietly = TRUE)) {
#   install.packages("pacman")
# }
# library(pacman)
library(glue)
library(patchwork)
library(data.table)
library(tidyverse)
library(dplyr)
# library(qs)
library(qs2)
library(ggplot2)
library(rstatix)

library(future)
library(future.apply)
library(future.callr)


options(scipen = 5) 

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

# 颜色 ----

# > palettes(category = "random",4,show_col=F)
# >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
library(RColorBrewer)

palettes_color<- c(
  
  "#5f75ae","#64a841","#e5486e","#de8e06","#b5aa0f",
  "#7ba39d","#b15928","#6a3d9a","#cab2d6","#374E55FF",
  "#80796BFF","#e31a1c","#fb9a99","#1f78b4","#a6cee3",
  "#008280FF","#8F7700FF","#A20056FF","#fdbf6f","#E78AC3",
  "#b2df8a","#CD534CFF","#008B45FF","#67001F","#00A087FF",
  "#A73030FF","#386CB0","#F0027F","#666666","#EFC000FF",
  "#003C67FF","#7AA6DCFF","#8F7700FF","#33a02c","#66C2A5",
  "#A6D854","#E5C494","#6A3D9A","#374E55FF","#DF8F44FF",
  "#8DA0CB","#80796BFF","#FFFF99","#E78AC3","#7FC97F",
  "#3B3B3BFF","#B24745FF","#3B4992FF","#631879FF","#7AA6DCFF",
  "#7ba39d","#b15928","#00A1D5FF","#a6a6a6","#386CB0","#F0027F",
  "#1B9E77","#7570B3","#67001F","#4DBBD5FF","#F39B7FFF","#7FC97F",
  "#BEAED4","#224444","#DF8F44FF","#B24745FF","#3B4992FF",
  "#631879FF","#7AA6DCFF","#003C67FF","#8F7700FF","#3B3B3BFF",
  "#984EA3","#a6a6a6","#8DA0CB","#E78AC3","#FFD92F",
  "#8DD3C7","#1F78B4","#66A61E","#D62728FF","#9467BDFF",
  "#8C564BFF","#E377C2FF","#7F7F7FFF","#17BECFFF","#FB9A99",
  "#FDBF6F","#33adff","#439373","#92C5DE","#CAB2D6"
)


good_color = c( '#57C3F3','#53A85F', '#F1BB72',"#D62728FF","#8DD3C7","#FFD92F",
                '#F3B1A0', "#BC80BD",'#D6E7A3', '#E5D2DD', '#476D87',"#003C67FF")

temp_color <- c(
  "#1976B6","#F47E21","#2BA348","#D4242C","#9167AB",
  "#F878AF","#8DD3C7","#BCBD31",'#57C3F3', '#F1BB72',
  "#FFD92F","#AA554A","#7F7F7F", '#C3B1A0', '#476D87',
  '#D6E7A3', '#E1A2DD',"#003C67FF","#7ba39d","#F3B1A0",
  "#631879FF","#B24745FF","#A6D854","#E5C494","#6A3D9A"
)



ramp_color<- colorRampPalette(c("seagreen", "white", "firebrick2"))(100)

group_color<- c("firebrick2", "#386CB0","orange","seagreen","#BC80BD","#17BECFFF")
risk_color<- c("#dc143c", "#6495ed")
volcano_color <- c("seagreen", "darkgray", "firebrick2")

# save_plot ----
save_plot <- function(filename, plot, outdir = ".", 
                      width = 7, height = 7, both = FALSE,
                      bg = 'white', family = "Arial", 
                      draw_func = 1, units = "in", res = 600) {
  
  # 自动识别是否是 ggplot 对象
  is_gg <- inherits(plot, "ggplot")
  
  # 确保输出目录存在
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # 构造完整的文件路径（去掉扩展名）
  base_name <- tools::file_path_sans_ext(filename)
  file_path_pdf <- file.path(outdir, paste0(base_name, ".pdf"))
  file_path_png <- file.path(outdir, paste0(base_name, ".png"))
  
  # 绘图执行器
  draw_plot <- function() {
    if (is_gg) {
      print(plot)  # 如果是 ggplot 对象，直接打印
    } else {
      if (is.function(plot)) {
        plot()  # 如果是绘图函数，执行该函数
      } else {
        print(plot)  # 否则直接打印传入的绘图对象
      }
    }
  }
  
  # 根据 draw_func 控制是否使用函数包装
  execute_plot <- function() {
    if (draw_func ==1 ) {
      draw_plot()
    }
    if (draw_func ==2 ){
      print(plot)  # 直接打印绘图对象
    }
    if (draw_func ==3 ){
      plot  # 直接打印绘图对象
      
    }
  }
  if (draw_func == 4){
    if (both) {
      # PDF 保存
      grDevices::cairo_pdf(file = file_path_pdf, width = width, height = height, bg = bg, family = family)
      plot
      dev.off()
      
      # PNG 保存
      png(file = file_path_png, width = width, height = height, 
          bg = bg, units = units, res = res, family = family)
      plot
      dev.off()
      
      cat("- Saved both PDF and PNG:\n", file_path_pdf, "\n", file_path_png, "\n")
      
    } else {
      
      ext <- tools::file_ext(filename)
      full_path <- file.path(outdir, filename)
      
      if (ext == "pdf") {
        grDevices::cairo_pdf(file = full_path, width = width, height = height, bg = bg, family = family)
        plot
        dev.off()
        
      } else if (ext == "png") {
        png(file = full_path, width = width, height = height, 
            bg = bg, units = units, res = res, family = family)
        plot
        dev.off()
        
      } else {
        stop("Unsupported file format. Use .pdf or .png")
        }
      }
  }
  if (both) {
    # PDF 保存
    grDevices::cairo_pdf(file = file_path_pdf, width = width, height = height, bg = bg, family = family)
    execute_plot()
    dev.off()
    
    # PNG 保存
    png(file = file_path_png, width = width, height = height, 
        bg = bg, units = units, res = res, family = family)
    execute_plot()
    dev.off()
    
    cat("- Saved both PDF and PNG:\n", file_path_pdf, "\n", file_path_png, "\n")
    
  } else {
    
    ext <- tools::file_ext(filename)
    full_path <- file.path(outdir, filename)
    
    if (ext == "pdf") {
      grDevices::cairo_pdf(file = full_path, width = width, height = height, bg = bg, family = family)
      execute_plot()
      dev.off()
      
    } else if (ext == "png") {
      png(file = full_path, width = width, height = height, 
          bg = bg, units = units, res = res, family = family)
      execute_plot()
      dev.off()
      
    } else {
      stop("Unsupported file format. Use .pdf or .png")
    }
    
    cat("- File saved to:", full_path, "\n")
  }
}


# head_view ----
head_view <- function(dat){View(head(dat))}
# venn_plot ----
# venn_plot <- function(data,save=FALSE,name="01.venn.pdf",width=8,height=7){
#   library(ggvenn)
#   color <-  group_color
#   suppressWarnings({
#     p<- ggvenn(data,
#                fill_color =color[1:length(data)],
#                set_name_size = 6,
#                # show_percentage = FALSE,
#                stroke_alpha = 0.5,
#                stroke_size = 0.3,
#                text_size = 5,
#                stroke_color="black",
#                stroke_linetype="solid",
#                set_name_color=color[1:length(data)],
#                text_color = 'black')
#   })
#   if (save){
#     save_plot(name,p,width=width,height=height,both=T)
#   }else{
#     return(p)
#   }
# }

venn_plot <- function(data, save = FALSE, name = "01.venn.pdf", width = 8, height = 7,
                      set_name_size=6,text_size = 5,mult=.2,num_sets=NULL) {
  # 获取数据集的数量
  if(is.null(num_sets)){num_sets <- length(data)}
  # 使用预定义的颜色组，假设 group_color 在环境中已定义
  color <- group_color
  
  if (num_sets >= 2 && num_sets <= 4) {
    # 当数据集数量在2到4之间时，使用 ggvenn 包
    library(ggvenn) # 加载 ggvenn 包
    suppressWarnings({ # 抑制警告信息
      p <- ggvenn(data,
                  fill_color = color[1:num_sets],       # 填充颜色
                  set_name_size = set_name_size,                     # 集合名称字体大小
                  stroke_alpha = 0.5,                    # 描边透明度
                  stroke_size = 0.3,                     # 描边粗细
                  text_size = text_size,                         # 文本（交集数量）字体大小
                  stroke_color = "black",                # 描边颜色
                  stroke_linetype = "solid",             # 描边线型
                  set_name_color = color[1:num_sets],    # 集合名称颜色
                  text_color = 'black')+                  # 文本（交集数量）颜色
        scale_x_continuous(expand = expansion(mult = mult))
    })
  } else if (num_sets > 4 && num_sets < 7) {
    # 当数据集数量在5到6之间时，使用 ggVennDiagram 包
    library(ggVennDiagram) # 加载 ggVennDiagram 包
    suppressWarnings({ # 抑制警告信息
      p<- venn::venn(mydata,ilabels = "counts",
                     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜色
                     opacity = 0.3,  # 调整颜色透明度
                     box = F,        # 是否添加边框
                     ilcs = 1,     # 数字大小
                     sncs = 1,ggplot = T,borders = T
      )
      
    })
  }else {
    stop("数据集的数量必须是2个或更多。")
  }
  
  # 如果是 ggvenn 或 ggVennDiagram 生成的图，则按常规方式保存或返回图表对象
  if (save && num_sets < 7) {
    save_plot(name, p, width = width, height = height, both = TRUE)
    intersect_res<- Reduce(intersect,data)
    return(intersect_res)
  } else if (!save && num_sets < 7) {
    return(p)
  }
}

# fast_geo ----
fast_geo <- function(geo_series_id,num_threads = 8,annot=T,simpd=F){
  options(timeout=100000000)
  temp_function<- function(x) {
    tryCatch(
      {
        gset <- tinyarray::geo_download(x,destdir = "/media/desk16/iyunlyl/project/GEO/",simpd = simpd)
        if(annot){
          gpl<-AnnoProbe::idmap(gset$gpl,destdir = "/media/desk16/iyunlyl/project/GEO/")
          rownames(gpl) <- gpl$probe_id
        }else{
          gpl<-gset$gpl
        }
        list(
          gset = gset,
          expr = gset$exp,
          pd = gset$pd,
          gpl = gpl)
      },error=function(e){print("Error in getGEO for ", x, ": ", e$message)}
    )
  }
  # 使用 mclapply 实现多线程处理
  geo_data_list <- parallel::mclapply(geo_series_id, temp_function,
                                      mc.cores = num_threads)
  names(geo_data_list) <- geo_series_id
  return(geo_data_list)
}

# get_geo ----
# get_geo <- function(geo_series_id, num_threads = 8,getGPL = T) {
#   library(GEOquery)
#   library(parallel)  #  用于并行处理
#   options(timeout = 100000000)
#   
#   # 使用 mclapply 实现多线程处理
#   geo_data_list <- mclapply(geo_series_id, function(x) {
#     tryCatch(
#       {
#         gset <- getGEO(x, GSEMatrix = TRUE,getGPL = getGPL, destdir = "/media/desk16/iyunlyl/project/GEO/")
#         dat_type <- gset[[1]]@experimentData@other$type
#         message(dat_type)
#         path <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", x)
#         d_path <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", x)
#         tpm_path <- ""
#         if (grepl("high throughput", dat_type)) {
#           urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
#           temp_file <- paste0(x, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
#           path <- paste0(urld, "&acc=", x, "&file=", temp_file)
#           
#           urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
#           temp_file1 <- paste0(x, "_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
#           path <- paste0(urld, "&acc=", x, "&file=", temp_file1)
#           
#           # count 
#           if (!file.exists(temp_file)) {
#             tryCatch({
#               curl::curl_download(path, destfile = temp_file, quiet = FALSE)
#               path <- temp_file
#             }, error = function(e) {
#               print(e$message)
#               path <<- e$message
#             })
#           } else {
#             path <- temp_file
#           }
#           
#           # tpm
#           if (!file.exists(temp_file1)) {
#             tryCatch({
#               curl::curl_download(path, destfile = temp_file1, quiet = FALSE)
#               tpm_path <- temp_file1
#             }, error = function(e) {
#               print(e$message)
#               path <<- e$message
#             })
#           } else {
#             tpm_path <- temp_file1
#           }
#         }
#         
#         if (length(gset) > 1) {
#           list(
#             gset = gset,
#             path = path,
#             d_path = d_path
#           )
#         } else {
#           list(
#             gset = gset,
#             expr = Biobase::exprs(gset[[1]]),
#             pd = Biobase::pData(gset[[1]]),
#             gpl = Biobase::fData(gset[[1]]),
#             path = path,
#             d_path = d_path,
#             tpm_path=tpm_path
#           )
#         }
#       },
#       error = function(e) {
#         print(paste("Error in getGEO for", x, ":", e$message))
#         return(NULL)
#       }
#     )
#   }, mc.cores = num_threads)  
#   
#   names(geo_data_list) <- geo_series_id
#   return(geo_data_list)
# }

get_geo_url <- function(geo_series_id,num_threads=1) {
  temp_func <- function(x){
    # 构造 series 目录，例如 GSE161nnn
    GEO <- toupper(x)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    
    if(grepl("GSE",GEO)){
      stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
      gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
      url <- sprintf(gdsurl, stub, GEO)
    }else{
      gplurl <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/"
      url <- sprintf(gplurl, stub, GEO)
    }
    
    temp_lines <- readLines(url)
    # 找到包含 .txt.gz 文件的那一行
    file_line <- temp_lines[grep("\\.gz", temp_lines)]
    # 使用正则表达式提取大小（如 4.6M）
    size_match <- trimws(gsub(".*\\.gz.*?([0-9.]+[KMGT]).*", "\\1", file_line, ignore.case = TRUE))
    message(paste0(url,",size:",size_match))
    return(url)
  }
  
  # 
  # # 定义清理函数：主进程退出时关闭multisession
  # on.exit({
  #   plan(sequential)  # 关键：自动清理子进程，避免残留
  # })
  # 
  # 启动 future 并行
  # library(future)
  # library(future.apply)
  # library(future.callr)
  # install.packages("future.callr")
  # plan(multisession, workers = 8)
  # plan(multicore, workers = 8)
  # plan(callr, workers = 8)
  
  # 配置future后端：优先multicore（类Unix），否则multisession（Windows）
  plan(multisession, workers = num_threads)
  
  # 使用future并发处理
  geo_data_list <- future_lapply(
    geo_series_id,
    temp_func,
    future.stdout = TRUE,
    future.conditions = c("message", "warning", "error")
  )
  names(geo_data_list) <- geo_series_id
  return(geo_data_list)
  
}



get_geo <- function(geo_series_id, num_threads = 16,by_geochina=T,force=F, destdir = "/media/desk16/iyunlyl/project/GEO/") {
  options(timeout = 100000000)
  library(AnnoProbe)
  library(GEOquery)
  temp_function <- function(x) {
    # 预定义，避免 not found 错误
    temp_file <- NULL
    temp_file1 <- NULL
    tpm <- NULL
    gset <- NULL
    
    tryCatch({

      if(force){
        file.remove(paste0(destdir, "/", x, "_eSet.qs2"))
        
      }
      
      # 如果未下载，尝试使用 geoChina
      if (!file.exists(paste0(destdir, "/", x, "_eSet.qs2"))) {
        
        if(file.exists(paste0(destdir, "/", x, "_eSet.Rdata"))){
          load(paste0(destdir, "/", x, "_eSet.Rdata"))
          file.remove(paste0(destdir, "/", x, "_eSet.Rdata"))
          qs_save(gset,paste0(destdir, "/", x, "_eSet.qs2"))
        }
        
        if (by_geochina){
          gset <- tryCatch({
            gset<- AnnoProbe::geoChina(x, destdir = destdir)
            file.remove(paste0(destdir, "/", x, "_eSet.Rdata"))
            qs_save(gset,paste0(destdir, "/", x, "_eSet.qs2"))
            gset
          }, error = function(e) {
            NULL
          })
        }
        
        # geoChina失败 → 使用 getGEO
        if (!is.list(gset)) {
          
          gset <- tryCatch({
          message(sprintf("%s not indexed by AnnoProbe. Downloading from GEO...", x))
          temp_url<- get_geo_url(x)
          gset<- getGEO(x, GSEMatrix = TRUE, getGPL = TRUE, destdir = destdir)
          gpl<- Biobase::fData(gset[[1]])
          qs_save(gset,paste0(destdir, "/", x, "_eSet.qs2"))
          qs_save(gpl,paste0(destdir, "/", gset[[1]]@annotation, "_bioc.qs2"))
          gset
          }, error = function(e) {
            NULL
          })
        }
      } else {
        gset<- qs_read(paste0(destdir, "/", x, "_eSet.qs2"))
      }

      if (!file.exists(paste0(destdir, "/", gset[[1]]@annotation, "_bioc.qs2"))) {
        gpl <- tryCatch({
          AnnoProbe::idmap(gset[[1]]@annotation,destdir = destdir)
          file.remove(paste0(destdir, "/", x, "_bioc.rda"))
        }, error = function(e) {
          NULL
        })
        rownames(gpl) <- gpl$probe_id
        if (is.null(gpl)) {
          gpl <- getGEO(gset[[1]]@annotation,destdir = destdir)
          gpl <- Table(gpl)
        }
        qs_save(gpl,paste0(destdir, "/", gset[[1]]@annotation, "_bioc.qs2"))
        
      }else{
        gpl<- qs_read(paste0(destdir, "/", gset[[1]]@annotation, "_bioc.qs2"))
        
      }
      
      # 实验类型
      dat_type <- gset[[1]]@experimentData@other$type
      message(sprintf("[%s] Data type: %s", x, dat_type))
      
      # RNA-seq 下载处理
      if (grepl("high throughput", dat_type)) {
        
        urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
        
        # raw count
        temp_file <- paste0(destdir, "/", x, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
        path <- paste0(urld, "&acc=", x, "&file=", temp_file)
        
        if (!file.exists(temp_file)) {
          tryCatch({
            curl::curl_download(path, destfile = temp_file, quiet = FALSE)
          }, error = function(e) {
            message(sprintf("Failed to download raw count for %s: %s", x, e$message))
          })
        }
        # expr <- data.table::fread(temp_file)
        
        # TPM
        temp_file1 <- paste0(destdir, "/", x, "_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
        path1 <- paste0(urld, "&acc=", x, "&file=", temp_file1)
        
        if (!file.exists(temp_file1)) {
          tryCatch({
            curl::curl_download(path1, destfile = temp_file1, quiet = FALSE)
          }, error = function(e) {
            message(sprintf("Failed to download TPM for %s: %s", x, e$message))
          })
        }
        # tpm <- data.table::fread(temp_file1)
        
      }
      
      # 返回内容
      if (!is.null(temp_file)){
        path = temp_file
        }else{
          path = paste0(destdir, "/", x, "_eSet.qs2")
          }
     list(
          gset = gset,
          expr = Biobase::exprs(gset[[1]]),
          pd = Biobase::pData(gset[[1]]),
          gpl = gpl,
          path = path,
          tpm_path = temp_file1
        )
    }, error = function(e) {
      message(sprintf("Error in getGEO for %s: %s", x, e$message))
      return(NULL)
    })
  }
  
  
  # 定义清理函数：主进程退出时关闭multisession
  on.exit({
    plan(sequential)  # 关键：自动清理子进程，避免残留
    message("已清理并行子进程")
  })
  
  # 启动 future 并行
  # library(future)
  # library(future.apply)
  # library(future.callr)
  # install.packages("future.callr")
  # plan(multisession, workers = 8)
  # plan(multicore, workers = 8)
  # plan(callr, workers = 8)
  
  # 配置future后端：优先multicore（类Unix），否则multisession（Windows）
  plan(multisession, workers = num_threads)
  
  # 使用future并发处理
  geo_data_list <- future_lapply(
    geo_series_id,
    temp_function,
    future.stdout = TRUE,
    future.conditions = c("message", "warning", "error")
  )

  names(geo_data_list) <- geo_series_id
  return(geo_data_list)

  
  # # 使用 mclapply 实现多线程处理
  # geo_data_list <- parallel::mclapply(geo_series_id, temp_function,
  #                                     mc.cores = num_threads)
  # names(geo_data_list) <- geo_series_id
  # return(geo_data_list)
  
  
  
}


get_data<- function(geo,sym_col=NULL,onlybitr=F,tpm=T,annot_file="/media/desk16/iyunlyl/project/Human.GRCh38.p13.annot.tsv.gz"){
  # 如果是高通量数据则读取raw_count,并转换ID
  if (nrow(geo$expr) == 0){
    if(tpm){
      geo$expr <- fread(geo$tpm_path,fill = T)
    }else{
      geo$expr <- fread(geo$path,fill = T)
    }
    }else{print(geo$path)}
  
  if (nrow(geo$gpl) == 0){
    temp_id_col <- colnames(geo$expr)[1]
    symbol <- NULL
    if(file.exists(annot_file) & onlybitr==F){
      annot<- fread(annot_file)
      symbol<- annot %>%select(c(temp_id_col,"Symbol"))%>%as.data.frame()
      colnames(symbol) <- c("ENTREZID","SYMBOL")
      no_annot<- geo$expr[[temp_id_col]][!geo$expr[[temp_id_col]] %in% symbol$ENTREZID]
    }
    
    if(onlybitr){no_annot<- geo$expr[[temp_id_col]]}
    if(length(no_annot)!=0){
      
      # 需要判断 temp_id_col 是什么类型的 根据类型更改fromType
      symbol1<- clusterProfiler::bitr(no_annot,
                                     fromType = "ENTREZID", # ENSEMBL
                                     toType = "SYMBOL",
                                     OrgDb = "org.Hs.eg.db")
      symbol <- rbind(symbol,symbol1)%>%unique()
    }
    
    expr_annot<- merge(symbol,geo$expr,by=1)
    
    expr_annot<- expr_annot %>% dplyr::select(-c("ENTREZID"))
    expr_annot <- expr_annot[expr_annot$SYMBOL!="",] %>% as.data.table()
    expr_annot <- expr_annot[,lapply(.SD, mean), by = "SYMBOL"]
    # expr_annot <- expr_annot[,lapply(.SD, round), by = "SYMBOL"]
  }else{
    ### 通过gpl注释表达矩阵
    if (is.null(sym_col)){
      sym_col<- colnames(geo$gpl)[grepl("SYMBOL",toupper(colnames(geo$gpl)))]
    }else{
      sym_col <- sym_col
    }
    if(length(sym_col)==0){
      message("no symbol")
      return(geo$expr)
      
    }
 
    if (all(grepl("^\\d+$", rownames(geo$gpl))) &  "ID" %in% colnames(geo$gpl) ) {
      # 如果所有行名都是由数字组成的字符串，则执行相应的操作
      rownames(geo$gpl) <- geo$gpl$ID
    }
    geo_annot<- geo$gpl[sym_col]
    expr_annot <- merge(geo_annot,geo$expr,by=0)
    expr_annot <- expr_annot[expr_annot[[sym_col]] !="",] %>% dplyr::select(-"Row.names") %>% as.data.table()
    colnames(expr_annot) <- gsub(sym_col,"SYMBOL",colnames(expr_annot))
    expr_annot <- expr_annot[,lapply(.SD, mean), by = "SYMBOL"]
  }
  return(expr_annot)
}

# 表达量验证 ----
get_exp_wilcox<- function(geo_expr,geo_group,gene,Control1="Control"){
  
  hubexp <- geo_expr[rownames(geo_expr) %in% gene, ]
  # 转置表达矩阵并合并分组信息
  hubexp <- t(hubexp) %>% as.data.frame()
  hubexp <- merge(hubexp, geo_group, by.x = 0, by.y = 1)
  rownames(hubexp) <- hubexp$Row.names
  colnames(hubexp)[colnames(hubexp) == "Row.names"] <- "sample"
  
  # 数据重塑为长格式
  boxdat <- reshape2::melt(hubexp, id.vars = c("sample", "group"))
  boxdat$value <- as.numeric(boxdat$value)
  boxdat$variable <- as.character(boxdat$variable)
  boxdat <- boxdat[order(boxdat$group, decreasing = TRUE), ]
  
  
  # 计算统计结果
  stat_res <- boxdat %>%
    na.omit() %>%
    group_by(variable) %>%
    wilcox_test(value ~ group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p") %>%
    as.data.frame()
  
  
  # 计算平均表达值
  temp_dat <- as.data.table(boxdat)
  mean_expr <- temp_dat[, .(Mean_Expression = mean(value)), by = .(variable, group)]
  mean_expr_wide <- dcast(mean_expr, variable ~ group, value.var = "Mean_Expression") %>%as.data.table()
  
  # 动态计算差值
  
  
  comparison <- c(Control1, unique(geo_group$group)[unique(geo_group$group) != Control1])
  
  mean_expr_wide <- mean_expr_wide[, Difference := get(comparison[2]) - get(comparison[1])]
  res<- merge(stat_res,mean_expr_wide[,c('variable','Difference')],by="variable")
  
  return(list(res=res,boxdat=boxdat))
}


plot_expression <- function(plot_data, expression_data,geo_name="",Control1="Control") {
  plot_data <- plot_data[plot_data$variable %in% as.character(expression_data$variable), ]
  plot_data <- arrange(plot_data, plot_data$variable)
  comparison <- c(Control1, unique(plot_data$group)[unique(plot_data$group) != Control1])
  plot_data$group <- factor(plot_data$group,levels = comparison)
  
  p <- ggplot(plot_data, aes(x = variable, y = value, fill = group)) +
    geom_boxplot(width = 0.5, alpha = 0.8, position = position_dodge(0.9), outlier.shape = NA) +
    scale_fill_manual(values =rev(group_color[1:2]))+
    # scale_fill_hue()+
    # scale_fill_manual(values =c('#33CCD0','#FA918A'))+
    # stat_compare_means(
    #   comparisons = list(c("T","N")),
    #   label = "p.signif")+
    annotate(geom = "text", x = as.character(expression_data$variable),
             y = max(plot_data$value), size = 5, 
             label = as.character(expression_data$p.signif)) +
    theme_classic() +
    theme(legend.position = "top",
          legend.title = element_text(face="bold"),
          legend.text = element_text(face="bold"),
    ) +
    theme(axis.title.x = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(angle = 90, size = 10,vjust = 0.5, hjust = 1, face = "bold"),
          axis.title.y = element_text(size = 15,face = "bold"),
          axis.text.y = element_text(size = 15,  face = "bold")) +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 14)) +
    labs(x = geo_name, y = "Expression") 
  return(p)
}

# plot_roc ----
plot_roc<- function(df,df.group,comparison,gene,title1="ROC"){
  
  library(ggplot2)
  library(pROC)
  library(patchwork)
  df <- as.data.frame(t(df))
  
  roc_obj_list <- list()
  for (hubgene in gene) {
    roc_obj <- roc(df.group$group, df[,hubgene], levels = comparison)
    temp<- paste0(hubgene,"(",round(auc(roc_obj),digits = 3),")")
    roc_obj_list[[temp]] <- roc_obj
    # pdf(paste0("0", i, ".", name, '_', colnames(hub_exp)[i], ".pdf"))
    # plot(roc_obj,
    #      print.auc=T,
    #      print.auc.x=0.4,print.auc.y=0.5,
    #      #auc.polygon=T,
    #      #auc.polygon.con="#fff7f7",
    #      grid=c(0.5,0.2),
    #      grid.col=c("black","black"),
    #      #print.thres=T,
    #      main=colnames(hub_exp)[i],
    #      col="#FF2E63",
    #      legacy.axes=T)
    # dev.off()
    
  }
  p<- ggroc(roc_obj_list,)+
    theme_light() + 
    ggtitle(title1) +
    scale_color_discrete(name="AUC")+
    geom_segment(data = data.frame(),aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color = "grey", linetype = "dashed") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          text = element_text(hjust = 0.5, face = 'bold'))
  return(p)
  
  
}

# 相关性分析 ----
get_corr <- function(data1=NULL,data2=NULL,method="pearson",cor_val=0.6,onlycor=F,drop=F,plot=T){
  library(ggcor)
  corr<- WGCNA::corAndPvalue(data1,data2,method=method)
  
  cor_r2 <- corr$cor %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene1") %>%
    gather(.,key = "gene2", value = "Correlation", -gene1)
  cor_p2 <- corr$p  %>% as.data.frame() %>% tibble::rownames_to_column(var = "gene1") %>%
    gather(.,key = "gene2", value = "Pvalue", -gene1)
  cor_dat <- cbind(cor_r2, cor_p2)[,c("gene1","gene2","Correlation","Pvalue")]
  cor_dat<- cor_dat[cor_dat$gene1 != cor_dat$gene2 ,]
  corr_plot<- cor_dat %>% add_significance("Pvalue") 
  corr_plot$label<- paste0(round(corr_plot$Correlation,2),corr_plot$Pvalue.signif)
  
  if(onlycor){
    corr_plot1<- corr_plot %>% filter(Pvalue.signif != "ns",abs(Correlation)>=cor_val)
    corr_plot <- corr_plot %>%filter(corr_plot$gene2 %in% corr_plot1$gene2)
    }
  if(plot){
  p1<- ggplot(corr_plot, aes(x=gene1,y=gene2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=label), color="black", size=4) + # 把星号添加进去
    scale_fill_gradientn(
      colors= colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      limits=c(-1,1),
      breaks = seq(-1, 1, by = 0.5),
      # name=paste0("*    p < 0.05","\n\n",
      #             "**  p < 0.01","\n\n",
      #             "*** p < 0.001","\n\n",
      #             "Correlation") # 把P值添加到图例
    ) +
    labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
    theme(axis.text.x = element_text(size=12,angle = 45,hjust = 1,color = "black"),
          axis.text.y = element_text(size=12,color = "black"),
          axis.ticks  = element_blank(),
          panel.background=element_blank()) # 做一些简单的美化
  
  message("max abs Correlation:",max(abs(cor_dat$Correlation)))
  plot_dat<- cor_dat %>% 
    mutate(
      rd = case_when(
        Correlation <= -cor_val ~ paste0("<= -",cor_val), # 小于等于 -0.6
        Correlation > -cor_val & Correlation < cor_val ~ paste0("-",cor_val," - ",cor_val), # 在 -0.6 到 0.6 之间
        Correlation >= cor_val ~ paste0(">= ",cor_val), # 大于等于 0.6
        TRUE ~ NA_character_ # 处理特殊情况
      ),
      pd = case_when(
        Pvalue <= 0.05 ~ "<= 0.05", 
        Pvalue > 0.05 ~ "> 0.05", 
        TRUE ~ NA_character_ # 处理特殊情况
      ),
      rpd=case_when(
        Pvalue <= 0.05 & Correlation <= -cor_val ~ 'bold',
        Pvalue <= 0.05 & Correlation >= cor_val ~ 'bold',
        Correlation > -cor_val & Correlation < cor_val ~ 'no_bold',
        TRUE ~ NA_character_ # 处理特殊情况
      )
    )
  if(onlycor){
    plot_dat1<- plot_dat %>% filter(rpd == "bold")
    data2 <- data2[,plot_dat1$gene2]
    plot_dat <- plot_dat%>%filter(plot_dat$gene2 %in% plot_dat1$gene2)
    }
  
  plot_dat$rd <- factor(plot_dat$rd,
                        levels= c(paste0("<= -",cor_val),
                                  paste0("-",cor_val," - ",cor_val),
                                  paste0(">= ",cor_val))
                        )
  plot_dat$pd <- factor(plot_dat$pd,levels= c("<= 0.05","> 0.05"))
  plot_dat$rpd <- factor(plot_dat$rpd,levels= c('bold',"no_bold"))
  
  if(!is.null(data2)){
  if (ncol(data2)>15 ){
    temp_size <- 8
    }else{
      temp_size <- 10
      }
  }
  
  col_map <- setNames(
    c("#4575B4", "#708090", "#D95F02"),
    c(paste0("<= -", cor_val),
      paste0("-", cor_val, " - ", cor_val),
      paste0(">= ", cor_val))
  )
  p_map <- setNames(
    c("solid","dotted"),
    c("<= 0.05","> 0.05")
  )
  
  p<- quickcor(data2, type = "upper") +
    geom_square() +
    anno_link(aes(colour = rd,
                  size = rpd,
                  linetype=pd,
                  alpha=pd, 
                  ), data = plot_dat,nudge_x = 0.5) +
    scale_fill_gradientn(
      colors= colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      breaks = seq(-1, 1, by = 0.5),
      limits=c(-1,1)
    )+
    scale_alpha_manual(values = c(0.7, 0.7),drop = FALSE)+  # 是否小于0.05 透明度
    scale_size_manual(values = c(1,0.3),drop = FALSE) + # 是否abs大于0.6 加粗
    # scale_linetype_manual(values = c("solid",'dashed','dotted'),drop = FALSE) +
    scale_linetype_manual(values = p_map,drop = drop) +
    scale_colour_manual(values = col_map, drop = drop)+
    theme(axis.text = element_text(size = temp_size))+
    guides(colour = guide_legend(title = paste0(tools::toTitleCase(method),"'s r"),
                                 # override.aes = list(
                                 #   colour = unname(col_map),   # 每个 key 的颜色
                                 #   linetype = "solid",
                                 #   size = 1.5,                 # 对点/线有效，视 ggplot 版本也可试 linewidth = 1.5
                                 #   alpha = 1
                                 # ),     没用 
                                 order = 2),
           # alpha = guide_legend(title = paste0(tools::toTitleCase(method),"'s p"),
           #                       # override.aes = list(size = 3),
           #                       order = 1),
           alpha="none",
           linetype = guide_legend(title = paste0(tools::toTitleCase(method),"'s p"),
                                # override.aes = list(size = 3),
                                order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3),
           size = "none"
           )
  }else{p=NULL;p1=NULL;plot_dat=NULL}
  
  # library(ggrepel)
  # dat<- cor_dat[cor_dat$gene1 == temp_gene,]
  # # 对相关系数和p值转换为分类变量
  # dat$cor1 <- cut(abs(dat$Correlation),# 绝对值
  #                 breaks = c(0, 0.3, 0.5, 0.7, 1),
  #                 labels = c("< 0.3","0.3 - 0.5","0.5 - 0.7","> 0.7"),
  #                 right=FALSE) # right=FALSE表示表示区间为左闭右开
  # dat$pvalue1 <- cut(dat$Pvalue,
  #                    breaks = c(0, 0.001, 0.01, 0.05, 1),
  #                    labels = c("< 0.001","< 0.01","< 0.05","> 0.05"),
  #                    right=FALSE) 
  # # 排序
  # dat = dat[order(dat$Correlation),]
  # dat$gene2 = factor(dat$gene2, levels = dat$gene2)
  # p1 = ggplot(dat, aes(x = Correlation, y = gene2, color = pvalue1)) +
  #   scale_size_manual(values = c(0.1,0.5,1,1.5, 2),drop = FALSE) +
  #   scale_color_manual(name="pvalue",
  #                      values = c("#E69F00", "#56B4E9", "#009E73", "gray"),drop=F)+
  #   geom_segment(aes(x = 0, y = gene2, xend = Correlation, yend = gene2),size = 1) +
  #   geom_point(aes(size = cor1))+
  #   geom_text_repel(aes(label=round(Pvalue,3)),nudge_y =0.1,show.legend = FALSE)+
  #   theme_bw()+
  #   labs(size = "Cor")
  # 
  
  return(list(corr=corr,res=cor_dat,p=p,p1=p1,plot_dat=plot_dat))
}

# dir_create ----
dir_create<- function(x,recursive=F){if (!dir.exists(x)) dir.create(x,recursive = recursive)}
# sci_number ----
sci_number<- function(x) sprintf("%.1e", x)  # 自定义科学计数法格式

# formatted_time ----
formatted_time <- function(){message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))}

# try_run ----
try_run <- function(abc){tryCatch({abc},error=function(e){formatted_time();print(e$message)})}

