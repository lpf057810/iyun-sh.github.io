# 初始化环境
rm(list = ls()); gc()


## 数据路径设置 ----
main_path <- "/media/desk16/iyunlyl/project/26.25YSH027F/"
num_path <- "08_Nomogram"
if (!dir.exists(main_path)) stop("错误：主路径不存在，请检查路径设置！")
work_path <- paste0(main_path, "/", num_path)
if (!dir.exists(work_path)) dir.create(work_path)
setwd(work_path)

## 加载基础包和自定义函数 ----
source("/media/desk16/iyunlyl/project/ikl_function.R")

set.seed(123)  # 设置随机种子以确保结果可重复

library(rms)

#run ----

gene <- read.csv("../07_Expression/05.gene.csv")$gene

name <- 'TCGA-LUAD'

if (! dir.exists(name)){
  dir.create(name)
}
df <- fread(paste0("../00_RawData/", name, ".tpm.csv")) %>% column_to_rownames("SYMBOL")
df.group <- fread(paste0("../00_RawData/", name, ".group.csv"))
comparison <-c("Normal",unique(df.group$group)[unique(df.group$group)!="Normal"])


tdf <- as.data.frame(t(df))
data <- tdf[df.group$sample,gene]
df.group$y <- ifelse(df.group$group == comparison[2],1,0)
#df.group$group <- ifelse(df.group$group == comparison[2],1,0)

d <- merge(data, df.group, by.x = "row.names", by.y = "sample")
ddist <- datadist(d)
options(datadist='ddist')

# 主程序
formula1 <- as.formula(paste0('y ~', paste(gene,sep = '', collapse = '+')))
fit1 <- lrm(formula1,  data = d, x = TRUE, y = TRUE,maxit=1000)
#fastbw(fit1, rule=c("aic"))

# write.csv(fit1$coefficients, file = paste0(name,"/01.lrm_coefficients.csv"))

regplot::regplot(fit1,  plots=c("boxplot","boxes"),
                 observation=F, title="Prediction Nomogram",
                 clickable=F, points=F,droplines=F,showP=F)


# nomogram 列线图
nom <- nomogram(fit1,  ##最佳模型
                fun=plogis, #进行logit转换 or function(x)1/(1+exp(-x))
                funlabel=paste0("Risk of ",comparison[2]),
                lp=F,  ##是否显示线性预测值
                #conf.int = F,##每个得分的置信度区间，用横线表示
                #abbrev = F,#是否用简称代表因子变量
                # fun.at=c(.1,.5,.9) ##风险坐标轴的范围和刻度
                fun.at=c(.001,.999) ##风险坐标轴的范围和刻度
                
)
plot(nom, cex.axis  = 1.2, cex.var = 1.7)


save_plot( paste0(name,"/01.nomogram.png"),
           plot(nom, cex.axis  = 1.2, cex.var = 1.7),
           height = 8, width = 10)

save_plot( paste0(name,"/01.nomogram.pdf"),
           plot(nom, cex.axis  = 1.2, cex.var = 1.7),
           height = 8, width = 10)


# 绘制列线图的校准曲线------
library(ResourceSelection)
dat2 <- d

cal1 <- calibrate(fit1, cmethod='KM', method='boot', B=1000)

hl1 <- ResourceSelection::hoslem.test(dat2$y,predict(fit1,dat2), g=10) 
hl1
hl2 <- stats::resid(fit1,"gof")
hl2
pval <- signif(hl2[5], 3)
pval    #0.391 需要大于0.05


#plot(cal1)

temp_cal1<- function(){
  par(mar = c(6,5,2,2))
  plot(cal1, lwd=2, lty=1,
       cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2,
       xlim=c(0, 1), ylim= c(0, 1),
       xlab="Predicted probability",
       ylab="Actual probability",
       col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
       legend=FALSE)
  lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
  abline(0, 1, lty=3, lwd=2)
  legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"),
         lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
  text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
  text(x = 0.4, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
}

save_plot(paste0(name,"/02.calibrate.png"),
          temp_cal1(),
          width=8,height=6)
save_plot(paste0(name,"/02.calibrate.pdf"),
          temp_cal1(),
          width=8,height=6)

# 决策曲线（DCA）----------------------------------------------------------------------

while (dev.cur() > 1) {dev.off()}
library(rmda)

model<- decision_curve(formula1,data=d)

dca_list <- lapply(gene, FUN = function(x){
  formula_str <- paste0("y~", x)
  gene_formula <- as.formula(formula_str)
  decision_curve(gene_formula, data=d)
})

temp_dca1 <- function(){
  plot_decision_curve(simple,
                      curve.names = 'model',
                      cost.benefit.axis = F,
                      confidence.intervals = F,
                      standardize = F,lty = c(1:3)
  )
}

temp_dca <- function(){
  plot_decision_curve(c(list(model),dca_list),
                      curve.names = c('model',gene),legend.position = "bottomleft",
                      cost.benefit.axis = F,
                      confidence.intervals = F,
                      standardize = F,lty = c(1,4:(length(dca_list)+4),2,3)
  )
}


save_plot(paste0(name,"/03.dca.png"),
          temp_dca()
)
save_plot(paste0(name,"/03.dca.pdf"),
          temp_dca()
)

# ROC---------------
library(pROC)
predicted<-predict(fit1,newdata = d)

#lrm_predict<-ifelse(predicted >0.5, comparison[2],comparison[1]) %>% as.factor()
#roc_curve = roc.curve(lrm_predict, weights.class0 =  d$Group == "DCM", curve = T)

roc_curve = roc(d$group,predicted)

temp_roc <- function(){
  par(pin = c(4,4), mar = c(6,6,6,1))
  plot(roc_curve, auc.main = T, legend = F, color = 'darkblue',
       xlab = "1-Specificity",
       print.auc=T,
       asp = 1 ,
       cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.main=2.0,   ##标题的缩放倍数
       main='model',
       col="#FF2E63",
       font.lab = 2,
       font.main = 2,
       font.sub =2)
}

save_plot(paste0(name,"/04.roc.pdf"),
          temp_roc()
)
save_plot(paste0(name,"/04.roc.png"),
          temp_roc()
)



# library(ggDCA)
# dca_data<- ggDCA:: dca(fit1,model.names =c(paste0(disese,' prediction nomogram')))
# 
# ## https://www.jianshu.com/p/a120f3f9ad78  
# ## 决策图美化
# library(ggprism)
# pdf(file=paste0(name,"/04.dca.pdf"), width =8, height = 7)
# par(pin = c(4,4), mar = c(6,6,6,1))
# ggplot(dca_data,linetype =T,lwd = 1.0)+
#   theme(legend.position="top")+
#   # scale_y_continuous(
#   #   limits = c(-0.01, 0.2),
#   #   guide = "prism_minor"
#   #   )+
#   scale_colour_prism(
#     palette = "floral",
#     # name = "Cylinders",
#     #label = c("模型1", "ALL", "None")
#   )+
#   theme(axis.title.x = element_text(size = 22, face = "bold"),
#         axis.title.y = element_text(size = 22, face = "bold"),
#         axis.text.x = element_text(size = 17, color='black',face = "bold",  vjust = 1, hjust = 1),
#         axis.text.y = element_text(size = 17, face = "bold"),
#         legend.text = element_text(size = 17,face = "bold"),
#         legend.title = element_blank(),
#         #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
#         text = element_text(family = "Times"),
#         panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()
# 
# png(file=paste0(name,"/04.dca.png"), width = 8, height = 7,units='in',res=600)
# par(pin = c(4,4), mar = c(6,6,6,1))
# ggplot(dca_data,linetype =T,lwd = 1.0)+
#   theme(legend.position="top")+
#   # scale_y_continuous(
#   #   limits = c(-0.01, 0.2),
#   #   guide = "prism_minor"
#   #   )+
#   scale_colour_prism(
#     palette = "floral",
#     # name = "Cylinders",
#     #label = c("模型1", "ALL", "None")
#   )+
#   theme(axis.title.x = element_text(size = 22, face = "bold"),
#         axis.title.y = element_text(size = 22, face = "bold"),
#         axis.text.x = element_text(size = 17, color='black',face = "bold",  vjust = 1, hjust = 1),
#         axis.text.y = element_text(size = 17, face = "bold"),
#         legend.text = element_text(size = 17,face = "bold"),
#         legend.title = element_blank(),
#         #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
#         text = element_text(family = "Times"),
#         panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()
# 

