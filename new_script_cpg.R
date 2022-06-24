# Title     : TODO
# Objective : TODO
# Created by: Guest0001
# Created on: 2021/9/2
# 文件主要是ratio 动态变化分析， 组学数据的质量控制及总体概览

### DEP_ratio Deseq2
#### "2021-09-03 11:25:09 CST"
#### V1

rm(list=ls())
library(readr)
library(DESeq2)
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
cts00 <- read.csv("D:/R_protemics/data/24to0ratio.CSV",row.names = 1)
#cts00 = as.numeric(cts00)
cts <- round(cts00)
mutant <- factor(rep(c('T26_4','15_102'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","15_102","T26_4"))
write.csv(as.data.frame(res01), file='res_DEPratio24to0.csv')
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("mutant"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=mutant)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_x_continuous(limits=c(-15, 15))+
  scale_y_continuous(limits=c(-15, 15))+
  coord_fixed()
ggsave("pro_ratio_PCA.pdf", width = 5, height = 3)

### 绘制火山图 ratio vocanol DEP
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
data <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)
View(data)
data$sig[data$pvalue > 0.05|data$log2FoldChange < 0.23|data$log2FoldChange > -0.23] <- "no"  ### 筛选并填充数据
data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= 0.23] <- "up"
data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -0.23] <- "down"
View(data)
data$sig =factor(data$sig, levels = c("down","up","no"))
down = table(data$sig)[1]
up = table(data$sig)[2]
nc = table(data$sig)[3]
p <- ggplot(data,aes(log2FoldChange,-log10(pvalue),color = sig))+
  geom_point(size = 0.5 )+
  xlim(-5,5) +
  scale_color_wsj('dem_rep') +
  labs( x= "log2FoldChange", y = "-log10(pvalue)",color="significance") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.23, 0.23), lty=4,col="grey",lwd=0.6) +
  scale_y_continuous(limits=c(0, 50)) +
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  annotate("text",x=-1.5,y=45,label=down)+
  annotate("text",x=1.5,y=45,label=up)+
  annotate("text",x=0,y=45,label=nc)+
  theme_hc()
p
ggsave("volcanol_DEP_ratio.pdf", width = 5, height = 5)


### ratio go enrichment bubble proteome
#"2021-09-04 09:58:12 CST"
#V1
rm(list=ls())
library(Rmisc)
library(clusterProfiler)
setwd("C:/Users/zhengxueao/OneDrive/Data&Fig/cpg/FigS/working_data")

df_ratio_p = read.csv("C:/Users/zhengxueao/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)%>%filter(log2FoldChange !="NA")
df_ratio_p$ID= str_extract(df_ratio_p$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= 0.23] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -0.23] <- "down"
ratio_up = df_ratio_p[df_ratio_p$sig == "up","ID"]
ratio_down = df_ratio_p[df_ratio_p$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df_ratio_p_up_res.csv")

ratio_up_res = read.csv("df_ratio_p_up_res_0.csv")
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
# ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))

p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="Enrichment"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  # scale_x_continuous(limits=c(1.5,5.5))+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ggsave("pro_ratio_up_0.pdf", width = 8.5, height = 9)
ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df_ratio_p_down_res.csv")

ratio_down_res = read.csv("df_ratio_p_down_res_0.csv")
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
# ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))
p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="Enrichment"
    #y="Pathway name",
    # title="Pathway enrichment")
  )+
  scale_x_continuous(limits=c(1.5,12))+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
ggsave("pro_ratio_down_0.pdf", width = 9, height = 9)

### DEG_ratio分析
#### "2021-09-03 11:25:09 CST"
##### V1
rm(list=ls())
library(readr)
library(DESeq2)
cts00 <- read.csv("D:/R_protemics/RNA_seq2/normalizedratio.CSV",row.names = 1)### change
View(cts00)
cts <- round(cts00)
View(cts)
mutant <- factor(rep(c('15_102','T26_4'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
View(colData)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","15_102","T26_4"))
res01
write.csv(as.data.frame(res01), file= 'res_DEGratio24to0.csv')
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("mutant"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=mutant)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_x_continuous(limits=c(-25, 25))+
  scale_y_continuous(limits=c(-25, 25))+
  coord_fixed()
ggsave("mRNA_ratio_PCA.pdf", width = 5, height = 3)

### ratio volcanol
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
data <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEGratio24to0.csv",header = T)
View(data)
data$sig[data$pvalue > 0.05|data$log2FoldChange < 0.584|data$log2FoldChange > -0.23] <- "no"  ### 筛选并填充数据
data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= 0.584] <- "up"
data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -0.584] <- "down"
View(data)
data$sig =factor(data$sig, levels = c("down","up","no"))
down = table(data$sig)[1]
up = table(data$sig)[2]
nc = table(data$sig)[3]
p <- ggplot(data,aes(log2FoldChange,-log10(pvalue),color = sig))+
  geom_point(size = 0.5 )+
  xlim(-5,5) +
  scale_color_wsj('dem_rep') +
  labs( x= "log2FoldChange", y = "-log10(pvalue)",color="significance") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.584, 0.584), lty=4,col="grey",lwd=0.6) +
  annotate("text",x=-2.5,y=205,label=down)+
  annotate("text",x=2.5,y=205,label=up)+
  annotate("text",x=0,y=205,label=nc)+
  #scale_y_continuous(limits=c(0, 50)) +
  #scale_x_continuous(limits=c(-2.5, 2.5)) +
  theme_hc()
p
ggsave("volcanol_DEG_ratio.pdf", width = 5, height = 5)


### Venn ratio
### "2021-09-03 11:45:41 CST"
# V1
rm(list=ls())
library(ggVennDiagram)
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
dfdeg <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEGratio24to0.csv",header = T)
dfdep <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)
dfdeg$ID= dfdeg$X
dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 0.584 & dfdeg$pvalue < 0.05,]
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -0.584 & dfdeg$pvalue < 0.05,]
length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)
x <- list(DG=dfdeg_down$ID,UG=dfdeg_up$ID,DP=dfdep_down$ID,UP=dfdep_up$ID)
ggVennDiagram(x)+ scale_fill_gradient(low="white",high = "plum1")+ scale_color_brewer()
ggvenn(x,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 4)
ggsave("Venn_DEGDEPratio.pdf", width = 6, height = 6)
now()


## correlation ratio
## "Wed Sep 01 16:43:35 2021"
## V1
setwd("C:/Users/31781/OneDrive/Data&Fig/Fig/working_data")
date()
rm(list=ls())
### ratio change (x= "mRNA ratio of (COM vs. cgb)", y="protein ration change (COM vs. cgb)", )
df005 = read.csv("D:/R_protemics/AAA_ogdata/res/res_DEGratio24to0.csv",header = T)
df006 = read.csv("D:/R_protemics/AAA_ogdata/res/res_DEPratio24to0.csv",header = T)
df056 = merge(df005,df006,by =  "ID")
p12 = ggplot(data=df056,
       aes(x=r.log2FoldChange, y=log2FoldChange)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm", color = "black", fill = "lightgray")+
  #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) +
    labs(
    x="mRNA ratio change of (cgb vs. COM)",
    y="protein ration change of (cgb vs. COM)",
  )+
  stat_cor(method = "pearson")+
  #scale_y_continuous(limits=c(0,6))+
  #scale_x_continuous(limits=c(0,75))+
  theme_classic2()
p12
ggsave("ratio change.pdf", width = 6.5, height = 5)

### mRNA total normal
# "2021-09-04 08:31:01 CST"
rm(list=ls())
library(readr)
library(DESeq2)
library(ggplot2)
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
cts00 <- read.csv("D:/R_protemics/RNA_seq2/gene_matrix02.csv",header = 1, row.names = 1)### change
cts <- round(cts00)
View(cts)
mutant <- factor(rep(c('cpg','COM'), each=6))
dpi <- factor(rep(c('0','24','0','24'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant,dpi)
View(colData)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant + dpi)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
write.csv(as.data.frame(normalized_counts), file= '../RNA_seq2/normalized_counts.csv')
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("mutant", "dpi"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=mutant, shape=dpi)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_x_continuous(limits=c(-40, 40))+
  scale_y_continuous(limits=c(-40, 40))+
  coord_fixed()
ggsave("mRNA_PCA_Total.pdf", width = 5, height = 3)
now()


### proteome PCA
"2021-09-04 08:41:43 CST"
rm(list=ls())
library(readr)
library(DESeq2)
library(ggplot2)
setwd("C:/Users/zhengxueao/OneDrive/Data&Fig/cpg/FigS/working_data")
cts00 <- read.csv("D:/R_protemics/protemics/PRO_normal_01.CSV",header = 1, row.names = 1)### change
cts = cts00*1000000
cts <- round(cts)

View(cts)
mutant <- factor(rep(c('cpg','COM'), each=6))
dpi <- factor(rep(c('0','24','0','24'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant,dpi)
View(colData)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant + dpi)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("mutant", "dpi"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=mutant, shape=dpi)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_x_continuous(limits=c(-15, 15))+
  scale_y_continuous(limits=c(-15, 15))+
  coord_fixed()
ggsave("pro_PCA_Total.pdf", width = 5, height = 3)
now()


### mRNA ratio go enrichment
#"2021-09-04 09:58:12 CST"
#V1
rm(list=ls())
library(Rmisc)
df_ratio_p = read.csv("C:/Users/zhengxueao/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEGratio24to0.csv",header = T)
df_ratio_p$ID= df_ratio_p$X  #### R 的 正则表达式 提取 .1  之前的
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= 0.584] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -0.584] <- "down"
table(df_ratio_p$sig)
ratio_up = df_ratio_p[df_ratio_p$sig == "up","ID"]
ratio_down = df_ratio_p[df_ratio_p$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df_ratio_mRNA_up_res.csv")

### 筛选过后从这里开始
ratio_up_res = read.csv("df_ratio_mRNA_up_res_0.csv")  ### 这里做过修改 更换成了 挑选后的 res
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))
p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="Enrichment"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  # scale_x_continuous(limits=c(1,4))+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ggsave("mRNA_ratio_up_0.pdf", width = 7.5, height = 9)

ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df_ratio_mrna_down_res.csv")

ratio_down_res = read.csv("df_ratio_mrna_down_res_0.csv") ### !!! 做过挑选后的 res
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))

p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="Enrichment"
    #y="Pathway name",
    # title="Pathway enrichment")
  )+
  #scale_x_continuous(limits=c(1.5,12))+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
ggsave("mRNA_ratio_down_0.pdf", width = 8.5, height = 9)


### vocanol of SA arabidopsis
### "2021-09-04 15:26:07 CST"
## V1
dfcount = read.csv("D:/AAA_spsm1/AAA_ogdata/AT_codon_count.csv",header = T)
dfdeg <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEGratio24to0.csv",header = T)
dfdep <- read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)

rm(list=ls())
setwd("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data")
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(Hmisc)
library(dplyr)
library(psych)
library(stringr)
library(ggpubr)
library(ggrepel)

df1 = read.csv("D:/AAA_spsm1/AAA_ogdata/At_go_matrix01.CSV",header = T)
df001 = read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEGratio24to0.csv",header = T)
df002 = read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)
df002$ID= str_extract(df002$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
df0011 = merge(df001,df1,by.x ="X", by.y = "ID")

df0 = df0011 %>% filter(stringr::str_detect(Term, 'salicylic'))
df0$sig <- "no"  ### 筛选并填充数据
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange >= 0.584] <- "up"
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange <= -0.584] <- "down"
df0$sig =factor(df0$sig, levels = c("down","up","no"))
df002 = read.delim("D:/AAA_spsm1/AAA_ogdata/gene_association_tair",header = F,sep="\t")
df002 = dplyr::distinct(df002,V3, .keep_all= TRUE)
df0 = merge(df002,df0,by.x = "V10",by.y ="X")
df0 = dplyr::distinct(df0,ID0, .keep_all= TRUE)

ggplot(df0,aes(log2FoldChange,-log10(pvalue), color = sig))+
  geom_point(size = 0.5 )+
  scale_color_wsj('dem_rep') +
  labs( x= "Log2(FoldChange)", y = "Log10(p-value)",color="significance") +
  theme_hc()+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.584, 0.584), lty=4,col="grey",lwd=0.6) +
  scale_y_continuous(limits=c(0, 75)) +
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  geom_text_repel(aes(log2FoldChange,-log10(pvalue), label = V3))
ggsave("RNA_sa_ratio_VOL.pdf", width = 8, height = 8)
ggsave("RNA_sa_ratio_VOL0.pdf", width = 14, height = 14)

dfcount = read.csv("D:/AAA_spsm1/AAA_ogdata/AT_codon_count2.csv",header = T)
df1 = merge(df0,dfcount, by.x = "V10",by.y = "ID")
library(ggpubr)
my_comparisons=list(c("down","no"))
la_p1 = 140
p1 <- ggplot(df1[df1$sig == "down"|df1$sig == "no",],aes(x=sig,y=s2codon,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "GAA&CAA&AAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p1
ggsave("mRNA_A_ending.pdf", width = 3, height = 5)

rm(list=ls())
df1 = read.csv("D:/AAA_spsm1/AAA_ogdata/At_go_matrix01.CSV",header = T)
df002 = read.csv("C:/Users/31781/OneDrive/Data&Fig/cpg/FigS/working_data/res_DEPratio24to0.csv",header = T)
df002$ID= str_extract(df002$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
df0012 = merge(df002,df1,by = "ID")
dfcount = read.csv("D:/AAA_spsm1/AAA_ogdata/AT_codon_count2.csv",header = T)
df0 = df0012 %>% filter(stringr::str_detect(Term, 'salicylic'))
df0$sig <- "no"  ### 筛选并填充数据
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange >= 0.23] <- "up"
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange <= -0.23] <- "down"
df0$sig =factor(df0$sig, levels = c("down","up","no"))
df002 = read.delim("D:/AAA_spsm1/AAA_ogdata/gene_association_tair",header = F,sep="\t")
df002 = dplyr::distinct(df002,V10, .keep_all= TRUE)
df0 = merge(df002,df0,by.x = "V10",by.y ="ID")
df0 = dplyr::distinct(df0,ID0, .keep_all= TRUE)
df1 = merge(df0,dfcount, by.x = "V10",by.y = "ID")

ggplot(df0,aes(log2FoldChange,-log10(pvalue), color = sig))+
  geom_point(size = 0.5 )+
  scale_color_wsj('dem_rep') +
  labs( x= "log2(FoldChange)", y = "Log10(pvalue)",color="significance") +
  theme_hc()+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.23, 0.23), lty=4,col="grey",lwd=0.6) +
  scale_y_continuous(limits=c(0, 30)) +
  scale_x_continuous(limits=c(-0.7, 0.7)) +
  geom_text_repel(aes(log2FoldChange,-log10(pvalue), label = V3))
ggsave("PRO_sa_ratio_VOL.pdf", width = 8, height = 8)

df0 = df0012 %>% filter(stringr::str_detect(Term, 'salicylic|immune'))
df0$sig <- "no"  ### 筛选并填充数据
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange >= 0.23] <- "up"
df0$sig[df0$pvalue <= 0.05 & df0$log2FoldChange <= -0.23] <- "down"
df0$sig =factor(df0$sig, levels = c("down","up","no"))
df002 = read.delim("D:/AAA_spsm1/AAA_ogdata/gene_association_tair",header = F,sep="\t")
df002 = dplyr::distinct(df002,V10, .keep_all= TRUE)
df0 = merge(df002,df0,by.x = "V10",by.y ="ID")
df0 = dplyr::distinct(df0,ID0, .keep_all= TRUE)
df1 = merge(df0,dfcount, by.x = "V10",by.y = "ID")

library(ggpubr)
my_comparisons=list(c("down","no"))
la_p1 = 100
p1 <- ggplot(df1[df1$sig == "down"|df1$sig == "no",],aes(x=sig,y=s2codon,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "GAA&CAA&AAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "t.test",label = "p.format",label.y = la_p1*3/4)
p1
ggsave("pro_A_ending.pdf", width = 3, height = 5)


### YEAST 2020 aggregation
setwd("C:/Users/zhengxueao/OneDrive/project/yeast_2021/working_data")
# setwd("C:/Users/31781/OneDrive/project/yeast_2021/working_data")

rm(list=ls())
library(Hmisc)
library(corrplot)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(psych)
library(stringr)
library(pheatmap)
library(Rmisc)
library(tidyverse)

df001 = read.csv("C:/Users/zhengxueao/OneDrive/og_data/Sc_amino_count.csv",header = T)
df002 = read.csv("res_ElpvsWT.csv",header = T)
df012 = merge(df001,df002,by.x = "ID", by.y=,"X")

#### Figure4A
df012 = df012%>% select_if(is.numeric)
df012 = df012[df012$pvalue <= 0.05 | df012$pvalue >= 0.95,]
list2 = rcorr(as.matrix(df012))
df01 = data.frame(list2[1])
df02 = data.frame(list2[3])
df01$ID = row.names(df01)

df03 = df01[1:20,c(23,28)]
df04 = df03 %>% arrange(-r.log2FoldChange)
list4 = df04$r.log2FoldChange
df05 = data.frame(list4[1:4],
                  list4[5:8],
                  list4[9:12],
                  list4[13:16],
                  list4[17:20]
)
df06 = apply(df05,2,function (fan) fan*(-1))
list5 = df04$ID
df005 = data.frame(list5[1:4],
                  list5[5:8],
                  list5[9:12],
                  list5[13:16],
                  list5[17:20]
)
pdf(file = "FigureA.pdf")
corrplot(as.matrix(df06),
         method = "circle",
         type = 'full',
         col = colorRampPalette(c("#003366","white" ))(8),
         #cl.lim = c(-0.16,0),
         tl.col = "black",
         addrect = 8,
         is.corr = FALSE )
dev.off()
write.csv(as.data.frame(df01),"cor_amnio_pro_yeast.csv")
write.csv(as.data.frame(df02),"cor_amnio_pro_yeast_p.csv")

###Figure4D
df034 = df012
df034$KRED = df034$K +df034$R+df034$E+df034$D
df034$KR = df034$K +df034$R
df034$ED = df034$E +df034$D
df034$FYW = df034$F + df034$Y + df034$W
df034$STNQ = df034$S +df034$T+ df034$N +df034$Q
df034$CGP = df034$C + df034$G + df034$P
df034$AVILM = df034$A + df034$V +df034$I +df034$L +df034$M
df034$STNQCH = df034$S +df034$T+ df034$N +df034$Q +df034$C +df034$H
df034$FAVILMGPYW = df034$F + df034$A + df034$V +df034$I +df034$L +df034$M + df034$G + df034$P + df034$Y + df034$W
df034$sig <- "NC"
df034$sig[df034$pvalue <= 0.05 & df034$log2FoldChange >= 0.23] <- "UP"
df034$sig[df034$pvalue <= 0.05 & df034$log2FoldChange <= -0.23] <- "DOWN"
df034$sig =factor(df034$sig, levels = c("DOWN","UP","NC"))
table(df034$sig)
write.csv(as.data.frame(df034),"amino_yeast.csv")
table(df034$sig)
df034$Total_scale[df034$Total > 0 & df034$Total < 100] <- "0-100"  ### 筛选并填充数据
df034$Total_scale[df034$Total >= 100 & df034$Total < 200] <- "100-200"  ### 筛选并填充数据
df034$Total_scale[df034$Total >= 200 & df034$Total < 300] <- "200-300"
df034$Total_scale[df034$Total >= 300 & df034$Total < 400] <- "300-400"
df034$Total_scale[df034$Total >= 400 & df034$Total < 500] <- "400-500"
df034$Total_scale[df034$Total >= 500 & df034$Total < 600] <- "500-600"
df034$Total_scale[df034$Total >= 600 & df034$Total < 700] <- "600-700"
df034$Total_scale[df034$Total >= 700 & df034$Total < 800] <- "700-800"
df034$Total_scale[df034$Total >= 800 & df034$Total < 900] <- "800-900"
df034$Total_scale[df034$Total >= 900 & df034$Total < 1000] <- "900-1000"
df034$Total_scale[df034$Total >= 1000 & df034$Total < 1100] <- "1000-1100"
df034$Total_scale[df034$Total >= 1100 ] <- ">1100"
df034$Total_scale =factor(df034$Total_scale, levels = c("0-100","100-200","200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000","1000-1100",">1100"))
View(df034)
df_scale = df034[,c("sig","Total_scale","Total","log2FoldChange")]
df_scale$sig =factor(df_scale$sig, levels = c("DOWN","UP","NC"))
df_scale <- subset(df_scale,df_scale$sig!="NC")
df_scale$Total_scale =factor(df_scale$Total_scale, levels = c("0-100","100-200","200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000","1000-1100",">1100"))

p <- ggplot(df_scale,aes(x=Total_scale,fill=sig)) +
  geom_bar(position="fill",width=0.6)+
  labs( x= "", y = "Percentage",fill="DEP") +
  scale_fill_wsj('dem_rep') +
  theme_classic2()+
  theme(axis.text.x = element_text(size = 15, vjust=0.6,angle = 45))+
  theme(axis.text.y = element_text(size = 15))
p
ggsave("FigureD.pdf", width = 10, height = 4)

###Figure4E
p <- ggplot(df_scale,aes(x=Total_scale,y=log2FoldChange,fill=Total_scale)) +
  geom_boxplot(position="dodge",outlier.colour = NA)+
  labs( x= "", y = "Log2FC",fill="DEP") +
  scale_fill_brewer(palette = "Spectral")+
  scale_y_continuous(limits=c(-6, 6))+
  guides(fill=FALSE)+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 15, vjust=0.6,angle = 45))+
  theme(axis.text.y = element_text(size = 15))
p
ggsave("FigureE.pdf", width = 10, height = 4)

### Figure4B
df012 = df034
my_comparisons=list(c("DOWN","UP"))
la_p3 = 1300
p3 <- ggplot(df012[df012$sig == "DOWN"|df012$sig == "UP",],aes(x=sig,y=Total,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "Protein length",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p3))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p3*2/3)
p3
ggsave("FigureB.pdf", width = 3, height = 5)

### Figure4C
x = mean(c(df012[df012$sig == "DOWN","Total"]))-mean(c(df012[df012$sig == "UP","Total"]))
table(df012$sig)
N1 = sum(df012$sig == "DOWN")
N2 = sum(df012$sig == "UP")
table(df012$sig)
df012 = subset(df012,df012$Total!="Inf")
sample_dovsup <- function(){
c1 = mean(c(sample(df012$Total,N1,replace = FALSE))) - mean(c(sample(df012$Total,N2,replace = FALSE)))
return(c1)
}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p14 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency") +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-280,280), breaks = c(-280,-240,-200,-160,-120,-80,-40,0,40,80,120,160,200,240,280))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p14
t.test(df001$c2,mu=x)
mean(c(df012[df012$sig == "DOWN","Total"]))-mean(c(df012[df012$sig == "UP","Total"]))
ggsave("FigureC.pdf", width = 6, height = 6)


### median  154aa
x = median(c(df012[df012$sig == "DOWN","Total"]))-median(c(df012[df012$sig == "UP","Total"]))
table(df012$sig)
N1 = sum(df012$sig == "DOWN")
N2 = sum(df012$sig == "UP")
table(df012$sig)
df012 = subset(df012,df012$Total!="Inf")
sample_dovsup <- function(){
c1 = median(c(sample(df012$Total,N1,replace = FALSE))) - median(c(sample(df012$Total,N2,replace = FALSE)))
return(c1)
}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p14 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency") +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-280,280), breaks = c(-280,-240,-200,-160,-120,-80,-40,0,40,80,120,160,200,240,280))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p14
t.test(df001$c2,mu=x)
median(c(df012[df012$sig == "DOWN","Total"]))-median(c(df012[df012$sig == "UP","Total"]))
ggsave("FigureC.pdf", width = 6, height = 6)




rm(list=ls())
# setwd("C:/Users/31781/OneDrive/Data&Fig/Fig2/working_data")
# library(ggVennDiagram)

setwd("C:/Users/zhengxueao/OneDrive/og_data/AT")
library(ggplot2)

df001 = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEP0H.csv",header = T)
df002 = read.csv("D:/AAA_spsm1/AAA_ogdata/AT_codon_count.csv",header = T)
df003 = read.csv("D:/AAA_spsm1/AAA_ogdata/AT_amino_count.csv",header = T)
df023 = merge(df002,df003,by.x = "X", by.y=,"X")

df023_0 = dplyr::distinct(df023,X,.keep_all= TRUE)
df023_0$ID= str_extract(df023_0$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
df023_0= arrange(df023_0,-df023_0$Total)
df023_0 = dplyr::distinct(df023_0,ID,.keep_all= TRUE)

df012 = merge(df001,df023_0,by.x = "ID", by.y=,"ID")
dfdeg <- read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEG0H.csv")### change
dfdep <- df012### change
dfdeg = na.omit(dfdeg)

df0_p = dfdep
df0_p$sig <- "no"  ### 筛选并填充数据
df0_p$sig[df0_p$pvalue <= 0.05 & df0_p$log2FoldChange >= 0.23] <- "up"
df0_p$sig[df0_p$pvalue <= 0.05 & df0_p$log2FoldChange <=  -0.23] <- "down"
table(df0_p$sig)

dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 0.584 & dfdeg$pvalue < 0.05,]
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -0.584 & dfdeg$pvalue < 0.05,]
length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)

x <- list(DG=dfdeg_down$ID,UG=dfdeg_up$ID,DP=dfdep_down$ID,UP=dfdep_up$ID)
ggVennDiagram(x)
ggsave("Venn_DEGDEP0H.pdf", width = 6, height = 6)

### Figure 2
setwd("C:/Users/31781/OneDrive/Data&Fig/Fig2/working_data")
rm(list=ls())
library(clusterProfiler)
library(ggplot2)
library(ggthemes)
library(org.At.tair.db)
library(Rmisc)
library(stringr)
library(dplyr)

##### df0 不分差异基因
rm(list=ls())
df = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEG0H.csv",header = T)
df$sig <- "no"  ### 筛选并填充数据
df$sig[df$pvalue <= 0.05 & df$log2FoldChange >= 0.584] <- "sig"
df$sig[df$pvalue <= 0.05 & df$log2FoldChange <=  -0.584] <- "sig"
sig_id = df[df$sig == "sig",]$ID
sig_res = enrichGO(sig_id, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(sig_res),"df_R0H_sig_res.csv")

emapplot(sig_res,
         showCategory = 500,
         pie="count",
         pie_scale=2)

ggsave("AT_ROH_BP.pdf", width = 40, height = 40,limitsize = FALSE)
ggsave("AT_ROH_BP144.pdf", width = 20, height = 20,limitsize = FALSE)

df = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEP0H.csv",header = T)
df$sig <- "no"  ### 筛选并填充数据
df$sig[df$pvalue <= 0.05 & df$log2FoldChange >= 0.23] <- "sig"
df$sig[df$pvalue <= 0.05 & df$log2FoldChange <=  -0.23] <- "sig"
sig_id = df[df$sig == "sig",]$ID
sig_res = enrichGO(sig_id, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(sig_res),"df_P0H_sig_res.csv")
emapplot(sig_res,
         showCategory = 500,
         pie="count",
         pie_scale=2)

ggsave("AT_POH_BP.pdf", width = 40, height = 40,limitsize = FALSE)
ggsave("AT_POH_BP144.pdf", width = 20, height = 20,limitsize = FALSE)

### df0_ratio
rm(list=ls())
df_ratio_p = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEPratio24to0.csv",header = T)
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= 0.23] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -0.23] <- "down"
ratio_up = df_ratio_p[df_ratio_p$sig == "up","ID"]
ratio_down = df_ratio_p[df_ratio_p$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")

write.csv(as.data.frame(ratio_up_res),"df_ratio_p_up_res.csv")
ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))

p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1

ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df_ratio_p_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))

p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
p=list(p1,p2)
pdf(file = "Figure_df_ratio_p.pdf",18,9)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()

##### df0_r
rm(list=ls())
df0_r = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEG0H.csv",header = T)
df0_r$sig <- "no"  ### 筛选并填充数据
df0_r$sig[df0_r$pvalue <= 0.01 & df0_r$log2FoldChange >= 0.584] <- "up"
df0_r$sig[df0_r$pvalue <= 0.01 & df0_r$log2FoldChange <=  -0.584] <- "down"
table(df0_r$sig)
ratio_up = df0_r[df0_r$sig == "up","ID"]
ratio_down = df0_r[df0_r$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df0_r_up_res.csv")
ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))
p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df0_r_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))
p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
p=list(p1,p2)
pdf(file = "Figure_df0_r.pdf",18,9)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()

##### df24_r
rm(list=ls())
df24_r = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEG24H.csv",header = T)
df24_r$sig <- "no"  ### 筛选并填充数据
df24_r$sig[df24_r$pvalue <= 0.05 & df24_r$log2FoldChange >= 0.584] <- "up"
df24_r$sig[df24_r$pvalue <= 0.05 & df24_r$log2FoldChange <=  -0.584] <- "down"
ratio_up = df24_r[df24_r$sig == "up","ID"]
ratio_down = df24_r[df24_r$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df24_r_up_res.csv")
ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))
p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df24_r_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))
p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
p=list(p1,p2)
pdf(file = "Figure_df24_r.pdf",20,12)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()

##### df0_p
rm(list=ls())
df0_p = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEP0H.csv",header = T)
df0_p$sig <- "no"  ### 筛选并填充数据
df0_p$sig[df0_p$pvalue <= 0.05 & df0_p$log2FoldChange >= 0.263] <- "up"
df0_p$sig[df0_p$pvalue <= 0.05 & df0_p$log2FoldChange <=  -0.263] <- "down"
table(df0_p$sig)
ratio_up = df0_p[df0_p$sig == "up","ID"]
ratio_down = df0_p[df0_p$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df0_p_up_res.csv")
ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))
p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df0_p_down_res.csv")
ratio_down_res = read.csv("C:/Users/31781/OneDrive/Data&Fig/Fig2/working_data/df0_p_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))
p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
ggsave("df0_p_down.pdf", width = 12, height = 9)

p=list(p1,p2)
pdf(file = "Figure_df0_p.pdf",18,9)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()

##### df24_p
rm(list=ls())
df24_p = read.csv("D:/AAA_spsm1/AAA_ogdata/res/res_DEP24H.csv",header = T)
df24_p$sig <- "no"  ### 筛选并填充数据
df24_p$sig[df24_p$pvalue <= 0.05 & df24_p$log2FoldChange >= 0.263] <- "up"
df24_p$sig[df24_p$pvalue <= 0.05 & df24_p$log2FoldChange <=  -0.263] <- "down"
ratio_up = df24_p[df24_p$sig == "up","ID"]
ratio_down = df24_p[df24_p$sig == "down","ID"]
ratio_up_res =  enrichGO(ratio_up, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_up_res),"df24_p_up_res.csv")
ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))

p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p1
ratio_down_res =  enrichGO(ratio_down, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
write.csv(as.data.frame(ratio_down_res),"df24_p_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))
p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x=""
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
p=list(p1,p2)
pdf(file = "Figure_df24_p.pdf",18,12)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()


### chloroplast
rm(list=ls())
library(pheatmap)

setwd("D:/project/chloroplast/working_data")
df0 <- read.csv("TOC_family_RNA1.CSV",header = T)

pheatmap(df0[,2:7],
         scale = "column",
         border=FALSE,
         cellwidth = 10,
         cellheight = 20,
         main = "RNA expression",
         color = colorRampPalette(c("#003366","white","#990000" ))(50), #### B2182B
         cutree_cols = 1,
         filename = "RNA.pdf",
         labels_row = df0$name)

df1 <- read.csv("TOC_famlily_PRO1.CSV") ### change
pheatmap(df1[,1:6],
         scale = "column",
         border=FALSE,
         cellwidth = 10,
         cellheight = 20,
         main = "protein expression",
         color = colorRampPalette(c("#003366","white","#990000" ))(50), #### B2182B
         cutree_cols = 1,
         filename = "protein.pdf",
         labels_row = df1$name)


### 完成热图 及 GO
### 差异表达的相关性
### 差异表达的S2 codon

#### HEATMAP
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)

df_dep15_list = df_dep_15 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_dep26_list = df_dep_26 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_deg_15_list = df_deg_15%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)
df_deg_26_list = df_deg_26%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)

dep_df = merge(df_dep26_list,df_dep15_list,all = TRUE)%>%
  left_join(df_pro[,2:14],by="ID")%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,1]

deg_df = merge(df_deg_15_list,df_deg_26_list,all = TRUE)%>%
  left_join(df_rna[,2:14],by="ID")%>%
  distinct(ID,.keep_all = T)
row.names(deg_df) = deg_df[,1]

dep_matrix = dep_df%>%
  dplyr::select(c(2:13))%>%
  apply(1,scale)
row.names(dep_matrix) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                                                "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")

deg_matrix = deg_df%>%
  dplyr::select(c(2:13))%>%
  apply(1,scale)
row.names(deg_matrix) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                                                "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")

dep_m_4=t(dep_matrix)%>%
  as.data.frame()%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(COM_0,COM_24,cgb_0,cgb_24)%>%
  as.matrix()
write.table(dep_m_4,file = "heatmatrix_dep.txt",row.names = T,col.names = T,quote = F, sep = "\t")

deg_ratio_m_4=t(deg_matrix)%>%
  as.data.frame()%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_0/cgb_24)%>%
  mutate(COM_ratio=COM_0/COM_24)%>%
  dplyr::select(cgb_ratio,COM_ratio)%>%as.matrix()
write.table(deg_ratio_m_4,file = "heatmatrix_ratio_deg.txt",row.names = T,col.names = T,quote = F, sep = "\t")

dep_ratio_m_4=t(dep_matrix)%>%
  as.data.frame()%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_0/cgb_24)%>%
  mutate(COM_ratio=COM_0/COM_24)%>%
  dplyr::select(cgb_ratio,COM_ratio)%>%as.matrix()
write.table(dep_ratio_m_4,file = "heatmatrix_ratio_dep.txt",row.names = T,col.names = T,quote = F, sep = "\t")


### dep/deg
library(clusterProfiler)
library(RColorBrewer)
library(viridis)
library(org.At.tair.db)

pheatmap(deg_m_4,
         cluster_cols = F,
         cluster_rows = TRUE,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = colorRampPalette(c("#003366","white","red"))(60), #### B2182B
         cutree_cols = 2,
         cutree_rows =3,
         filename = "H2.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 5,
         height = 15)

list=pheatmap(deg_m_4,
         cluster_cols = F,
         cluster_rows = TRUE,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 2,
         cutree_rows =2,
         filename = "H3.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 3,
         height = 15)

####
row_cluster=cutree(list$tree_row,k=2)
ID = row.names(as.data.frame(row_cluster))

cluster_go_1=as.data.frame(row_cluster)%>%
  mutate(ID=ID)%>%
  filter(row_cluster == 1)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_1),"cluster_go_1.csv")

library(ggplot2)
df1 = as.data.frame(cluster_go_1)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust1.pdf",width = 6,height =6)

cluster_go_2=as.data.frame(row_cluster)%>%
  mutate(ID=ID)%>%
  filter(row_cluster == 2)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_2),"cluster_go_2.csv")

df1 = as.data.frame(cluster_go_2)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust2.pdf",width = 6.5,height =9)


### ratio 数据分析
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 找到所有 DEP/DEG
df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)

df_dep15_list = df_dep_15 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_dep26_list = df_dep_26 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_deg_15_list = df_deg_15%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)
df_deg_26_list = df_deg_26%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)

### 合并所有 DEG/DEP
dep_df = merge(df_dep26_list,df_dep15_list,all = TRUE)%>%
  left_join(df_pro[,2:14],by="ID")%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,1]

deg_df = merge(df_deg_15_list,df_deg_26_list,all = TRUE)%>%
  left_join(df_rna[,2:14],by="ID")%>%
  distinct(ID,.keep_all = T)
row.names(deg_df) = deg_df[,1]

### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("ID","cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")
colnames(deg_df) = c("ID","cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")

dep_ratio_m = dep_df%>%
  dplyr::select(c(2:13))%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_24/cgb_0)%>%
  mutate(COM_ratio=COM_24/COM_0)%>%
  select(cgb_ratio,COM_ratio)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  mutate(order=cgb_ratio-COM_ratio)%>%
  mutate(mean0=(cgb_ratio+COM_ratio)/2)%>%
  mutate(cgb_ratio=cgb_ratio/mean0)%>%
  mutate(COM_ratio=COM_ratio/mean0)%>%
  arrange(order)
write.table(dep_ratio_m,file = "heatmatrix_ratio_dep.txt",row.names = T,col.names = T,quote = F, sep = "\t")

deg_ratio_m = deg_df%>%
  dplyr::select(c(2:13))%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_24/cgb_0)%>%
  mutate(COM_ratio=COM_24/COM_0)%>%
  select(cgb_ratio,COM_ratio)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  mutate(order=cgb_ratio-COM_ratio)%>%
  mutate(mean0=(cgb_ratio+COM_ratio)/2)%>%
  mutate(cgb_ratio=cgb_ratio/mean0)%>%
  mutate(COM_ratio=COM_ratio/mean0)%>%
  arrange(order)
write.table(deg_ratio_m,file = "heatmatrix_ratio_deg.txt",row.names = T,col.names = T,quote = F, sep = "\t")

### DEP ratio
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
library(clusterProfiler)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(org.At.tair.db)

pheatmap(dep_ratio_m[,c(2,1)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = colorRampPalette(c("#003366","white","red"))(60), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H2.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)
pheatmap(dep_ratio_m[,c(2,1)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H3.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)

cluster_go_1=dep_ratio_m %>%
  mutate(ID=rownames(dep_ratio_m))%>%
  filter(order>0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.53,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_1),"cluster_go_1.csv")

library(ggplot2)
df1 = as.data.frame(cluster_go_1)[1:20,]%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust1.pdf",width = 6,height =8)

cluster_go_2=dep_ratio_m %>%
  mutate(ID=rownames(dep_ratio_m))%>%
  filter(order<0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_2),"cluster_go_2.csv")

df1 = as.data.frame(cluster_go_2)[1:20,]%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust2.pdf",width = 6,height =8)


#### DEG ratio
pheatmap(deg_ratio_m[,c(2,1)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = colorRampPalette(c("#003366","white","red"))(60), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H2.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)

pheatmap(deg_ratio_m[,c(2,1)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H3.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)

cluster_go_1=deg_ratio_m %>%
  mutate(ID=rownames(deg_ratio_m))%>%
  filter(order>0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.53,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_1),"cluster_go_1.csv")

library(ggplot2)
df1 = as.data.frame(cluster_go_1)[1:20,]%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust1.pdf",width = 6,height =8)

cluster_go_2=deg_ratio_m %>%
  mutate(ID=rownames(deg_ratio_m))%>%
  filter(order<0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_2),"cluster_go_2.csv")

df1 = as.data.frame(cluster_go_2)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust2.pdf",width = 6,height =8)


### pro/mRNA数据分析
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 找到所有 DEP/DEG
df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)
df_dep15_list = df_dep_15 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_dep26_list = df_dep_26 %>%
  filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
  select(ID)
df_deg_15_list = df_deg_15%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)
df_deg_26_list = df_deg_26%>%
  filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
  select(ID)

### 合并所有 DEG/DEP
dep_df = merge(df_dep26_list,df_dep15_list,all = TRUE)%>%
  left_join(df_pro[,2:14],by="ID")%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,1]

deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("ID","cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

dep_ratio_m = dep_df%>%
  dplyr::select(c(2:13))%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_24/cgb_0)%>%
  mutate(COM_ratio=COM_24/COM_0)%>%
  select(cgb_ratio,COM_ratio)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  # mutate(order=cgb_ratio-COM_ratio)%>%
  # mutate(mean0=(cgb_ratio+COM_ratio)/2)%>%
  # mutate(cgb_ratio=cgb_ratio/mean0)%>%
  # mutate(COM_ratio=COM_ratio/mean0)%>%
  # arrange(order)%>%
  mutate(pro_cgb_ratio = cgb_ratio,
         pro_COM_ratio = COM_ratio)%>%
  mutate(ID = rownames(dep_ratio_m)) %>%
  select(pro_COM_ratio,pro_cgb_ratio,ID)

deg_ratio_m = deg_df%>%
  dplyr::select(c(1:12))%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_24/cgb_0)%>%
  mutate(COM_ratio=COM_24/COM_0)%>%
  select(cgb_ratio,COM_ratio)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  mutate(rna_cgb_ratio = cgb_ratio, rna_COM_ratio = COM_ratio)%>%
  mutate(ID = rownames(deg_ratio_m)) %>%
  select(rna_COM_ratio,rna_cgb_ratio,ID)%>%
  right_join(dep_ratio_m, by="ID")

dep_deg_te = deg_ratio_m%>%
  right_join(dep_ratio_m, by="ID")%>%
  mutate(TE_COM = pro_COM_ratio/rna_COM_ratio)%>%
  mutate(TE_cgb = pro_cgb_ratio/rna_cgb_ratio)%>%
  mutate(order=TE_COM-TE_cgb)%>%
  mutate(mean0=(TE_COM+TE_cgb)/2)%>%
  mutate(TE_cgb=TE_cgb/mean0)%>%
  mutate(TE_COM=TE_COM/mean0)%>%
  arrange(-order)%>%
  dplyr::select(TE_COM,TE_cgb,order,ID)%>%
  filter(TE_COM != "NA")

### DEP ratio
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
library(clusterProfiler)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(org.At.tair.db)

pheatmap(dep_deg_te[,c(1,2)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = colorRampPalette(c("#003366","white","red"))(60), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H2.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)
pheatmap(dep_deg_te[,c(1,2)],
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H3.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)

cluster_go_1=dep_deg_te %>%
  mutate(order=TE_COM-TE_cgb) %>%
  filter(order>0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.53,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_1),"cluster_go_1.csv")

library(ggplot2)
df1 = as.data.frame(cluster_go_1)[1:20,]%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust1.pdf",width = 6,height =8)

cluster_go_2=dep_deg_te %>%
  mutate(order=TE_COM-TE_cgb) %>%
  filter(order<0)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_2),"cluster_go_2.csv")

df1 = as.data.frame(cluster_go_2)[1:20,]%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust2.pdf",width = 8,height =8)




### TE degp analysis 数据log
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

dep_ratio_m = dep_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,COM_0=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  mutate(cgb_24V0_1=log(cgb_24_1/cgb_0),
              cgb_24V0_2=log(cgb_24_2/cgb_0),
              cgb_24V0_3=log(cgb_24_3/cgb_0),
              COM_24V0_1=log(COM_24_1/COM_0),
              COM_24V0_2=log(COM_24_2/COM_0),
              COM_24V0_3=log(COM_24_3/COM_0))%>%
  dplyr::select(ID,cgb_24V0_1,cgb_24V0_2,cgb_24V0_3,COM_24V0_1,COM_24V0_2,COM_24V0_3)

dep_deg_te = deg_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=log(cgb_24/cgb_0))%>%
  mutate(COM_ratio=log(COM_24/COM_0))%>%
  dplyr::select(cgb_ratio,COM_ratio,ID)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  dplyr::mutate(rna_cgb_ratio = cgb_ratio, rna_COM_ratio = COM_ratio)%>%
  dplyr::select(rna_COM_ratio,rna_cgb_ratio,ID)%>%
  right_join(dep_ratio_m, by="ID")%>%
  mutate(te_cgb_1=2^(cgb_24V0_1/rna_cgb_ratio)*100000,
              te_cgb_2=2^(cgb_24V0_2/rna_cgb_ratio)*100000,
              te_cgb_3=2^(cgb_24V0_3/rna_cgb_ratio)*100000,
              te_COM_1=2^(COM_24V0_1/rna_COM_ratio)*100000,
              te_COM_2=2^(COM_24V0_2/rna_COM_ratio)*100000,
              te_COM_3=2^(COM_24V0_3/rna_COM_ratio)*100000)%>%
  dplyr::select(ID,te_cgb_1,te_cgb_2,te_cgb_3,
         te_COM_1,te_COM_2,te_COM_3)%>%
  filter(te_cgb_1 != "NA",te_cgb_1 != Inf,te_COM_3!= Inf)
rownames(dep_deg_te)=dep_deg_te$ID
write.csv(as.data.frame(dep_deg_te), file='TE_normalized.csv')

library(DESeq2)
cts <- round(dep_deg_te[,2:7])
mutant <- factor(rep(c('cgb','COM'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","cgb","COM"))
write.csv(as.data.frame(res01), file='res_TE_cgbvsCOM.csv')

### 火山图作图
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
data <- read.csv("res_TE_cgbvsCOM.csv",header = T)
View(data)
data$sig<- "no"  ### 筛选并填充数据
data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= 0.584] <- "up"
data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -0.584] <- "down"
View(data)
data$sig =factor(data$sig, levels = c("down","up","no"))
down = table(data$sig)[1]
up = table(data$sig)[2]
nc = table(data$sig)[3]
p <- ggplot(data,aes(log2FoldChange,-log10(pvalue),color = sig))+
  geom_point(size = 0.5 )+
  xlim(-5,5) +
  scale_color_wsj('dem_rep') +
  labs( x= "log2FoldChange", y = "-log10(pvalue)",color="significance") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.23, 0.23), lty=4,col="grey",lwd=0.6) +
  scale_y_continuous(limits=c(0, 50)) +
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  annotate("text",x=-1.5,y=45,label=down)+
  annotate("text",x=1.5,y=45,label=up)+
  annotate("text",x=0,y=45,label=nc)+
  theme_hc()
p
ggsave("volcanol_TE_ratio.pdf", width = 5, height = 5)

### TE down GO 图形生成
rm(list=ls())
setwd("D:/BigData/cgb")
library(clusterProfiler)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)
library(Rmisc)
df_ratio_p = read.csv("res_TE_cgbvsCOM.csv",header = T)
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= 0.23] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -0.23] <- "down"

df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df_ratio_p%>%filter(sig == "down")%>%
  left_join(df_anno,by=c("X"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("TE_down_TableS7.csv",row.names = F)

cgb_UP_go=df_ratio_p%>%filter(sig == "down")%>%pull(X)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.50,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)
write.csv(cgb_UP_go,"TE_down_tableS7.csv")

df1 = as.data.frame(cgb_UP_go)%>%
  arrange(pvalue)%>%dplyr::slice(1:25)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "less upregulated"

ggplot(df1,aes(x=-log(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size=""
      # x="Enrichment"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("TE_DOWN.pdf",width = 5.5,height =6)



ratio_up_res =  df_ratio_p%>%
  filter(sig == "up")%>%
  pull(X)%>%
  enrichGO( OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%simplify(cutoff = 0.55,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)
write.csv(as.data.frame(ratio_up_res),"df_ratio_p_up_res.csv")

ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:30,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))

p1 = ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="",
    # y="Pathway name",
    title="cgb vs COM upregulated")+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_text(size=rel(1.3))
  )
p1

ratio_down_res =  df_ratio_p%>%
  filter(sig == "down")%>%
  pull(X)%>%
  enrichGO( OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%simplify(cutoff = 0.56,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

write.csv(as.data.frame(ratio_down_res),"df_ratio_p_down_res.csv")
ratio_down_res = as.data.frame(ratio_down_res)
ratio_down_res$B1 = as.numeric(str_extract(ratio_down_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$B2 = as.numeric(str_extract(ratio_down_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$BR = ratio_down_res$B1/ratio_down_res$B2
ratio_down_res$G1 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$G2 = as.numeric(str_extract(ratio_down_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_down_res$GR = ratio_down_res$G1/ratio_down_res$G2
ratio_down_res$EnrichmentFold = ratio_down_res$GR/ratio_down_res$BR
ratio_down_res = arrange(ratio_down_res,-ratio_down_res$Count)
ratio_down_res = ratio_down_res[1:30,]
ratio_down_res = arrange(ratio_down_res,ratio_down_res$EnrichmentFold)
ratio_down_res$Description = factor(ratio_down_res$Description,levels = c(ratio_down_res$Description))

p2 = ggplot(ratio_down_res,aes(x=EnrichmentFold,y=Description)) +
  geom_point(aes(size=Count,color=pvalue))+
  scale_colour_gradient(low="red",high="blue",n.breaks = 7)+
  labs(
    size="",
    x="",
    title="cgb vs COM downregulated"
  )+
  theme_gdocs()+
  theme(
    axis.text.y = element_text(size = rel(1.3)),
    axis.title.x = element_text(size=rel(1.3)),
    axis.title.y = element_blank()
  )
p2
p=list(p2,p1)
pdf(file = "Figure_df_ratio_p.pdf",18,9)
multiplot(plotlist = p[1:2], cols = 2)
dev.off()

### arabidopsis 注释模块
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df_te =  df_ratio_p%>%
  left_join(df_anno,by=c("X"="TAIR"))%>%
  write.csv("df_te_anno.csv",row.names = F)

### 从文字中提取list
x =ratio_down_res[5,"geneID"]
xl= as.data.frame(unlist(strsplit(x,"/")),mode="list")%>%
  mutate(TAIR = unlist(strsplit(x,"/")))%>%
  left_join(df_anno,by="TAIR")%>%
  dplyr::select(-1)%>%
  write.csv("df_bacter_te_anno.csv",row.names = F)

### vven 复杂Venn 图
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "-",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

# dfdeg=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEGratio24to0.csv",header = T)
# dfdep=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEPratio24to0.csv",header = T)
df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)
dfdeg = df_deg_15
dfdep = df_dep_15
# dfdeg$ID= dfdeg$X
# dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 1 & dfdeg$pvalue < 0.05,]
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -1 & dfdeg$pvalue < 0.05,]
length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)
cgb_UG = dfdeg_up
cgb_DG = dfdeg_down
cgb_UP = dfdep_up
cgb_DP = dfdep_down
dfdeg = df_deg_26
dfdep = df_dep_26
# dfdeg$ID= dfdeg$X
# dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 1 & dfdeg$pvalue < 0.05,]
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -1 & dfdeg$pvalue < 0.05,]
length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)
COM_UG = dfdeg_up
COM_DG = dfdeg_down
COM_UP = dfdep_up
COM_DP = dfdep_down
library(ggVennDiagram)
library(ggvenn)
X_COM_up = list(COM_UP=COM_UP$ID,COM_UG = COM_UG$ID)
X_cgb_up = list(cgb_UP=cgb_UP$ID,cgb_UG = cgb_UG$ID)
X_COM_down = list(COM_DP=COM_DP$ID,COM_DG = COM_DG$ID)
X_cgb_down = list(cgb_DP=cgb_DP$ID,cgb_DG = cgb_DG$ID)

p1=ggVennDiagram(X_COM_up)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p2=ggVennDiagram(X_cgb_up)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p3=ggVennDiagram(X_COM_down)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p4=ggVennDiagram(X_cgb_down)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

X_COM_cgb_up_pro = list(COM_UP=COM_UP$ID,cgb_UP = cgb_UP$ID)
X_COM_cgb_up_rna = list(COM_UG=COM_UG$ID,cgb_UG = cgb_UG$ID)

p1=ggVennDiagram(X_COM_cgb_up_pro)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p2=ggVennDiagram(X_COM_cgb_up_rna)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

library(cowplot)
pdf(file = "ggvenn.pdf",20,17)
plot_grid(p1,p2,nrow = 2)
dev.off()

### 这个部分数据要仔细想好
### TE degp analysis
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

# ### 找到所有 DEP/DEG
# df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
# df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
# df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
# df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)
# df_dep15_list = df_dep_15 %>%
#   filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
#   select(ID)
# df_dep26_list = df_dep_26 %>%
#   filter(pvalue<0.05,abs(log2FoldChange)>0.263)%>%
#   select(ID)
# df_deg_15_list = df_deg_15%>%
#   filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
#   select(ID)
# df_deg_26_list = df_deg_26%>%
#   filter(pvalue<0.05,abs(log2FoldChange)>1)%>%
#   select(ID)

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
dep_ratio_m = dep_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,COM_0=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  mutate(cgb_24V0_1=100000*cgb_24_1/cgb_0,
              cgb_24V0_2=100000*cgb_24_2/cgb_0,
              cgb_24V0_3=100000*cgb_24_3/cgb_0,
              COM_24V0_1=100000*COM_24_1/COM_0,
              COM_24V0_2=100000*COM_24_2/COM_0,
              COM_24V0_3=100000*COM_24_3/COM_0)%>%
  dplyr::select(ID,cgb_24V0_1,cgb_24V0_2,cgb_24V0_3,COM_24V0_1,COM_24V0_2,COM_24V0_3)
dep_deg_te = deg_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  mutate(cgb_ratio=cgb_24/cgb_0)%>%
  mutate(COM_ratio=COM_24/COM_0)%>%
  dplyr::select(cgb_ratio,COM_ratio,ID)%>%
  filter(cgb_ratio != 0,cgb_ratio != Inf)%>%
  filter(COM_ratio != 0,COM_ratio != Inf)%>%
  mutate(rna_cgb_ratio = cgb_ratio, rna_COM_ratio = COM_ratio)%>%
  dplyr::select(rna_COM_ratio,rna_cgb_ratio,ID)%>%
  right_join(dep_ratio_m, by="ID")%>%
  mutate(te_cgb_1=cgb_24V0_1/rna_cgb_ratio,
              te_cgb_2=cgb_24V0_2/rna_cgb_ratio,
              te_cgb_3=cgb_24V0_3/rna_cgb_ratio,
              te_COM_1=COM_24V0_1/rna_COM_ratio,
              te_COM_2=COM_24V0_2/rna_COM_ratio,
              te_COM_3=COM_24V0_3/rna_COM_ratio)%>%
  dplyr::select(ID,te_cgb_1,te_cgb_2,te_cgb_3,
         te_COM_1,te_COM_2,te_COM_3)%>%
  filter(te_cgb_1 != "NA")
rownames(dep_deg_te)=dep_deg_te$ID
write.csv(as.data.frame(dep_deg_te), file='TE_normalized.csv')

library(DESeq2)
cts <- round(dep_deg_te[,2:7])
mutant <- factor(rep(c('cgb','COM'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","cgb","COM"))
write.csv(as.data.frame(res01), file='res_TE_cgbvsCOM.csv')

library(ggplot2)
library(RColorBrewer)
library(ggthemes)
data <- read.csv("res_TE_cgbvsCOM.csv",header = T)
View(data)
data$sig<- "no"  ### 筛选并填充数据
data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= 0.23] <- "up"
data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -0.23] <- "down"
View(data)
data$sig =factor(data$sig, levels = c("down","up","no"))
down = table(data$sig)[1]
up = table(data$sig)[2]
nc = table(data$sig)[3]
p <- ggplot(data,aes(log2FoldChange,-log10(pvalue),color = sig))+
  geom_point(size = 0.5 )+
  xlim(-5,5) +
  scale_color_wsj('dem_rep') +
  labs( x= "log2FoldChange", y = "-log10(pvalue)",color="significance") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.23, 0.23), lty=4,col="grey",lwd=0.6) +
  scale_y_continuous(limits=c(0, 50)) +
  scale_x_continuous(limits=c(-2.5, 2.5)) +
  annotate("text",x=-1.5,y=45,label=down)+
  annotate("text",x=1.5,y=45,label=up)+
  annotate("text",x=0,y=45,label=nc)+
  theme_hc()
p
ggsave("volcanol_TE_ratio.pdf", width = 8, height = 8)

### 比较s2 codon数量
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
data <- read.csv("res_TE_cgbvsCOM.csv",header = T)
data$sig<- "no"  ### 筛选并填充数据
data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= 0.23] <- "up"
data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -0.23] <- "down"
View(data)
data$sig =factor(data$sig, levels = c("down","up","no"))

dfcount = read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/AT_codon_count2.csv",header = T)
df1 = merge(data,dfcount, by.x = "X",by.y = "ID")
library(ggpubr)
down = median(df1[df1$sig=="down","s2codon"])
up = median(df1[df1$sig=="up","s2codon"])
nc = median(df1[df1$sig=="no","s2codon"])

my_comparisons=list(c("down","up"))
la_p1 = 140
p1 <- ggplot(df1[df1$sig == "down"|df1$sig == "up",],aes(x=sig,y=s2codon,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "GAA&CAA&AAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p1
ggsave("mRNA_A_ending.pdf", width = 3, height = 5)


### GO 整体
rm(list=ls())
setwd("D:/BigData/cgb")
library(clusterProfiler)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)
library(Rmisc)
df_ratio_p = read.csv("res_TE_cgbvsCOM.csv",header = T)
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 ] <- "sig"
ratio_up_res =  df_ratio_p%>%
  filter(sig == "sig")%>%
  pull(X)%>%
  enrichGO( OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(cutoff = 0.6,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

ratio_up_res = as.data.frame(ratio_up_res)
ratio_up_res$B1 = as.numeric(str_extract(ratio_up_res$BgRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$B2 = as.numeric(str_extract(ratio_up_res$BgRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$BR = ratio_up_res$B1/ratio_up_res$B2
ratio_up_res$G1 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(.*)(?=(\\/))")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$G2 = as.numeric(str_extract(ratio_up_res$GeneRatio,"(?<=(\\/))(.*)")) #### R 的 正则表达式 提取 | 后的文本
ratio_up_res$GR = ratio_up_res$G1/ratio_up_res$G2
ratio_up_res$EnrichmentFold = ratio_up_res$GR/ratio_up_res$BR
ratio_up_res = arrange(ratio_up_res,-ratio_up_res$Count)
ratio_up_res = ratio_up_res[1:40,]
ratio_up_res = arrange(ratio_up_res,ratio_up_res$EnrichmentFold)
ratio_up_res$Description = factor(ratio_up_res$Description,levels = c(ratio_up_res$Description))
ggplot(ratio_up_res,aes(x=EnrichmentFold,y=Description, fill =-log(pvalue))) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)",
      #y="Pathway name",
      title="cgb vs COM")+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("Figure_TE_sig_col.pdf",width = 9,height= 9)

library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
library(clusterProfiler)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(org.At.tair.db)
dep_ratio_m = read.csv("TE_normalized.csv",header = T,row.names = 1)%>%
  as.data.frame()%>%
  dplyr::select(c(2:7))
dep_ratio_m1 = read.csv("TE_normalized.csv",header = T,row.names = 1)%>%
  as.data.frame()%>%
  dplyr::select(c(2:7))%>%
  apply(1,scale)%>%
  t()

colnames(dep_ratio_m1)=colnames(dep_ratio_m)
pheatmap(dep_ratio_m1,
         cluster_cols = F,
         cluster_rows = T,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = colorRampPalette(c("#003366","white","red"))(60), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "H2.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 2,
         height = 15)
pheatmap(dep_ratio_m1,
         cluster_cols = F,
         cluster_rows = T,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =1,
         filename = "H3.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 6,
         height = 15)


### correlation analysis
### TE degp analysis
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))
### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
### 整理表格生成heatmap 可以读取得 matrix
dep_ratio_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
dep_ratio_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)%>%
  left_join(dep_ratio_rna, by="ID")%>%
  mutate(cgb_pro_ratio=log(cgb_24_pro/cgb_0_pro),
              COM_pro_ratio=log(COM_24_pro/COM_0_pro),
              cgb_rna_ratio=log(cgb_24_rna/cgb_0_rna),
              COM_rna_ratio=log(COM_24_rna/COM_0_rna),
              cgb_pro_minus=(cgb_24_pro-cgb_0_pro),
              COM_pro_minus=(COM_24_pro-COM_0_pro),
              cgb_rna_minus=(cgb_24_rna-cgb_0_rna),
              COM_rna_minus=(COM_24_rna-COM_0_rna))


library(Rmisc)
library(ggpubr)
library(ggpmisc)
px1=ggplot(dep_ratio_pro,
       aes(x=cgb_rna_ratio, y=cgb_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  theme_classic2() +
  # geom_abline(slope = 0.162,intercept=0.107)+
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="cgb_ratio")
px1

px2=ggplot(dep_ratio_pro,
       aes(x=COM_rna_ratio, y=COM_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  # geom_abline(slope = 0.165,intercept=0.0803)+
  theme_classic2() +
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="COM_ratio")
px2

# px3=ggplot(dep_ratio_pro,
#        aes(x=log(cgb_0_rna), y=log(cgb_0_pro))) +
#   geom_hex(bins = 70) +
#   scale_fill_continuous(type = "viridis") +
#   stat_smooth(method="lm",color = "black", fill = "lightgray")+
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
#   # stat_cor(method = "spearman")+
#   # scale_y_continuous(limits=c(-8, 2))+
#   # scale_x_continuous(limits=c(-5, 13))+
#   theme_classic2() +
#   labs(
#     x="mRNA abundance",
#     y="protein abundance",
#     title="cgb_minus")
# px3
#
# px4=ggplot(dep_ratio_pro,
#        aes(x=log(COM_0_rna), y=log(COM_0_pro))) +
#   geom_hex(bins = 70) +
#   scale_fill_continuous(type = "viridis") +
#   stat_smooth(method="lm",color = "black", fill = "lightgray")+
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
#   # stat_cor(method = "spearman")+
#   # scale_y_continuous(limits=c(-8, 2))+
#   # scale_x_continuous(limits=c(-5, 13))+
#   theme_classic2() +
#   labs(
#     x="mRNA abundance",
#     y="protein abundance",
#     title="COM_minus")
# px4
library(cowplot)
pdf(file = "ggcorelation.pdf",10,4.5)
plot_grid(px1,px2, nrow = 1)
dev.off()



### 所有的差异因数量的描述性统计 饼图
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
setwd("D:/BigData/cgb/geom_col/og_data")

temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEG.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEG= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
  data <- read.csv(fil,header = T)%>%filter(log2FoldChange !="NA")
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= log2(1.5)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.5)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    group_by(sig) %>%
    summarize(count=n())%>%
    mutate(file = fil)%>%
    mutate(perc = count / sum(count)) %>%
    mutate(labels = scales::percent(perc))
  data_final_DEG = rbind(data_final_DEG,data)
}

temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEP.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEP= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
  data <- read.csv(fil,header = T)%>%filter(log2FoldChange !="NA")
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >=log2(1.2)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.2)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    group_by(sig) %>%
    summarize(count=n())%>%
    mutate(file = fil)%>%
    mutate(perc = count / sum(count)) %>%
    mutate(labels = scales::percent(perc))
  data_final_DEP = rbind(data_final_DEP,data)
}

data_final = rbind(data_final_DEG,data_final_DEP)

temp_0=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*.csv")
data_final$file=factor(data_final$file,
                       levels = c("res0H_DEG.csv",
  "res24H_DEG.csv",
  "res_T26_24hvs0h_DEG.csv",
  "res_15102_24hvs0h_DEG.csv",
  "res_ratio24to0_DEG.csv",
  "res0H_DEP.csv",
  "res24H_DEP.csv",
  "res_T26_24hvs0h_DEP.csv",
  "res_15102_24hvs0h_DEP.csv",
   "res_ratio24to0_DEP.csv"))

### 饼图
ggplot(data_final,aes(x="",y=perc, fill =sig))+
  facet_grid(~file,)+
  geom_col() +
    scale_fill_wsj('dem_rep')+
    coord_polar(theta = "y")+
    geom_text(aes(label = paste0(round(perc * 100, 0), "%")),
            position = position_fill(vjust = 0.5)) +
    # labs(
    #   # size="",
    #   # x=""
    #   #y="Pathway name",
    #   # title="Pathway enrichment")
    # )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
    )
ggsave("total_pie.pdf",width = 15,height = 3)


rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(dplyr)
setwd("D:/BigData/cgb/geom_col/og_data")

temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEG.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEG= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
  data <- read.csv(fil,header = T)%>%filter(log2FoldChange !="NA")
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= log2(1.5)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.5)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    group_by(sig) %>%
    summarize(count=n())%>%
    filter(sig != "no")%>%
    mutate(file = fil)
  data_final_DEG = rbind(data_final_DEG,data)
}

temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEP.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEP= data.frame(sig=list(),count=list(),file=list())

for (fil in temp) {
  data <- read.csv(fil,header = T)%>%filter(log2FoldChange !="NA")
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >=log2(1.2)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.2)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    group_by(sig) %>%
    summarize(count=n())%>%
    filter(sig != "no")%>%
    mutate(file = fil)
  data_final_DEP = rbind(data_final_DEP,data)
}
data_final = rbind(data_final_DEG,data_final_DEP)

temp_0=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*.csv")
data_final$file=factor(data_final$file,
                       levels = c("res0H_DEG.csv",
  "res24H_DEG.csv",
  "res_T26_24hvs0h_DEG.csv",
  "res_15102_24hvs0h_DEG.csv",
  "res_ratio24to0_DEG.csv",
  "res0H_DEP.csv",
  "res24H_DEP.csv",
  "res_T26_24hvs0h_DEP.csv",
  "res_15102_24hvs0h_DEP.csv",
   "res_ratio24to0_DEP.csv"))

ggplot(data_final,aes(x=sig,y=count, fill =sig ))+
    geom_col(aes())+
    facet_grid(~file,)+
    scale_fill_wsj('dem_rep')+
    labs(
      size="",
      x=""
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("total_col.pdf",width = 15,height = 3)


### boxplot logFC
rm(list=ls())
temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEG.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEG= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
  data <- read.csv(fil,header = T)
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >= log2(1.5)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.5)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    filter(sig != "no")%>%
    mutate(file = fil)
  data_final_DEG = rbind(data_final_DEG,data)
}

temp=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*DEP.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEP= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
  data <- read.csv(fil,header = T)
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >=log2(1.2)] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= -log2(1.2)] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))
  data = data%>%
    filter(sig != "no")%>%
    mutate(file = fil)
  data_final_DEP = rbind(data_final_DEP,data)
}
data_final = rbind(data_final_DEG,data_final_DEP)

temp_0=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*.csv")
data_final$file=factor(data_final$file,
                       levels = c("res0H_DEG.csv",
  "res24H_DEG.csv",
  "res_T26_24hvs0h_DEG.csv",
  "res_15102_24hvs0h_DEG.csv",
  "res_ratio24to0_DEG.csv",
  "res0H_DEP.csv",
  "res24H_DEP.csv",
  "res_T26_24hvs0h_DEP.csv",
  "res_15102_24hvs0h_DEP.csv",
   "res_ratio24to0_DEP.csv"))

ggplot(data_final,aes(x=sig,y=abs(log2FoldChange), fill =sig ))+
    geom_boxplot(position="dodge",outlier.colour = NA)+
    facet_grid(~file,)+
    scale_fill_wsj('dem_rep')+
    labs(
      size="",
      x=""
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    scale_y_continuous(limits=c(0,3))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("total_box_fc.pdf",width = 15,height = 3)


### 24H
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(dplyr)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

data_24_r <- read.csv("res24H_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(log2FoldChange>0)

data_24_p <- read.csv("res24H_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.263)%>%
  inner_join(data_24_r, by = "ID")%>%
  pull(ID)%>%
  enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(
  cutoff = 0.88,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df1 = as.data.frame(data_24_p)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 0),
      axis.title.y = element_blank()
    )
ggsave("24h_0525.pdf",width = 6,height =6)


### 24H
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(dplyr)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

data_24_r <- read.csv("res_ratio24to0_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(log2FoldChange>0)

data_24_p <- read.csv("res_ratio24to0_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.263)%>%
  inner_join(data_24_r, by = "ID")%>%
  pull(ID)%>%
  enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(
  cutoff = 0.88,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df1 = as.data.frame(data_24_p)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 0),
      axis.title.y = element_blank()
    )
ggsave("24h_0525.pdf",width = 6,height =6)




### correlation analysis
### TE degp analysis

rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
setwd("D:/BigData/cgb")

df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))
### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

# library(DESeq2)
# df= as.data.frame(dep_df)%>%dplyr::select(c(4,5,6,10,11,12))*1000000
# cts <- round(df)
# mutant <- factor(rep(c('cgb','COM'), each=3))
# colData <- data.frame(row.names=colnames(cts),mutant)
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = colData,
#                               design= ~ mutant)
# dds <- DESeq(dds)
# res01 <- results(dds, contrast=c("mutant","cgb","COM"))
# write.csv(as.data.frame(res01), file='res_DEP_24h_new.csv')
#
# library(DESeq2)
# df= as.data.frame(deg_df)%>%dplyr::select(c(4,5,6,10,11,12))
# cts <- round(df)
# mutant <- factor(rep(c('cgb','COM'), each=3))
# colData <- data.frame(row.names=colnames(cts),mutant)
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = colData,
#                               design= ~ mutant)
# dds <- DESeq(dds)
# res01 <- results(dds, contrast=c("mutant","cgb","COM"))
# write.csv(as.data.frame(res01), file='res_DEG_24h_new.csv')


### 整理表格生成heatmap 可以读取得 matrix
dep_ratio_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)

dep_ratio_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)%>%
  left_join(dep_ratio_rna, by="ID")%>%
  mutate(cgb_pro_ratio=log(cgb_24_pro/cgb_0_pro),
              COM_pro_ratio=log(COM_24_pro/COM_0_pro),
              cgb_rna_ratio=log(cgb_24_rna/cgb_0_rna),
              COM_rna_ratio=log(COM_24_rna/COM_0_rna),
              cgb_pro_minus=(cgb_24_pro-cgb_0_pro),
              COM_pro_minus=(COM_24_pro-COM_0_pro),
              cgb_rna_minus=(cgb_24_rna-cgb_0_rna),
              COM_rna_minus=(COM_24_rna-COM_0_rna))

### 筛选出上调表达COM 的差异蛋白 用于correlation
df001=read.csv("./geom_col/og_data/res_T26_24hvs0h_DEP.csv",header = T)%>%
filter(pvalue<0.05)%>%left_join(dep_ratio_pro)

library(Rmisc)
library(ggpubr)
library(ggpmisc)
px1=ggplot(df001,
       aes(x=cgb_rna_ratio, y=cgb_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  theme_classic2() +
  # geom_abline(slope = 0.162,intercept=0.107)+
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="cgb_ratio_dep")
px1

px2=ggplot(df001,
       aes(x=COM_rna_ratio, y=COM_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  # geom_abline(slope = 0.165,intercept=0.0803)+
  theme_classic2() +
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="COM_ratio_dep")
px2

px3=ggplot(dep_ratio_pro,
       aes(x=cgb_rna_ratio, y=cgb_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  theme_classic2() +
  # geom_abline(slope = 0.162,intercept=0.107)+
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="cgb_ratio_total")
px3

px4=ggplot(dep_ratio_pro,
       aes(x=COM_rna_ratio, y=COM_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-2, 2))+
  scale_x_continuous(limits=c(-5, 5))+
  # geom_abline(slope = 0.165,intercept=0.0803)+
  theme_classic2() +
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="COM_ratio_total")
px4

library(cowplot)
pdf(file = "ggcorelation.pdf",20,4.5)
plot_grid(px1,px2,px3,px4, nrow = 1)
dev.off()


### 所有的差异因数量的描述性统计及GO
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)

X = list(cgb_24vs0_deg=cgb_24vs0_deg,COM_24vs0_deg=COM_24vs0_deg)
p1=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_deg.pdf",width = 5,height = 2.5)

X = list(cgb_24vs0_dep=cgb_24vs0_dep,COM_24vs0_dep=COM_24vs0_dep)
p2=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_dep.pdf",width = 5,height = 2.5)

X = list(setdiff_com_deg=setdiff_com_deg,setdiff_com_dep=setdiff_com_dep)
p3=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_diff_com_dep_deg.pdf",width = 5,height = 2.5)

library(cowplot)
pdf(file = "venn_super.pdf",15,3.5)
plot_grid(p1,p2,p3, nrow = 1)
dev.off()

### 取差异集合
# setdiff(a, b) ### a的差异集合
# setdiff(b, a) ### b的差异集合
diff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)%>%
enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(
  cutoff = 0.68,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
diff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)%>%
enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(
  cutoff = 0.5,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
diff_final=setdiff(setdiff_com_dep,setdiff_com_deg)%>%
enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
  simplify(
  cutoff = 0.45,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
df1 = as.data.frame(diff_com_deg)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
p4=ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 0),
      axis.title.y = element_blank()
    )
df2 = as.data.frame(diff_com_dep)%>%arrange(-pvalue)
df2$Description=factor(df2$Description,levels=df2$Description)
p5=ggplot(df2,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 0),
      axis.title.y = element_blank()
    )
df3 = as.data.frame(diff_final)%>%arrange(-pvalue)
df3$Description=factor(df3$Description,levels=df3$Description)
p6=ggplot(df3,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 0),
      axis.title.y = element_blank()
    )
library(cowplot)
pdf(file = "go_com_setdiff.pdf",22,9)
plot_grid(p4,p5,p6, nrow = 1)
dev.off()

### 1149基因的相关性分析
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))

setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))
### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

### 整理表格生成heatmap 可以读取得 matrix
dep_ratio_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
dep_ratio_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)%>%
  right_join(dep_ratio_rna, by="ID")%>%COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%

  mutate(cgb_pro_ratio=log(cgb_24_pro/cgb_0_pro),
              COM_pro_ratio=log(COM_24_pro/COM_0_pro),
              cgb_rna_ratio=log(cgb_24_rna/cgb_0_rna),
              COM_rna_ratio=log(COM_24_rna/COM_0_rna),
              cgb_pro_minus=(cgb_24_pro-cgb_0_pro),
              COM_pro_minus=(COM_24_pro-COM_0_pro),
              cgb_rna_minus=(cgb_24_rna-cgb_0_rna),
              COM_rna_minus=(COM_24_rna-COM_0_rna))
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))%>%
  left_join(dep_ratio_pro,by="ID")

px3=ggplot(df_1149,
       aes(x=cgb_rna_ratio, y=cgb_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-1.5, 1.5))+
  scale_x_continuous(limits=c(-5, 5))+
  theme_classic2() +
  # geom_abline(slope = 0.162,intercept=0.107)+
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="cgb_ratio_total")
px3
px4=ggplot(df_1149,
       aes(x=COM_rna_ratio, y=COM_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-1.5, 1.5))+
  scale_x_continuous(limits=c(-5, 5))+
  # geom_abline(slope = 0.165,intercept=0.0803)+
  theme_classic2() +
  labs(
    x="mRNA abundance",
    y="protein abundance",
    title="COM_ratio_total")
px4
library(cowplot)
pdf(file = "correlation_1149.pdf",9,4)
plot_grid(px3,px4, nrow = 1)
dev.off()




### 1149 蛋白s2codon的比较
### 1149基因的相关性分析
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))

df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))%>%mutate(type="1149")
# df_2537= data.frame(ID= intersect(COM_24vs0_dep,cgb_24vs0_dep))%>%mutate(type="2537")
# df_new = rbind(df_1149,df_2537)

data <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  mutate(sig ="no")%>%left_join(df_1149,by="ID")%>%mutate(type = gsub("na", "other",type) )

data$sig[data$type == "NA"] <- "up"
data$sig[data$type == "1149"] <- "1149"

dfcount = read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/AT_codon_count2.csv",header = T)%>%
  right_join(data,by="ID")%>%filter(AAA != "NA",CAA!="NA",GAA!="NA")

library(ggpubr)
my_comparisons=list(c("no","1149"))
la_p1 = 50
p1 <- ggplot(dfcount,aes(x=sig,y=AAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "AAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p1

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","AAA"]))-median(c(df034[df034$sig == "down","AAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$AAA,N1,replace = FALSE))) - mean(c(sample(df034$AAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p11 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p11

la_p1 =30
p2 <- ggplot(dfcount,aes(x=sig,y=CAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "CAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p2

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","CAA"]))-median(c(df034[df034$sig == "down","CAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$CAA,N1,replace = FALSE))) - mean(c(sample(df034$CAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p12 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p12

la_p1 = 60
p3 <- ggplot(dfcount,aes(x=sig,y=GAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "GAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p3

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","GAA"]))-median(c(df034[df034$sig == "down","GAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$GAA,N1,replace = FALSE))) - mean(c(sample(df034$GAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p13 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p13

library(cowplot)
pdf(file = "mRNA_A_ending.pdf",7,4)
plot_grid(p3,p2,p1, nrow = 1)
dev.off()

library(cowplot)
pdf(file = "sampling.pdf",7,4)
plot_grid(p13,p12,p11, nrow = 1)
dev.off()

### 1149的差异表达GO及热图
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
library(stringr)

### 获得1149 蛋白
setwd("D:/BigData/cgb/geom_col/og_data")
cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))
data <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%right_join(df_1149,by="ID")
  data$sig<- "no"  ### 筛选并填充数据
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange >=0] <- "up"
  data$sig[data$pvalue <= 0.05 & data$log2FoldChange <= 0] <- "down"
  data$sig =factor(data$sig, levels = c("down","up","no"))

### 1149生成蛋白matrix
setwd("D:/BigData/cgb")
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

dep_df_com = dep_df%>%
mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,COM_0_pro,COM_24_pro,cgb_0_pro,cgb_24_pro)%>%right_join(data,by="ID")
row.names(dep_df_com) = dep_df_com$ID

up_1149_GO= dep_df_com %>%
  filter(sig == "up")%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
 simplify(
  cutoff = 0.45,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)
df1 = as.data.frame(up_1149_GO)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank())
ggsave("up_1149_GO.pdf",width = 6,height =6)

down_1149_GO= dep_df_com %>%
  filter(sig == "down")%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
 simplify(
  cutoff = 0.55,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
df1 = as.data.frame(down_1149_GO)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("down_1149_GO.pdf",width = 4,height =6)

library(pheatmap)
library(viridis)
matrix = dep_df_com%>%dplyr::select(c(2:5))%>%
arrange(COM_0_pro)%>%
apply(1,scale)%>%
t()%>%
pheatmap(
         cluster_cols = F,
         cluster_rows = T,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =2,
         filename = "1149_heatmap.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 3.5,
         height = 12)
table(data$sig)


### 叶绿体蛋白组的变化
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(stringr)
setwd("D:/BigData/cgb")
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)%>%
  mutate(chrom = str_extract(ID,"...."))%>%
  filter(chrom =="ATCG")
row.names(dep_df) = dep_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID","chrom")
dep_df_chlo = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)

### org注释更好 20220526
library(org.At.tair.db)
df_name = select(org.At.tair.db, keys=dep_df_chlo$ID, columns=c("GENENAME","SYMBOL","TAIR"), keytype="TAIR")
dep_df_chlo =  dep_df_chlo%>%
  left_join(df_name,by=c("ID"="TAIR"))%>%
  distinct(ID,.keep_all = T)
row.names(dep_df_chlo) = dep_df_chlo$SYMBOL

library(pheatmap)
library(viridis)
dep_df_chlo%>%
dplyr::select(c(2:5))%>%
apply(1,scale)%>%
pheatmap(
         cluster_cols = T,
         cluster_rows = F,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =1,
         filename = "chloro.pdf",   ### ***
         display_numbers = F,
         show_colnames = T,
         show_rownames = T,
         width = 12,
         height = 6)

### Human S2 codon
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(DOSE)
setwd("D:/BigData/cgb")
df = read.csv("D:/BigData/length/og_data_go/24Homo_sapiens.csv",header = T)%>%
  group_by(GeneID) %>%
  summarize(count=n(),
            mea=mean(Length,na.rm=T),
            med=median(Length,na.rm=T))%>%
  arrange(-mea)%>%
  dplyr::slice(c(1:975))

res=df%>%
  pull(GeneID)%>%
  enrichGO()
write.csv(as.data.frame(res),"GO_Hs_S2.csv")

ggsave("Hs.pdf", width = 40, height = 40,limitsize = FALSE)
ggsave("Hs_20.pdf", width = 20, height = 20,limitsize = FALSE)


### cluster 多分类
#### HEATMAP
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)

setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]

### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

### 整理表格生成heatmap 可以读取得 matrix
deg_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)

dep_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)


### 获得差异基因蛋白
setwd("D:/BigData/cgb/geom_col/og_data")
cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
deg = as.data.frame(union(cgb_24vs0_deg,COM_24vs0_deg))
colnames(deg)="ID"

cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
dep = as.data.frame(union(cgb_24vs0_dep,COM_24vs0_dep))
colnames(dep)="ID"

deg = deg%>%
  left_join(deg_rna,by="ID")
dep = dep%>%
  left_join(dep_pro,by="ID")

row.names(deg)=deg$ID
deg=deg%>%
  dplyr::select(cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
row.names(dep)=dep$ID
dep=dep%>%
  dplyr::select(cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)

library(factoextra)
library(cluster)
setwd("D:/BigData/cgb")
setwd("D:/BigData/cgb")
deg1=pam(deg, 4, metric = "euclidean", stand = FALSE)
fviz_nbclust(deg, pam, method = "wss")
clust = deg1$clustering%>%
  as.data.frame()

clust$ID=row.names(clust)
deg = deg_rna%>%
  right_join(clust,by ="ID")
write.csv(deg,"cluster_deg.csv")

###manhattan cluster
rm(list=ls())
library(dplyr)
library(tidyr)
library(clusterProfiler)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_cluster = read.csv("cluster_deg.csv",header=T)%>%select(ID,COM_0_rna,COM_24_rna,cgb_0_rna,cgb_24_rna,cluster)
table(df_cluster$cluster)
library(pheatmap)
library(viridis)

df=df_cluster%>%
dplyr::select(c(2:5))%>%
apply(1,scale)%>%
t()
rownames(df)=df_cluster$ID
colnames(df)=c("COM_0_rna","COM_24_rna","cgb_0_rna","cgb_24_rna")

list=pheatmap(df,
         cluster_cols = F,
         cluster_rows = T,
         scale = "none",
         border=FALSE,
         # kmeans_k=4,
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =6,
         clustering_distance_rows = "manhattan",
         filename = "cluster_mean_manhattan.pdf",
         display_numbers = F,
         show_colnames = T,
         show_rownames = F,
         width = 4,
         height = 15)

row_cluster=cutree(list$tree_row,k=6)
ID = row.names(as.data.frame(row_cluster))
table(cluster_go_1$row_cluster)

cluster_go_1=as.data.frame(row_cluster)%>%
  mutate(ID=ID)%>%
  filter(row_cluster == 6)%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.48,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
write.csv(as.data.frame(cluster_go_1),"cluster_go_6.csv")

library(ggplot2)
df1 = as.data.frame(cluster_go_1)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =-log(pvalue) )) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="Log(Pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("clust1.pdf",width = 6,height =6)

### 1149的差异表达GO及热图
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
library(stringr)
### 获得1149 蛋白
setwd("D:/BigData/cgb/geom_col/og_data")
cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))
### 1149生成蛋白matrix
setwd("D:/BigData/cgb")
up_1149_GO= df_1149 %>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
 simplify(
  cutoff = 0.55,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)


library(ggplot2)
df1 = as.data.frame(up_1149_GO)%>%
  arrange(-pvalue)%>%
  filter(ID!="GO:0044282")%>%
  filter(ID!="GO:0090066")%>%
  filter(ID!="GO:0009657")%>%
  filter(ID!="GO:0044087")%>%
  filter(ID!="GO:0010243")%>%
  filter(ID!="GO:1901361")%>%
  filter(ID!="GO:0018208")%>%
  filter(ID!="GO:1902600")%>%
  filter(ID!="GO:0052546")%>%
  filter(ID!="GO:0018208")%>%
  filter(ID!="GO:0019439")%>%
  filter(ID!="GO:0031503")%>%
  filter(ID!="GO:0009150")%>%
  filter(ID!="GO:0006091")%>%
  filter(Description !="positive regulation of transcription initiation from RNA polymerase II promoter")%>%
  filter(Description !="regulation of transcription initiation from RNA polymerase II promoter")
df1$Description=factor(df1$Description,levels=df1$Description)
ggplot(df1,aes(x=-log(pvalue),y=Description)) +
    geom_col(aes())+
    scale_fill_gradient(low="#158bb8",high="#cc163a")+
    labs(
      size="",
      x="pvalue"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("1149_GO.pdf",width = 5.5,height =6)
write.csv(df1,"1149_go.csv")

### SA相关蛋白的词云
setwd("D:/BigData/cgb/geom_col/og_data")
data <- read.csv("res_TE_cgbvsCOM.csv",header = T)
str_SA = unlist(strsplit(df1["GO:0009627","geneID"], '/', fixed=TRUE))
str_Chloro = unlist(strsplit(df1["GO:0045037","geneID"], '/', fixed=TRUE))
### 注释
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
str_Chloro =  as.data.frame(str_Chloro)%>%
  left_join(df_anno,by=c("str_Chloro"="TAIR"))
str_immune = unlist(strsplit(df1["GO:0002376","geneID"], '/', fixed=TRUE))

### 注释
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df_immune =  as.data.frame(str_immune)%>%
  left_join(df_anno,by=c("str_immune"="TAIR"))%>%
  left_join(data,by=c("str_immune"="X"))%>%
  dplyr::select("name","log2FoldChange")%>%
  filter(name !="")%>%
  mutate(value=abs(log2FoldChange))%>%
  dplyr::select("name","value")
my_graph=wordcloud2(df_immune, size = 1, minSize = 0, gridSize =  0,
           fontFamily = "微软雅黑", fontWeight = 'normal',
           backgroundColor = "white",
           minRotation = -pi/4, maxRotation = pi/4, rotateRatio = 0.4,
           shape = 'circle', ellipticity = 0.65, widgetsize = NULL)

library(webshot)
library(htmlwidgets)
saveWidget(my_graph,"tmp.html",selfcontained = F)
webshot("tmp.html","world_cloud_immune.pdf", delay =5, vwidth = 720, vheight=480)

library(ggwordcloud)
library(webshot)

# Make the graph
my_graph= wordcloud2(df_immune, size = 1, minSize = 0, gridSize =  0,
           fontFamily = "微软雅黑", fontWeight = 'normal',
           backgroundColor = "white",
           minRotation = -pi/4, maxRotation = pi/4, rotateRatio = 0.4,
           shape = 'circle', ellipticity = 0.65, widgetsize = NULL)

# save it in html
library("htmlwidgets")
saveWidget(my_graph,"tmp.html",selfcontained = F)
# and in pdf
webshot("tmp.html","world_cloud.pdf", delay =5, vwidth = 720, vheight=480)


### 1149的相关蛋白的词云
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
library(clusterProfiler)
library(stringr)

setwd("D:/BigData/cgb/geom_col/og_data")
cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  pull(ID)
setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)

data <- read.csv("res_TE_cgbvsCOM.csv",header = T)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)

df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))%>%
  left_join(data,by=c("ID"="X"))%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select("name","log2FoldChange")%>%
  filter(name !="")%>%
  mutate(value=abs(log2FoldChange)*100)%>%
  dplyr::select("name","value")

my_graph=wordcloud2(df_1149, size = 1, minSize = 0, gridSize =  0,
           fontFamily = "微软雅黑", fontWeight = 'normal',
           backgroundColor = "white",
           minRotation = -pi/4, maxRotation = pi/4, rotateRatio = 0.4,
           shape = 'circle', ellipticity = 0.65, widgetsize = NULL)

library(webshot)
library(htmlwidgets)
saveWidget(my_graph,"tmp.html",selfcontained = F)
webshot("tmp.html","world_cloud_total.pdf", delay =5, vwidth = 720, vheight=480)










### ratio new
### ratio 数据分析
rm(list=ls())
setwd("D:/BigData/cgb")
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(org.At.tair.db)
library(ggplot2)

dfdeg=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEGratio24to0.csv",header = T)
dfdep=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEPratio24to0.csv",header = T)

# dfdeg$ID= dfdeg$X
# dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
cgb_UG <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
cgb_DG <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
cgb_UP <- dfdeg[dfdeg$log2FoldChange > 0.584 & dfdeg$pvalue < 0.05,]
cgb_DP <- dfdeg[dfdeg$log2FoldChange < -0.584 & dfdeg$pvalue < 0.05,]

###1
cgb_UG_go=cgb_UG %>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.51,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df1 = as.data.frame(cgb_UG_go)%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "up"

###2
cgb_DG_go=cgb_DG %>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.56,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df2 = as.data.frame(cgb_DG_go)%>%
  arrange(-pvalue)%>%
  filter(Description != "immune response")
df2$Description=factor(df2$Description,levels=df2$Description)
df2$sig = "down"
df = rbind (df1,df2)

ggplot(df,aes(x=-log(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size="",
      x="Enrichment"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("cgb_DEG_go.pdf",width = 6,height =12)

###3
cgb_UP_go=cgb_UP %>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.5,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df1 = as.data.frame(cgb_UP_go)%>%
  filter(Description != "response to wounding")%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "less upregulated"


###4
cgb_DP_go=cgb_DP %>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.62,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)

df2 = as.data.frame(cgb_DP_go)%>%
  arrange(-pvalue)
df2$Description=factor(df2$Description,levels=df2$Description)
df2$sig = "less downregulated"
df = rbind (df1,df2)

ggplot(df,aes(x=-log(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size="",
      x="Enrichment"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("cgb_DEP_go.pdf",width = 6,height =12)



### 0H TE 数据分析
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]

### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

### 整理表格生成heatmap 可以读取得 matrix
deg_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,COM_0_rna)

dep_pro= dep_df%>%
  dplyr::select("cgb_0_1","cgb_0_2","cgb_0_3",
                  "COM_0_1","COM_0_2","COM_0_3","ID")%>%
  left_join(deg_rna,by="ID")%>%
  filter(cgb_0_rna !=0,COM_0_rna!=0)%>%
  mutate(cgb_0_1 = cgb_0_1/cgb_0_rna*100000000,
         cgb_0_2 = cgb_0_2/cgb_0_rna*100000000,
         cgb_0_3 = cgb_0_3/cgb_0_rna*100000000)%>%
  mutate(COM_0_1 = COM_0_1/COM_0_rna*100000000,
           COM_0_2 = COM_0_2/COM_0_rna*100000000,
           COM_0_3 = COM_0_3/COM_0_rna*100000000)%>%
  dplyr::select("cgb_0_1","cgb_0_2","cgb_0_3",
                  "COM_0_1","COM_0_2","COM_0_3","ID")

rownames(dep_pro) = dep_pro$ID
library(DESeq2)
cts00 <- dep_pro[,c(1:6)]
cts <- round(cts00)
mutant <- factor(rep(c('cgb','COM'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","cgb","COM"))
write.csv(as.data.frame(res01), file='res_0H_TE.csv')

data = as.data.frame(res01)
cgb_UP <- data[data$log2FoldChange > 0.23 & data$pvalue < 0.05,]
cgb_DP <- data[data$log2FoldChange < -0.23 & data$pvalue < 0.05,]


###
cgb_UP_go=cgb_UP %>%
  mutate(ID=rownames(cgb_UP))%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.6,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
df1 = as.data.frame(cgb_UP_go)%>%
  filter(Description != "response to wounding")%>%
  arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "TE up"


###
cgb_DP_go=cgb_DP %>%
  mutate(ID=rownames(cgb_DP))%>%
  dplyr::pull(ID)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.5,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
df2 = as.data.frame(cgb_DP_go)%>%
  arrange(-pvalue)
df2$Description=factor(df2$Description,levels=df2$Description)
df2$sig = "TE down"
df = rbind (df1,df2)
df = rbind (df1,df2)
df$sig =factor(df$sig, levels = c("TE up","TE down"))

ggplot(df2,aes(x=-log10(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size="",
      x="-log10(pvalue)"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("cgb_0H_TE_go.pdf",width = 6,height =12)

str_Chloro = unlist(strsplit(df["GO:0033365","geneID"], '/', fixed=TRUE))
### 注释
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
str_Chloro =  as.data.frame(str_Chloro)%>%
  left_join(df_anno,by=c("str_Chloro"="TAIR"))



#### ICS1, EDS5, PBS3, SARD1, CBP60g
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(ggthemes)

setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]

### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

### 整理表格生成heatmap 可以读取得 matrix
deg_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
dep_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)

deg_df_SA=deg_df %>%
  filter(ID == "AT4G39030"|ID == "AT1G74710"|ID =="AT5G26920"|ID =="AT5G13320"|ID =="AT1G73805")%>%
  mutate(name=c("EDS5","ICS1","CBP60G","PBS3","SARD1"))
rownames(deg_df_SA) = c("EDS5","ICS1","CBP60G","PBS3","SARD1")
deg_df_SA=deg_df_SA%>%t()%>%as.data.frame()%>%dplyr::slice(c(1:12))%>%
mutate(geno= rep(c("cgb_0","cgb_24","COM_0","COM_24"), each=3))%>%
mutate(type = "RNA")%>%
gather(key="Gene", value="Reads",c(1:5))%>%
  mutate(Reads=as.numeric(Reads))

dep_df_SA=dep_df%>%
  filter(ID == "AT4G39030"|ID == "AT1G74710"|ID =="AT5G26920"|ID =="AT5G13320"|ID =="AT1G73805")%>%
  mutate(name=c("EDS5","ICS1","CBP60G","PBS3","SARD1"))
rownames(dep_df_SA) = c("EDS5","ICS1","CBP60G","PBS3","SARD1")
dep_df_SA=dep_df_SA%>%t()%>%as.data.frame()%>%dplyr::slice(c(1:12))%>%
mutate(geno= rep(c("cgb_0","cgb_24","COM_0","COM_24"), each=3))%>%
mutate(type = "PRO")%>%
gather(key="Gene", value="Reads",c(1:5))%>%
  mutate(Reads=as.numeric(Reads))

library(ggplot2)
ggplot(deg_df_SA,aes(x=geno,y=Reads,fill=geno))+
    geom_boxplot(position="dodge",outlier.colour = NA)+
  facet_grid(~Gene,scales = "free")+
    labs(
      size="",
      x=""
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    theme_classic()+
    scale_fill_wsj() +
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )

ggsave("deg_SA.pdf",width = 6,height = 3.5)

library(ggplot2)
ggplot(dep_df_SA,aes(x=geno,y=Reads,fill=geno))+
    geom_boxplot(position="dodge",outlier.colour = NA)+
  facet_grid(~Gene,scales = "free")+
    labs(
      size="",
      x=""
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    theme_classic()+
    scale_fill_wsj() +
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1),angle = 65),
      axis.title.y = element_blank()
    )
ggsave("dep_SA.pdf",width = 6,height = 3.5)

### 叶绿体全蛋白组
### 从数据库中提取序列计算不同 BP/CC transcript length
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggthemes)
library(clusterProfiler)
library(org.EcK12.eg.db)
library(org.Sc.sgd.db)
library(org.At.tair.db)
library(org.Hs.eg.db)
library(org.Ce.eg.db)
library(GO.db)
columns(GO.db)
setwd("D:/BigData/length/og_data_go")
rm(list=ls())
db = org.At.tair.db
#1 extract all go terms
df01 = select(db,
              keys=keys(db, keytype="ENTREZID"),
              columns=c("GOALL","TAIR"),
              keytype="ENTREZID")%>%
  filter(ONTOLOGYALL == "CC")%>%
  distinct(GOALL, .keep_all= TRUE)
df06 = merge(select(GO.db, keys=df01$GOALL, columns=c("TERM","ONTOLOGY"), keytype="GOID"),
  df01,
  by.x = "GOID",
  by.y ="GOALL")
write.csv(df06,sav_file)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df_chl = select(db,
              keys=keys(db, keytype="ENTREZID"),
              columns=c("GOALL","TAIR"),
              keytype="ENTREZID")%>%
  filter(ONTOLOGYALL == "CC")%>%
  filter(GOALL=="GO:0009507")%>%
  distinct(TAIR, .keep_all= TRUE)%>%
  left_join(df_anno,by="TAIR")%>%
  dplyr::select(c(5:9))
write.csv(df_chl,"chloroplast_all_anno.csv")

### DEP/DEG matrix 模块
rm(list=ls())
setwd("D:/BigData/cgb")
df_chl = read.csv("chloroplast_all_anno.csv",row.names = 1)
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
deg_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
dep_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)

df_chl_deg = df_chl%>%
  left_join(deg_rna,by=c("TAIR"="ID"))%>%
  filter(cgb_0_rna !="NA")%>%
  filter(cgb_24_rna !="NA")%>%
  filter(COM_0_rna !="NA")%>%
  filter(COM_24_rna !="NA")%>%
  filter(cgb_0_rna !="NaN")%>%
  filter(cgb_24_rna !="NaN")%>%
  filter(COM_0_rna !="NaN")%>%
  filter(COM_24_rna !="NaN")%>%
  filter(cgb_0_rna !=Inf)%>%
  filter(cgb_24_rna !=Inf)%>%
  filter(COM_0_rna !=Inf)%>%
  filter(COM_24_rna !=Inf)%>%
  dplyr::select(cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)%>%
  apply(1,scale)%>%
  t()
 df_chl_deg <- df_chl_deg[complete.cases(df_chl_deg), ]## complete case
 colnames(df_chl_deg)=c("cgb_0_rna","cgb_24_rna","COM_0_rna","COM_24_rna")

 df_chl_dep = df_chl%>%
  left_join(dep_pro,by=c("TAIR"="ID"))%>%
  filter(cgb_0_pro !="NA")%>%
  dplyr::select(cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)%>%
  apply(1,scale)%>%
  t()
colnames(df_chl_dep)=c("cgb_0_pro","cgb_24_pro","COM_0_pro","COM_24_pro")

library(pheatmap)
library(viridis)
pheatmap(df_chl_deg,
         cluster_cols = F,
         cluster_rows = T,
         scale = "none",
         border=FALSE,
         # cellwidth = 50,
         # cellheight = 25,         main = "",
         color = inferno(20), #### B2182B
         cutree_cols = 1,
         cutree_rows =1,
         filename = "df_chl_deg.pdf",   ### ***
         display_numbers = F,
         show_colnames =T,
         show_rownames = F,
         width = 6,
         height = 12)

### 叶绿体全蛋白组的差异数量
rm(list=ls())
setwd("D:/BigData/cgb/geom_col/og_data")
deg_0H <- read.csv("res0H_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)
dep_0H <- read.csv("res0H_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)

setwd("D:/BigData/cgb")
deg_df_chl = read.csv("chloroplast_all_anno.csv",row.names = 1)%>%
  left_join(deg_0H,by=c("TAIR"="ID"))%>%
  filter(baseMean !="NA")
deg_df_chl$sig<- "no"  ### 筛选并填充数据
deg_df_chl$sig[deg_df_chl$pvalue <= 0.05 & deg_df_chl$log2FoldChange >= 0] <- "up"
deg_df_chl$sig[deg_df_chl$pvalue <= 0.05 & deg_df_chl$log2FoldChange <= 0] <- "down"
table(deg_df_chl$sig)
dep_df_chl = read.csv("chloroplast_all_anno.csv",row.names = 1)%>%
  left_join(dep_0H,by=c("TAIR"="ID"))%>%
  filter(baseMean !="NA")
dep_df_chl$sig<- "no"  ### 筛选并填充数据
dep_df_chl$sig[dep_df_chl$pvalue <= 0.05 & dep_df_chl$log2FoldChange >= 0] <- "up"
dep_df_chl$sig[dep_df_chl$pvalue <= 0.05 & dep_df_chl$log2FoldChange <= 0] <- "down"
table(dep_df_chl$sig)

### 处理表格并添加注释
rm(list=ls())
library(dplyr)
setwd("D:/BigData/cgb/geom_col/og_data")
df <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.584] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.584] <- "down"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("cgb_24hvs0h_DEG.csv",row.names = F)

df <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.584] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.584] <- "down"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("COM_24hvs0h_DEG.csv",row.names = F)

df <- read.csv("res_ratio24to0_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.584] <- "less downregulated"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.584] <- "less upregulated"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("DEGratio24to0.csv",row.names = F)

df <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.23] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.23] <- "down"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("cgb_24hvs0h_DEP.csv",row.names = F)

df <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.23] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.23] <- "down"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("COM_24hvs0h_DEP.csv",row.names = F)

df <- read.csv("res_ratio24to0_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= 0.23] <- "less downregulated"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -0.23] <- "less upregulated"
df%>%filter(sig != "no")%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("DEPratio24to0.csv",row.names = F)

###1149个蛋白注释
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggVennDiagram)
setwd("D:/BigData/cgb/geom_col/og_data")

df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(abs(log2FoldChange) >= 0.584)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(abs(log2FoldChange) >= 0.584)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(abs(log2FoldChange) >= 0.23)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  filter(abs(log2FoldChange) >= 0.23)%>%
  pull(ID)
setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_870= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("870_anno.csv",row.names = F)
X = list(COM_24vs0_deg=COM_24vs0_deg,cgb_24vs0_deg=cgb_24vs0_deg)
p1=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_deg.pdf",width = 5,height = 2.5)
X = list(COM_24vs0_dep=COM_24vs0_dep,cgb_24vs0_dep=cgb_24vs0_dep)
p2=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_dep.pdf",width = 5,height = 2.5)
X = list(setdiff_com_deg=setdiff_com_deg,setdiff_com_dep=setdiff_com_dep)
p3=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
ggsave("venn_diff_com_dep_deg.pdf",width = 5,height = 2.5)

library(cowplot)
pdf(file = "venn_super.pdf",15,3.5)
plot_grid(p1,p2,p3, nrow = 1)
dev.off()


### Venn 图  up regulated
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
# library(ComplexHeatmap)
# library(circlize)
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "-",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))

# dfdeg=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEGratio24to0.csv",header = T)
# dfdep=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEPratio24to0.csv",header = T)
df_deg_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_15102_24hvs0h.CSV",header = T)
df_deg_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEG_T26_24hvs0h.CSV",header = T)
df_dep_15=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_15102_24hvs0h.csv",header = T)
df_dep_26=read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/res_DEP_T26_24hvs0h.csv",header = T)
dfdeg = df_deg_15
dfdep = df_dep_15
# dfdeg$ID= dfdeg$X
# dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]%>%
  filter(!is.na(log2FoldChange))
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]%>%
  filter(!is.na(log2FoldChange))
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 0.584 & dfdeg$pvalue < 0.05,]%>%
  filter(!is.na(log2FoldChange))
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -0.584 & dfdeg$pvalue < 0.05,]%>%
  filter(!is.na(log2FoldChange))

length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)
cgb_UG = dfdeg_up
cgb_DG = dfdeg_down
cgb_UP = dfdep_up
cgb_DP = dfdep_down

dfdeg = df_deg_26
dfdep = df_dep_26
# dfdeg$ID= dfdeg$X
# dfdep$ID= str_extract(dfdep$X,"(.*)(?=\\.\\d)")  #### R 的 正则表达式 提取 .1  之前的
dfdep_up <- dfdep[dfdep$log2FoldChange > 0.23 & dfdep$pvalue < 0.05,]
dfdep_down <- dfdep[dfdep$log2FoldChange < -0.23 & dfdep$pvalue < 0.05,]
dfdeg_up <- dfdeg[dfdeg$log2FoldChange > 1 & dfdeg$pvalue < 0.05,]
dfdeg_down <- dfdeg[dfdeg$log2FoldChange < -1 & dfdeg$pvalue < 0.05,]
length(dfdep_up$ID)
length(dfdep_down$ID)
length(dfdeg_up$ID)
length(dfdeg_down$ID)
COM_UG = dfdeg_up
COM_DG = dfdeg_down
COM_UP = dfdep_up
COM_DP = dfdep_down
library(ggVennDiagram)
library(ggvenn)
X_COM_up = list(COM_UP=COM_UP$ID,COM_UG = COM_UG$ID)
X_cgb_up = list(cgb_UP=cgb_UP$ID,cgb_UG = cgb_UG$ID)
X_COM_down = list(COM_DP=COM_DP$ID,COM_DG = COM_DG$ID)
X_cgb_down = list(cgb_DP=cgb_DP$ID,cgb_DG = cgb_DG$ID)

p1=ggVennDiagram(X_COM_up)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p2=ggVennDiagram(X_cgb_up)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p3=ggVennDiagram(X_COM_down)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p4=ggVennDiagram(X_cgb_down)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

X_COM_cgb_up_pro = list(COM_UP=COM_UP$ID,cgb_UP = cgb_UP$ID)
X_COM_cgb_up_rna = list(COM_UG=COM_UG$ID,cgb_UG = cgb_UG$ID)

p1=ggVennDiagram(X_COM_cgb_up_pro)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
p2=ggVennDiagram(X_COM_cgb_up_rna)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

library(cowplot)
pdf(file = "ggvenn.pdf",20,17)
plot_grid(p1,p2,nrow = 2)
dev.off()



### up gene 442做图五
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggpubr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.584)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.584)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.23)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.23)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)

X = list(COM_24vs0_deg=COM_24vs0_deg,cgb_24vs0_deg=cgb_24vs0_deg)
p1=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

X = list(COM_24vs0_dep=COM_24vs0_dep,cgb_24vs0_dep=cgb_24vs0_dep)
p2=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

X = list(setdiff_com_dep=setdiff_com_dep,setdiff_com_deg=setdiff_com_deg)
p3=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()

library(cowplot)
pdf(file = "venn_super_UP.pdf",15,3.5)
plot_grid(p1,p2,p3, nrow = 1)
dev.off()


### 442 correlation
setwd("D:/BigData/cgb")
df_rna=read.csv("mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))
### 合并所有 DEG/DEP
dep_df = df_pro[,2:14]%>%
  distinct(ID,.keep_all = T)
row.names(dep_df) = dep_df[,13]
deg_df = df_rna[,2:14]
row.names(deg_df) = deg_df[,13]
### 整理表格生成heatmap 可以读取得 matrix
colnames(dep_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")
colnames(deg_df) = c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3",
                  "COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3","ID")

### 整理表格生成heatmap 可以读取得 matrix
dep_ratio_rna = deg_df%>%
  mutate(cgb_0_rna=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_rna=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_rna=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_rna=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_rna,cgb_24_rna,COM_0_rna,COM_24_rna)
dep_ratio_pro = dep_df%>%
  mutate(cgb_0_pro=(cgb_0_1+cgb_0_2+cgb_0_3)/3,
              cgb_24_pro=(cgb_24_1+cgb_24_2+cgb_24_3)/3,
              COM_0_pro=(COM_0_1+COM_0_2+COM_0_3)/3,
              COM_24_pro=(COM_24_1+COM_24_2+COM_24_3)/3)%>%
  dplyr::select(ID,cgb_0_pro,cgb_24_pro,COM_0_pro,COM_24_pro)%>%
  left_join(dep_ratio_rna, by="ID")%>%
  mutate(cgb_pro_ratio=log(cgb_24_pro/cgb_0_pro),
              COM_pro_ratio=log(COM_24_pro/COM_0_pro),
              cgb_rna_ratio=log(cgb_24_rna/cgb_0_rna),
              COM_rna_ratio=log(COM_24_rna/COM_0_rna),
              cgb_pro_minus=(cgb_24_pro-cgb_0_pro),
              COM_pro_minus=(COM_24_pro-COM_0_pro),
              cgb_rna_minus=(cgb_24_rna-cgb_0_rna),
              COM_rna_minus=(COM_24_rna-COM_0_rna))

library(ggpmisc)

df_442= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))%>%
  left_join(dep_ratio_pro,by="ID")%>%
  filter(cgb_rna_ratio>-2.5&cgb_rna_ratio<2.5)
px3=ggplot(df_442,
       aes(x=cgb_rna_ratio, y=cgb_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-1.5, 1.5))+
  scale_x_continuous(limits=c(-2.5, 2.5))+
  theme_classic2() +
  # geom_abline(slope = 0.162,intercept=0.107)+
  labs(
    x="mRNA change",
    y="protein change",
    title="cgb_ratio_total")
px3

px4=ggplot(df_442,
       aes(x=COM_rna_ratio, y=COM_pro_ratio)) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  stat_smooth(method="lm",color = "black", fill = "lightgray")+
  stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T, label.x = "middle")+
  # stat_cor(method = "spearman")+
  scale_y_continuous(limits=c(-1.5, 1.5))+
  scale_x_continuous(limits=c(-2.5, 2.5))+
  # geom_abline(slope = 0.165,intercept=0.0803)+
  theme_classic2() +
  labs(
    x="mRNA change",
    y="protein change",
    title="COM_ratio_total")
px4

library(cowplot)
pdf(file = "correlation_df_442.pdf",9,4)
plot_grid(px3,px4, nrow = 1)
dev.off()

### up gene 442做图五
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggpubr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

cgb_24vs0_deg <- read.csv("res_15102_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_deg <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  pull(ID)
cgb_24vs0_dep <- read.csv("res_15102_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  pull(ID)
COM_24vs0_dep <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  pull(ID)

setdiff_com_deg=setdiff(COM_24vs0_deg,cgb_24vs0_deg)
setdiff_com_dep=setdiff(COM_24vs0_dep,cgb_24vs0_dep)
df_1149= data.frame(ID=setdiff(setdiff_com_dep,setdiff_com_deg))

data <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  mutate(sig ="no")%>%right_join(df_1149,by="ID")

data$sig[data$pvalue <=0.05&data$log2FoldChange>=0 ] = "up"
data$sig[data$pvalue <=0.05&data$log2FoldChange<=0 ] = "down"

DEP=read.csv("res_T26_24hvs0h_DEP.csv",header = T)
dfcount00 = read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/AT_codon_count2.csv",header = T)%>%
  filter(!is.na(s2codon))
dfcount0 = read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/AT_codon_count2.csv",header = T)%>%
  right_join(DEP,by="ID")%>%filter(!is.na(s2codon))
dfcount = read.csv("C:/Users/zhengxueao/OneDrive/project/og_data/AT/AT_codon_count2.csv",header = T)%>%
  right_join(data,by="ID")%>%filter(!is.na(s2codon))%>%filter(sig=="up")

median(dfcount00$s2codon)
median(dfcount0$s2codon)
median(dfcount$s2codon)

library(ggpubr)
my_comparisons=list(c("no","1149"))
la_p1 = 50
p1 <- ggplot(dfcount,aes(x=sig,y=AAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "AAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p1

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","AAA"]))-median(c(df034[df034$sig == "down","AAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$AAA,N1,replace = FALSE))) - mean(c(sample(df034$AAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p11 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p11

la_p1 =30
p2 <- ggplot(dfcount,aes(x=sig,y=CAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "CAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p2

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","CAA"]))-median(c(df034[df034$sig == "down","CAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$CAA,N1,replace = FALSE))) - mean(c(sample(df034$CAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p12 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p12

la_p1 = 60
p3 <- ggplot(dfcount,aes(x=sig,y=GAA,fill=sig)) +
  geom_boxplot(outlier.colour = NA)+
  labs( x= "", y = "GAA",fill="") +
  theme_classic() +
  scale_fill_wsj('dem_rep') +
  scale_y_continuous(limits=c(0, la_p1))+
  guides(fill=FALSE)+
  stat_compare_means(aes(group = sig),comparisons = my_comparisons , method= "wilcox.test",label = "p.format",label.y = la_p1*2/3)
p3

df034 = dfcount
y_limit = 4
N1 = sum(df034$sig == "up")
N2 = sum(df034$sig == "down")
x = median(c(df034[df034$sig == "up","GAA"]))-median(c(df034[df034$sig == "down","GAA"]))
table(df034$sig)
sample_dovsup <- function(){
c1 = mean(c(sample(df034$GAA,N1,replace = FALSE))) - mean(c(sample(df034$GAA,N2,replace = FALSE)))
     return(c1)}
c2 = replicate(10000,sample_dovsup())
df001 <- data.frame(c2)
p13 <- ggplot(df001,aes(x=c2)) +
  geom_density()+
  labs( x= "", y = "Frenquency",title = x) +
  scale_fill_wsj('dem_rep')+
  #facet_grid(~Total_aa_scale,) +
  scale_x_continuous(limits=c(-y_limit,y_limit),breaks = c(-y_limit,-y_limit*3/4,-y_limit*1/2,-y_limit*1/4,0,y_limit*1/4,y_limit*1/2,y_limit*3/4,y_limit))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p13
library(cowplot)
pdf(file = "mRNA_A_ending.pdf",7,4)
plot_grid(p3,p2,p1, nrow = 1)
dev.off()




### 20220621
### 上调不上调

rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(dplyr)
library(ggpubr)
library(ggVennDiagram)
library(clusterProfiler)
setwd("D:/BigData/cgb/geom_col/og_data")

COM_24vs0_deg_up <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.584)%>%
  pull(ID)
COM_24vs0_deg_down <- read.csv("res_T26_24hvs0h_DEG.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.584)%>%
  pull(ID)
COM_24vs0_dep_up <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.23)%>%
  pull(ID)
COM_24vs0_dep_down <- read.csv("res_T26_24hvs0h_DEP.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.23)%>%
  pull(ID)
ratio_deg_up <- read.csv("res_ratio24to0_DEG.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.584)%>%
  pull(ID)
ratio_deg_down <- read.csv("res_ratio24to0_DEG.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.584)%>%
  pull(ID)
ratio_dep_up <- read.csv("res_ratio24to0_DEP.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > 0.23)%>%
  pull(ID)
ratio_dep_down <- read.csv("res_ratio24to0_DEP.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -0.23)%>%
  pull(ID)

deg_changeup_ratiodown=intersect(COM_24vs0_deg_up,ratio_deg_down)
deg_changeup_ratiodown_dif1=setdiff(COM_24vs0_deg_up,ratio_deg_down)
deg_changeup_ratiodown_dif2=setdiff(ratio_deg_down,COM_24vs0_deg_up)

l1= length(deg_changeup_ratiodown)
l2= length(deg_changeup_ratiodown_dif1)
l3= length(deg_changeup_ratiodown_dif2)

deg_changedown_ratioup=intersect(COM_24vs0_deg_down,ratio_deg_up)
dep_changeup_ratiodown=intersect(COM_24vs0_dep_up,ratio_dep_down)
dep_changedown_ratioup=intersect(COM_24vs0_dep_down,ratio_dep_up)

X = list(COM_24vs0_deg_up=COM_24vs0_deg_up,ratio_deg_down=ratio_deg_down)
p1=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
X = list(COM_24vs0_deg_down=COM_24vs0_deg_down,ratio_deg_up=ratio_deg_up)
p2=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
X = list(COM_24vs0_dep_up=COM_24vs0_dep_up,ratio_dep_down=ratio_dep_down)
p3=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
X = list(COM_24vs0_dep_down=COM_24vs0_dep_down,ratio_dep_up=ratio_dep_up)
p4=ggVennDiagram(X)+ scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")+ scale_color_brewer()
library(cowplot)
pdf(file = "Venn.pdf",15,3.5)
plot_grid(p1,p2,p3,p4, nrow = 1)
dev.off()

length(COM_24vs0_deg_up)

df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)

as.data.frame(deg_changeup_ratiodown)%>%rename(ID=deg_changeup_ratiodown)%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("TableS3.csv",row.names = F)
as.data.frame(dep_changeup_ratiodown)%>%rename(ID=dep_changeup_ratiodown)%>%
  left_join(df_anno,by=c("ID"="TAIR"))%>%
  select(-type)%>%
  write.csv("TableS4.csv",row.names = F)

###1
cgb_UG_go=deg_changeup_ratiodown%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.55,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)

df1 = as.data.frame(cgb_UG_go)%>%
  arrange(pvalue)%>%dplyr::slice(1:15)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "less upregulated"
write.csv(df1,"TableS5.csv")
# ###2
# cgb_DG_go=deg_changedown_ratioup %>%
#  enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
# simplify(
#   cutoff = 0.52,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang",
#   semData = NULL)
# df2 = as.data.frame(cgb_DG_go)%>%
#   arrange(-pvalue)%>%
#   filter(Description != "immune response")
# df2$Description=factor(df2$Description,levels=df2$Description)
# df2$sig = "less downregulated"
# df = rbind (df2,df1)
ggplot(df1,aes(x=-log(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size=""
      #y="Pathway name",
      # title="Pathway enrichment")
      # )+
    # scale_x_continuous(limits=c(1.5,12)
      )+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("cgb_DEG_go.pdf",width = 5.3,height =6)

###3
cgb_UP_go=dep_changeup_ratiodown%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.50,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)


df1 = as.data.frame(cgb_UP_go)%>%
  arrange(pvalue)%>%dplyr::slice(1:15)%>%arrange(-pvalue)
df1$Description=factor(df1$Description,levels=df1$Description)
df1$sig = "less upregulated"
write.csv(df1,"TableS6.csv")
# ###4
# cgb_DP_go=dep_changedown_ratioup%>%
#  enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
# simplify(
#   cutoff = 0.55,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang",
#   semData = NULL)
# df2 = as.data.frame(cgb_DP_go)%>%
#   arrange(-pvalue)
# df2$Description=factor(df2$Description,levels=df2$Description)
# df2$sig = "less downregulated"
# df = rbind (df2,df1)

ggplot(df1,aes(x=-log(pvalue),y=Description, fill =sig)) +
    geom_col(aes())+
    labs(
      size=""
      # x="Enrichment"
      #y="Pathway name",
      # title="Pathway enrichment")
    )+
    # scale_x_continuous(limits=c(1.5,12))+
    theme_classic()+
    theme(
      axis.text.y = element_text(size = rel(1)),
      axis.text.x = element_text(size=rel(1)),
      axis.title.y = element_blank())
ggsave("cgb_DEP_go.pdf",width = 5,height =6)

