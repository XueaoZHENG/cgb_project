# Title     : TODO
# Objective : TODO
# Created by: Guest0001
# Created on: 2021/9/2

### mRNA total normal
### "2021-09-04 08:31:01 CST"
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


### 所有差异分析
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
setwd("D:/BigData/cgb/20220625_new")
df_rna=read.csv("./og/mRNA_normal.csv",header = T)%>%
  mutate(ID = gsub("gene:", "",X))
df_pro = read.csv("./og/pro_normal.CSV",header = T)%>%
  mutate(ID = gsub("\\..*", "",Accession))%>%
  mutate_if(is.numeric, funs(1000000*.))
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

### DEG_0Hvs24H_cgb
library(DESeq2)
cts <- round(deg_df[,c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3")])
mutant <- factor(rep(c('0H','24H'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","24H","0H"))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.5)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.5)] <- "down"
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("DEG_0Hvs24H_cgb0.csv",row.names = F)
df%>%filter(sig != "no")%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  write.csv("DEG_0Hvs24H_cgb.csv",row.names = F)

### DEG_0Hvs24H_COM
library(DESeq2)
cts <- round(deg_df[,c("COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")])
mutant <- factor(rep(c('0H','24H'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","24H","0H"))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.5)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.5)] <- "down"
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("DEG_0Hvs24H_COM0.csv",row.names = F)
df%>%filter(sig != "no")%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  write.csv("DEG_0Hvs24H_COM.csv",row.names = F)

### DEP_0Hvs24H_cgb
library(DESeq2)
cts <- round(dep_df[,c("cgb_0_1","cgb_0_2","cgb_0_3","cgb_24_1","cgb_24_2","cgb_24_3")])
mutant <- factor(rep(c('0H','24H'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","24H","0H"))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.2)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.2)] <- "down"

df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("DEP_0Hvs24H_cgb0.csv",row.names = F)
df%>%filter(sig != "no")%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  write.csv("DEP_0Hvs24H_cgb.csv",row.names = F)

### DEP_0Hvs24H_COM
library(DESeq2)
cts <- round(dep_df[,c("COM_0_1","COM_0_2","COM_0_3","COM_24_1","COM_24_2","COM_24_3")])
mutant <- factor(rep(c('0H','24H'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant","24H","0H"))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.2)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.2)] <- "down"
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("DEP_0Hvs24H_COM0.csv",row.names = F)
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%filter(sig != "no")%>%
  write.csv("DEP_0Hvs24H_COM.csv",row.names = F)

### Ch_DEP
library(DESeq2)
dep_ratio_m = dep_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,COM_0=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  mutate(cgb_24V0_1=1000000*cgb_24_1/cgb_0,
              cgb_24V0_2=1000000*cgb_24_2/cgb_0,
              cgb_24V0_3=1000000*cgb_24_3/cgb_0,
              COM_24V0_1=1000000*COM_24_1/COM_0,
              COM_24V0_2=1000000*COM_24_2/COM_0,
              COM_24V0_3=1000000*COM_24_3/COM_0)%>%
  dplyr::select(ID,cgb_24V0_1,cgb_24V0_2,cgb_24V0_3,COM_24V0_1,COM_24V0_2,COM_24V0_3)
cts <- round(dep_ratio_m[,c("cgb_24V0_1","cgb_24V0_2","cgb_24V0_3","COM_24V0_1","COM_24V0_2","COM_24V0_3")])
mutant <- factor(rep(c('cgb','COM'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant",'cgb','COM'))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.2)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.2)] <- "down"
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("Ch_DEP0.csv",row.names = F)
df%>%filter(sig != "no")%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  write.csv("Ch_DEP.csv",row.names = F)

### Ch_DEG
library(DESeq2)
deg_ratio_m = deg_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,COM_0=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  mutate(cgb_24V0_1=1000000*cgb_24_1/cgb_0,
              cgb_24V0_2=1000000*cgb_24_2/cgb_0,
              cgb_24V0_3=1000000*cgb_24_3/cgb_0,
              COM_24V0_1=1000000*COM_24_1/COM_0,
              COM_24V0_2=1000000*COM_24_2/COM_0,
              COM_24V0_3=1000000*COM_24_3/COM_0)%>%
  dplyr::select(ID,cgb_24V0_1,cgb_24V0_2,cgb_24V0_3,COM_24V0_1,COM_24V0_2,COM_24V0_3)%>%
  filter(cgb_24V0_1 != 0,cgb_24V0_1 != Inf,
           cgb_24V0_2 != 0,cgb_24V0_2 != Inf,
           cgb_24V0_3 != 0,cgb_24V0_3 != Inf,
           COM_24V0_1 != 0,COM_24V0_1 != Inf,
           COM_24V0_2 != 0,COM_24V0_2 != Inf,
           COM_24V0_3 != 0,COM_24V0_3 != Inf)
cts <- round(deg_ratio_m[,c("cgb_24V0_1","cgb_24V0_2","cgb_24V0_3","COM_24V0_1","COM_24V0_2","COM_24V0_3")])
mutant <- factor(rep(c('cgb','COM'), each=3))
colData <- data.frame(row.names=colnames(cts),mutant)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design= ~ mutant)
dds <- DESeq(dds)
res01 <- results(dds, contrast=c("mutant",'cgb','COM'))
df <-as.data.frame(res01)%>%
  filter(pvalue != "NA")%>%
  filter(pvalue <0.05)%>%
  mutate(sig = "no")
df$ID=row.names(df)
df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df$sig[df$pvalue < 0.05 & df$log2FoldChange >= log2(1.5)] <- "up"
df$sig[df$pvalue < 0.05 & df$log2FoldChange <= -log2(1.5)] <- "down"
df%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("Ch_DEG0.csv",row.names = F)
df%>%filter(sig != "no")%>%left_join(df_anno,by=c("ID"="TAIR"))%>%
  write.csv("Ch_DEG.csv",row.names = F)


### TE 分析
dep_ratio_m = dep_df%>%
  mutate(cgb_0=(cgb_0_1+cgb_0_2+cgb_0_3)/3,COM_0=(COM_0_1+COM_0_2+COM_0_3)/3)%>%
  mutate(cgb_24V0_1=1000000*cgb_24_1/cgb_0,
              cgb_24V0_2=1000000*cgb_24_2/cgb_0,
              cgb_24V0_3=1000000*cgb_24_3/cgb_0,
              COM_24V0_1=1000000*COM_24_1/COM_0,
              COM_24V0_2=1000000*COM_24_2/COM_0,
              COM_24V0_3=1000000*COM_24_3/COM_0)%>%
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
  mutate(te_cgb_1=cgb_24V0_1/(rna_cgb_ratio*100),
              te_cgb_2=cgb_24V0_2/(rna_cgb_ratio*100),
              te_cgb_3=cgb_24V0_3/(rna_cgb_ratio*100),
              te_COM_1=COM_24V0_1/(rna_COM_ratio*100),
              te_COM_2=COM_24V0_2/(rna_COM_ratio*100),
              te_COM_3=COM_24V0_3/(rna_COM_ratio*100))%>%
  dplyr::select(ID,te_cgb_1,te_cgb_2,te_cgb_3,
         te_COM_1,te_COM_2,te_COM_3)%>%
  filter(te_cgb_1 != "NA",te_cgb_1 !=Inf,te_cgb_2 != "NA",te_cgb_2 !=Inf,te_cgb_3 != "NA",te_cgb_3 !=Inf,
         te_COM_1 != "NA",te_COM_1 !=Inf,te_COM_2 != "NA",te_COM_2 !=Inf,te_COM_3 != "NA",te_COM_3 !=Inf)
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
write.csv(as.data.frame(res01), file='TE_cgbvsCOM.csv')

### TE down GO 图形生成
rm(list=ls())
setwd("D:/BigData/cgb/20220625_new")
library(clusterProfiler)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)
library(Rmisc)
df_ratio_p = read.csv("TE_cgbvsCOM.csv",header = T)
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= log2(1.2)] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -log2(1.2)] <- "down"

df_anno=read.csv("D:/BigData/arabidopsis_go/Arabidopsis_anno_zxa_v1.csv",header = T)
df_ratio_p%>%filter(sig == "down")%>%
  left_join(df_anno,by=c("X"="TAIR"))%>%
  dplyr::select(-type)%>%
  write.csv("TE_down_geneTableS7.csv",row.names = F)

cgb_UP_go=df_ratio_p%>%filter(sig == "down")%>%pull(X)%>%
 enrichGO(OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")%>%
simplify(
  cutoff = 0.6,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL)
write.csv(cgb_UP_go,"TE_down_tableS8.csv")

df1 = read.csv("TE_down_tableS7.csv",row.names = 1)%>%
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

###TE 热图
rm(list=ls())
library(dplyr)
library(tidyr)
library(pheatmap)
setwd("D:/BigData/cgb/20220625_new")
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
  t()%>%
  as.data.frame()%>%
  mutate(ID=row.names(dep_ratio_m))
colnames(dep_ratio_m1)=c(colnames(dep_ratio_m),"ID")

df_ratio_p = read.csv("TE_cgbvsCOM.csv",header = T)
df_ratio_p$sig <- "no"  ### 筛选并填充数据
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange >= log2(1.2)] <- "up"
df_ratio_p$sig[df_ratio_p$pvalue <= 0.05 & df_ratio_p$log2FoldChange <=  -log2(1.2)] <- "down"

df=df_ratio_p%>%filter(sig == "down")%>%
  left_join(dep_ratio_m1,by=c("X"="ID"))%>%
  dplyr::select(c(12:14),c(9:11))

pheatmap(df,
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
         width = 5,
         height = 6)


### 所有的差异因数量的描述性统计 饼图
rm(list=ls())
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
setwd("D:/BigData/cgb/20220625_new")

temp=list.files(path="D:/BigData/cgb/20220625_new",pattern="*0.csv")
name_pdf=paste(temp,".pdf",sep="")
data_final_DEG= data.frame(sig=list(),count=list(),file=list())
for (fil in temp) {
data = read.csv(fil,header = T)
data$sig =factor(data$sig, levels = c("down","up","no"))
data = data%>%filter(log2FoldChange !="NA")%>%
    group_by(sig) %>%
    summarize(count=n())%>%
    mutate(file = fil)%>%
    mutate(perc = count / sum(count)) %>%
    mutate(labels = scales::percent(perc))
  data_final_DEG = rbind(data_final_DEG,data)
}

data_final = data_final_DEG
temp_0=list.files(path="D:/BigData/cgb/geom_col/og_data",pattern="*0.csv")
data_final$file=factor(data_final$file,
                       levels = c("DEG_0Hvs24H_cgb0.csv",
  "DEG_0Hvs24H_COM0.csv",
  "Ch_DEG0.csv",
  "DEP_0Hvs24H_cgb0.csv",
  "DEP_0Hvs24H_COM0.csv",
  "Ch_DEP0.csv"))

### 饼图
ggplot(data_final,aes(x="",y=perc, fill =sig))+
  facet_grid(~file,)+
  geom_col() +
    scale_fill_wsj('dem_rep')+
    coord_polar(theta = "y")+
    geom_text(aes(label = paste0(round(perc * 100, 0), "%", count)),
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
setwd("D:/BigData/cgb/20220625_new")

COM_24vs0_deg_up <- read.csv("DEG_0Hvs24H_COM0.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > log2(1.5))%>%
  pull(ID)
COM_24vs0_deg_down <- read.csv("DEG_0Hvs24H_COM0.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -log2(1.5))%>%
  pull(ID)
COM_24vs0_dep_up <- read.csv("DEP_0Hvs24H_COM0.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > log2(1.2))%>%
  pull(ID)
COM_24vs0_dep_down <- read.csv("DEP_0Hvs24H_COM0.csv",header = T)%>%
  filter(!is.na(pvalue))%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -log2(1.2))%>%
  pull(ID)

ratio_deg_up <- read.csv("Ch_DEG0.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > log2(1.5))%>%
  pull(ID)
ratio_deg_down <- read.csv("Ch_DEG0.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -log2(1.5))%>%
  pull(ID)
ratio_dep_up <- read.csv("Ch_DEP0.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange > log2(1.2))%>%
  pull(ID)
ratio_dep_down <- read.csv("Ch_DEP0.csv",header = T)%>%
  filter(pvalue <0.05)%>%
  filter(log2FoldChange < -log2(1.2))%>%
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

