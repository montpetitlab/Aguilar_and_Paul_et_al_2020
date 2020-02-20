library(ggplot2)
library(reshape2)
library(data.table)
library(cluster)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(plyr)
setwd("~/Paul_et_al_2019/Tag_seq/")
WT_pBM5vscsl4ph_pBM5 <- read.csv("output_edgeR/WT_pBM5vscsl4ph_pBM5.tsv", sep = "\t")[,c(6,1)]
WT_pBM5vscsl4ph_pBM766 <- read.csv("output_edgeR/WT_pBM5vscsl4ph_pBM766.tsv", sep = "\t")[,c(6,1)]
WT_pBM5vsenp1_1_pBM5 <- read.csv("output_edgeR/WT_pBM5vsenp1_1_pBM5.tsv", sep = "\t")[,c(6,1)]
WT_pBM5vsenp1_1_pBM766 <- read.csv("output_edgeR/WT_pBM5vsenp1_1_pBM766.tsv", sep = "\t")[,c(6,1)]
WT_pBM5vsWT_pBM766 <- read.csv("output_edgeR/WT_pBM5vsWT_pBM766.tsv", sep = "\t")[,c(6,1)]
WT_pBM766vscsl4ph_pBM766 <- read.csv("output_edgeR/WT_pBM766vscsl4ph_pBM766.tsv", sep = "\t")[,c(6,1)]
WT_pBM766vsenp1_1_pBM766 <- read.csv("output_edgeR/WT_pBM766vsenp1_1_pBM766.tsv", sep = "\t")[,c(6,1)]
WT_pBM5vsWT_pBM766 <- read.csv("output_edgeR/WT_pBM5vsWT_pBM766.tsv", sep = "\t")[,c(6,1)]

temp <- merge(WT_pBM5vscsl4ph_pBM5, WT_pBM5vscsl4ph_pBM766, by="mrna", all = T)
temp <- merge(temp, WT_pBM5vsenp1_1_pBM5, by="mrna", all = T)
temp <- merge(temp, WT_pBM5vsenp1_1_pBM766, by="mrna", all = T)
temp <- merge(temp, WT_pBM5vsWT_pBM766, by="mrna", all = T)
colnames(temp) <- c("RNA", "WT_pBM5vscsl4ph_pBM5", "WT_pBM5vscsl4ph_pBM766", "WT_pBM5vsenp1_1_pBM5", "WT_pBM5vsenp1_1_pBM766", "WT_pBM5vsWT_pBM766")
rna_class <- read.csv("~/Downloads/files/mRNA_list.txt", sep = "\t", header = F)[,c(1,3)]
colnames(rna_class) <- c("RNA", "RNA_Class")
df <- merge(temp, rna_class, by="RNA", all=F)

temp <- merge(temp, WT_pBM766vscsl4ph_pBM766, by="mrna", all = T)
temp <- merge(temp, WT_pBM766vsenp1_1_pBM766, by="mrna", all = T)
temp <- merge(temp, WT_pBM5vsWT_pBM766, by="mrna", all = T)
temp <- merge(temp, csl4ph_pBM5vscsl4ph_pBM766, by="mrna", all = T)
temp <- merge(temp, enp1_1_pBM5vsenp1_1_pBM766, by="mrna", all = T)
colnames(temp) <- c("RNA", "WT_pBM5vscsl4ph_pBM5", "WT_pBM5vsenp1_1_pBM5", "WT_pBM766vscsl4ph_pBM766", "WT_pBM766vsenp1_1_pBM766", "WT_pBM5vsWT_pBM766", "csl4ph_pBM5vscsl4ph_pBM766", "enp1_1_pBM5vsenp1_1_pBM766")
temp <- temp[complete.cases(temp),]
rna_class <- read.csv("~/Desktop/rna3.txt", sep="\t")
colnames(rna_class) <- c("RNA", "RNA_Class")
levels(rna_class$RNA_Class)[levels(rna_class$RNA_Class) %in%  c("CUT","SUT","XUT", 'NUT')] <- "Pervasive transcript"
log2fc <- merge(temp, rna_class, by='RNA', all=FALSE)
df <- as.matrix(log2fc[,c(2:6)])
hc <- hclust(df)
pa = pam(df, k = 6)
png(file="~/Desktop/log2fc_dt.png")
Heatmap(df, name = "log2FC", split = paste0("C", pa$clustering), gap = unit(1, "mm"),row_title_gp = gpar(font = 1:2), 
        width = unit(30, "mm"), show_row_names = FALSE, column_names_gp = gpar(fontsize = 8), col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
        cluster_columns = F, show_row_dend = F) + Heatmap(log2fc$RNA_Class, name = "RNA Class", width = unit(5, "mm"))
dev.off()
Heatmap(log2fc$RNA_Class, name = "RNA Class", width = unit(5, "mm"))

df <- melt(temp)
CUT <- temp[grep("CUT", temp$RNA),]
df <- melt(CUT)
p0 = ggplot(df, aes(value)) + stat_ecdf(aes(color = variable)) +scale_color_viridis_d()+ xlim(-4,4)
p0 = ggplot(df, aes(variable, value)) + geom_violin() + geom_boxplot(width=0.05) + theme(axis.text.x = element_text(angle = 45)) + labs(title = "CUTs", x="Strains", y="log2FC") 
SUT <- temp[grep("SUT", temp$RNA),]
df <- melt(SUT)
p0 = ggplot(df, aes(value)) + stat_ecdf(aes(color = variable)) +scale_color_viridis_d()+ xlim(-5,5)
p0 = ggplot(df, aes(variable, value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45)) + labs(title = "SUTs", x="Strains", y="log2FC")

SUT <- temp[grep("NUT", temp$RNA),]
df <- melt(SUT)
p0 = ggplot(df, aes(value)) + stat_ecdf(aes(color = variable)) +scale_color_viridis_d()+ xlim(-5,5)
p0 = ggplot(df, aes(variable, value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45)) + labs(title = "NUTs", x="Strains", y="log2FC")

SUT <- temp[grep("XUT", temp$RNA),]
df <- melt(SUT)
p0 = ggplot(df, aes(value)) + stat_ecdf(aes(color = variable)) +scale_color_viridis_d()+ xlim(-5,5)
p0 = ggplot(df, aes(variable, value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45)) + labs(title = "XUTs", x="Strains", y="log2FC")


RNA_Class <- read.csv("~/DATA2/Tag_Seq_Analysis/snoRNA.txt", sep = "\t")
colnames(RNA_Class) <- c("RNA", "RNA_Class")
temp1 <- merge(temp, RNA_Class, by="RNA", all = F)

enp1 <- temp1[,c(1, 2, 3, 4, 5, 6)]
df <- melt(enp1)
p0 = ggplot(df, aes(value)) + stat_ecdf(aes(color = variable)) +scale_color_viridis_d()+ xlim(-2,5) + xlab("log2FC") + ylab("Percentage") +
  labs(title = "snoRNA")
csl4 <- temp1[,c(1,2,4,7)]
df <- melt(csl4)
p0 = ggplot(df, aes(x = value)) + 
  geom_density(aes(group = variable, colour = variable)) + facet_wrap(~RNA_Class)
df1 <- melt(df)
ggplot(df1, aes(x = variable, y = value)) + 
  geom_boxplot() + facet_wrap(~RNA_Class) 
