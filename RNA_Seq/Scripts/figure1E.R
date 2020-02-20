library(ggplot2)
setwd("~/Desktop/Paul_et_al_2019/")
fc_df_rb <- read.csv("./results/log2FC_all_rb.tsv", sep='\t')
rna_class <- read.csv("./annotation/rna3.csv", sep = "\t")
xut <- read.csv("./annotation/xut_list.txt", header = F, sep = ";")
colnames(xut) <- c("gene_id", "RNAclass")
rna_class <- rbind(rna_class, xut)
df <- merge(fc_df_rb, rna_class, by = 'gene_id', all=FALSE)
df <- df[apply(df!=0, 1, all),]
df <- df[,-1]
df1 <- melt(df)
png("./Figures/figure1e.png")
ggplot(df1, aes(x=value, color=RNAclass)) + geom_line(stat="density", size = 1.0) + facet_wrap(~variable) + xlim(-4, 8) + theme_bw()
dev.off()
