library(ggplot2)
library(tidyverse)
fc_df_rb <- read.csv("./results/log2FC_all_rb.tsv", sep='\t')
rna_class <- read.csv("./annotation/tollarvy_rna_classes2.csv", sep = "\t")
colnames(rna_class) <- c("gene_id", "RNAclass")
df <- merge(fc_df_rb, rna_class, by = 'gene_id', all=FALSE)

df <- df[,-1]
df1 <- melt(df)
png("./Figures/figure2a.png")
ggplot(df1, aes(x = RNAclass, y = value, fill=variable)) + facet_wrap( ~ variable) + xlab("Classes of mRNA") +
  ylab("log2FC") + ylim(-6,6) + geom_boxplot(outlier.size = 0.1, notch=TRUE) + theme(panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.7, color = "black"),
                                                                                     axis.text=element_text(size=8), axis.title=element_text(size=14,face="bold"), 
                                                                                     strip.text.x = element_text(size = 15, face ="bold.italic"), legend.position = c(0.8, 0.2))
dev.off()
