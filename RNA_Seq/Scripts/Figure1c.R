#Code for making heatmap of correlation co-effient
library(data.table)
library(corrplot)
library(dplyr)
setwd("~/Desktop/Paul_et_al_2019/")
fc_df_rb <- read.csv("./results/log2FC_all_rb.tsv", sep='\t')
log2fc <- fc_df_rb[complete.cases(fc_df_rb), ]
df <- log2fc[, c("rrp6", "dis3", "csl4","enp1","srm1")]
M <- cor(df)
png("./Figures/corr_plot_logfc_rb.png")
corrplot(M, method = "color", type="upper", addCoef.col = "white")
dev.off()
