library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggdendro)
library(ggpubr)
setwd("~/Paul_et_al_2019/")

fc_df_rb <- read.csv("./RNA_Seq/results/log2FC_all_rb.tsv", sep='\t')
transcript = fread('./RNA_Seq/annotation/yeast_all_annot.gtf')

#----------------------------------------------------------------------------------
type_rb = c('rrp6', 'dis3', 'csl4', 'enp1', 'srm1')
data_rb = fc_df_rb[,1:6]
colnames(data_rb) <- append(c('gene_id'), type_rb)

data_rb <- as.data.frame(data_rb)
biotype = c()
slc_row = c()
for (row in 1:nrow(data_rb)){
    current_id <- data_rb[row, 'gene_id']
    to_add = transcript[transcript$gene_id == current_id]$gene_biotype
    gene = transcript[transcript$gene_id == current_id]$gene_id
    biotype <- c(biotype, to_add)
    slc_row <- c(slc_row, gene)
}

data_rb <- data_rb[data_rb$gene_id %in% slc_row,]
   

data_rb$biotype <- biotype
data_rb <- data_rb[,c(1,6,5,2,3,4,7)]
#---------------------------------------------------------

ha = HeatmapAnnotation(df = data.frame(type = type_rb))

tmp <- data_rb[data_rb$biotype %in% c('pervasive', 'protein_coding'),]

colnames(tmp) <- c('gene_id', 'rrp6', 'dis3', 'csl4', 'enp1', 'srm1', 'RNAclass')

tmp$RNAclass[tmp$RNAclass == 'pervasive'] <- "Pervasive"
tmp$RNAclass[tmp$RNAclass == 'protein_coding'] <- "mRNA"
to_plot = as.matrix(tmp[,2:6])
rownames(to_plot) <- tmp$gene_id

kclus <- kmeans(to_plot, 6)
split <- paste0("C", kclus$cluster)

ht = Heatmap(to_plot, split=split, name = 'log2FC', gap = unit(2, "mm"),
       col = colorRamp2(c(-5, 0, 5), c('blue', '#FFFFFF', 'red')), width = unit(6,'cm'),
       top_annotation_height = unit(10, "mm"),row_title='Kmean clusters',
       show_row_names = FALSE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE)

ha = Heatmap(tmp$RNAclass, name = "", show_row_names = FALSE,
        width = unit(5, "mm"), col = c('#808080','#70B2FF'))

ht_list = ht + ha
png("./Figures/figure1d.png")
draw(ht_list, row_title='Kmean clusters')
dev.off()

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(to_plot, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
