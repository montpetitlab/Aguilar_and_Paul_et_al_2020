library(pheatmap)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(RColorBrewer)

# heatmap -----------------------------------------------------------------

samples <- read_csv("RNA_Seq/inputs/ribo-samples.csv")
files <- files <- list.files("RNA_Seq/outputs/deseq2", ".csv$", full.names = T)
files <- files[!grepl(pattern = "*sig*", x = files)]
res <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("RNA_Seq\\/outputs\\/deseq2\\/", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  filter(padj < .01) %>%
  select(source, X1, log2FoldChange) 

# pivot wide so each experiment has a col of log2fold changes
res_wide <- pivot_wider(data = res, id_cols = X1, names_from = source, values_from = log2FoldChange)
# remove rows with NAs
res_wide <- res_wide %>%
  filter(!is.na(res_srm1)) %>%
  filter(!is.na(res_dis3)) %>%
  filter(!is.na(res_csl4ph)) %>%
  filter(!is.na(res_enp1)) %>%
  filter(!is.na(res_rrp6))
# make gene names rownames
res_wide <- as.data.frame(res_wide)
rownames(res_wide) <- res_wide$X1
res_wide <- res_wide[ , -1]
res_wide <- as.matrix(res_wide)

# set a label for row to label transcript type
mat_row <- data.frame(transcript = rownames(res_wide))
mat_row$transcript <- ifelse(grepl(pattern = "mRNA", x = mat_row$transcript), 
                             "mRNA", "pervasive_transcript")
rownames(mat_row) <- rownames(res_wide)
# set colors
pal <- brewer.pal(9, "Set1")
mat_colors <- list(transcript = pal[c(2, 9)])
names(mat_colors$transcript) <- unique(mat_row$transcript)

pdf(snakemake@output[["plot"]], width = 6, height = 6)
pheatmap(mat                  = res_wide,
         show_rownames        = F,
         show_colnames        = T,
         cluster_rows         = T,
         cluster_cols         = T,
         annotation_row       = mat_row,
         annotation_colors    = mat_colors,
         annotation_names_row = F)
dev.off()