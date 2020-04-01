library(pheatmap)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)

# heatmap -----------------------------------------------------------------

files <- unlist(snakemake@input)
res <- files[1:5] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  filter(padj < .05) %>%
  select(source, X1, log2FoldChange) %>%
  filter(!grepl("sn", X1)) %>%
  filter(!grepl("RDN", X1))

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

pdf(snakemake@output[["heatmap"]], width = 6, height = 6)
# reset values to a min and max; needs to be paired with a heatmap
# volcano plot like seen below
tmp <- ifelse(res_wide > 7, 7, res_wide)
tmp <- ifelse(tmp < -7, -7, tmp)
pheatmap(mat                  = tmp,
         show_rownames        = F,
         show_colnames        = T,
         cluster_rows         = T,
         cluster_cols         = T,
         annotation_row       = mat_row,
         annotation_colors    = mat_colors,
         annotation_names_row = F)
dev.off()

res_wide2 <- as.data.frame(res_wide)
res_wide2$transcripts <- rownames(res_wide2)
mat_row2 <- mat_row
mat_row2$transcripts <- rownames(mat_row2)

res_long <- pivot_longer(res_wide2, cols = -transcripts, 
                         names_to = "sample", values_to = "log2FC")
res_long <- left_join(res_long, mat_row2)

pdf(snakemake@output[['density']], height = 5, width = 8)
ggplot(res_long, aes(x = log2FC, fill = transcript, color = transcript)) +
  geom_density(alpha = .7) +
  theme_minimal() +
  facet_wrap(~sample) +
  scale_fill_manual(values = c(pervasive_transcript = "#377EB8",
                               mRNA = "#999999")) +
  scale_color_manual(values = c(pervasive_transcript = "#377EB8",
                               mRNA = "#999999")) +
  theme(legend.position = "bottom")
dev.off()