library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# read in differential expression results
files <- list.files("outputs/deseq2", full.names = T)
files <- files[!grepl(pattern = "sig\\.csv$", x = files)]

log2fc <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub(".csv", "", source))

# reset levels so control comes first
log2fc$source <- factor(log2fc$source, 
                        levels = c("wt", "csl4", "enp1"))


# make plots --------------------------------------------------------------

cut <- ggplot(log2fc %>% filter(grepl("CUT", X1)), 
       aes(x = source, y = log2FoldChange)) +
  geom_violin() + 
  theme_minimal() + 
  labs(title = "CUTs", x = "", y = "Log2FC")

nut <-ggplot(log2fc %>% filter(grepl("NUT", X1)), 
       aes(x = source, y = log2FoldChange)) +
  geom_violin() + 
  theme_minimal() + 
  labs(title = "NUTs", x = "", y = "Log2FC")

sut <- ggplot(log2fc %>% filter(grepl("ST", X1)), 
       aes(x = source, y = log2FoldChange)) +
  geom_violin() + 
  theme_minimal() + 
  labs(title = "SUTs", x = "", y = "Log2FC")

# XUTs are not named with "XUT" in the name; create plot using 
# names in original bed file. 
xuts <- read_tsv("inputs/pervasive_transcripts/XUT.sorted.bed", 
                 col_names = c("chr", "start", "end", "name", "class", "strand"))
log2fc <- separate(data = log2fc, col = X1, into = c("name", "position"),
                   sep = "::", remove = T)
xut <- ggplot(log2fc %>% filter(name %in% xuts$name), 
       aes(x = source, y = log2FoldChange)) +
  geom_violin() + 
  theme_minimal() + 
  labs(title = "XUTs", x = "", y = "Log2FC")

pdf(snakemake@output[['fig']], height = 6.6, width = 8)
ggarrange(cut, nut, sut, xut, ncol = 2, nrow = 2)
dev.off()

# anova -------------------------------------------------------------------

sut_log2fc <- log2fc %>% 
  filter(grepl("ST", name)) %>%
  filter(!is.na(log2FoldChange))

sut_aov <- aov(log2FoldChange ~ source, data = sut_log2fc)
summary(sut_aov)

nut_log2fc <- log2fc %>% 
  filter(grepl("NUT", name)) %>%
  filter(!is.na(log2FoldChange))

nut_aov <- aov(log2FoldChange ~ source, data = nut_log2fc)
summary(nut_aov)

cut_log2fc <- log2fc %>% 
  filter(grepl("CUT", name)) %>%
  filter(!is.na(log2FoldChange))

cut_aov <- aov(log2FoldChange ~ source, data = cut_log2fc)
summary(cut_aov)

xut_log2fc <- log2fc %>% 
  filter(name %in% xuts$name) %>%
  filter(!is.na(log2FoldChange))

xut_aov <- aov(log2FoldChange ~ source, data = xut_log2fc)
summary(xut_aov)
