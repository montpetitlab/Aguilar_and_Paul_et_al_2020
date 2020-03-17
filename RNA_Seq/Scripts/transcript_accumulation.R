library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(tidyr)

samples <- snakemake@input[["samples"]]
files <- unlist(snakemake@input)
res <- files[1:5] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  filter(padj < .01) %>%
  select(source, X1, log2FoldChange) 

## define transcript class based on transcript name
res <- separate(data = res, col = X1, into = c("class", "coord"), 
                sep = "::", remove = F)
## edit class
## CUT == CUT
## NUT == NUT
## ST  == SUT 
## *_mRNA == mRNA
## ..-... == XUT

res$class <- ifelse(grepl(pattern = "mRNA", x = res$class),
                    "mRNA", res$class)
res$class <- ifelse(grepl(pattern = "NME1", x = res$class),
                    "mRNA", res$class)
res$class <- ifelse(grepl(pattern = "ST", x = res$class),
                    "SUT", res$class)
snRNA <- c("LSR1", "snR14", "snR7-L", "snR7-S", "snR6", "snR19")
res$class <- ifelse(grepl(pattern = "sn", x = res$class),
                           ifelse(!res$class %in% snRNA,
                                  "snoRNA", "snRNA"), res$class)
res$class <- ifelse(grepl(pattern = "LSR1", x = res$class),
                           "snRNA", res$class)
res$class <- ifelse(grepl(pattern = "[FR]-", x = res$class),
                    "XUT", res$class)
res$class <- ifelse(grepl(pattern = "NUT", x = res$class),
                    "NUT", res$class)
res$class <- ifelse(grepl(pattern = "CUT", x = res$class),
                    "CUT", res$class)

# filter to just pervasive transcripts
pervasive <- res %>%
  filter(class %in% c("NUT", "CUT", "SUT", "XUT"))

pervasive <- pervasive %>% 
  group_by(source, class) %>% 
  arrange(log2FoldChange) %>% 
  mutate(rn = row_number()) %>%
  mutate(prop = rn/length(unique(rn)))

pdf(snakemake@output[["accum_plot"]], height = 6, width = 7)
ggplot(pervasive) + 
  geom_line(aes(x=log2FoldChange, y= prop, color=source)) +
  theme_classic() + 
  facet_wrap(~class) + 
  xlim(c(-8, 10)) +
  theme(legend.position = "bottom") +
  labs(y = "Total Fraction of Transcripts in Class",
       x = "log2FC (mutant vs. control @ 37C)")
dev.off()

