library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

files <- unlist(snakemake@input)

# files <- c("outputs/deseq2/res_wt_bm766_v_wt.csv", "outputs/deseq2/res_csl4_bm5_v_wt.csv",
#            "outputs/deseq2/res_csl4_bm766_v_wt.csv", "outputs/deseq2/res_enp1_bm5_v_wt.csv",
#            "outputs/deseq2/res_enp1_bm766_v_wt.csv")

res <- files[1:5] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("_v_wt.csv", "", source)) %>%
  dplyr::select(source, X1, log2FoldChange, padj) 

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


pervasive$plasmid <- ifelse(grepl("bm5", pervasive$source), "bm5", "bm766")
pervasive$strain <- ifelse(grepl("wt", pervasive$source), "wt", NA)
pervasive$strain <- ifelse(grepl("enp1", pervasive$source), "enp1", pervasive$strain)
pervasive$strain <- ifelse(grepl("csl4", pervasive$source), "csl4", pervasive$strain)


pdf(snakemake@output[["fig"]], height = 6, width = 7)
ggplot(pervasive) + 
  geom_line(aes(x=log2FoldChange, y= prop, color=strain, linetype = plasmid), size = 1.1, alpha = .7) +
  theme_minimal() + 
  facet_wrap(~class, scales = "free") + 
  xlim(c(-5, 10)) +
  theme(legend.position = "bottom", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Total Fraction of Transcripts in Class",
       x = "log2FC (mutant vs. control with bm5)") +
  scale_color_manual(values = c(csl4 = "#47a4e4", wt = "black", 
                                enp1 = "#129060"))

dev.off()
