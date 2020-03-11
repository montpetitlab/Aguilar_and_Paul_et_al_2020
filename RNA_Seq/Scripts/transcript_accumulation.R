library(dplyr)
library(readr)
library(purrr)
library(ggpubr)


# heatmap -----------------------------------------------------------------

samples <- read_csv("inputs/ribo-samples.csv")
files <- files <- list.files("outputs/deseq2", ".csv$", full.names = T)
files <- files[!grepl(pattern = "*sig*", x = files)]
res <- files %>%
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
## sn* == snNRA
## *_mRNA == mRNA
## ..-... == XUT

res$class <- ifelse(grepl(pattern = "mRNA", x = res$class),
                    "mRNA", res$class)
res$class <- ifelse(grepl(pattern = "NME1", x = res$class),
                    "mRNA", res$class)
res$class <- ifelse(grepl(pattern = "ST", x = res$class),
                    "SUT", res$class)
res$class <- ifelse(grepl(pattern = "sn", x = res$class),
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

sut <- pervasive %>% filter(class == "SUT")
sut_plt <- ggplot(sut %>% 
         group_by(source) %>% 
         arrange(log2FoldChange) %>% 
         mutate(rn = row_number())) + 
  geom_line(aes(x=log2FoldChange, y=rn/length(unique(rn)), color=source)) +
  theme_classic() + 
  labs(x = "", y = "", title = "SUTs")

nut <- pervasive %>% filter(class == "NUT")
nut_plt <- ggplot(nut %>% 
                    group_by(source) %>% 
                    arrange(log2FoldChange) %>% 
                    mutate(rn = row_number())) + 
  geom_line(aes(x=log2FoldChange, y=rn/length(unique(rn)), color=source)) +
  theme_classic() + 
  labs(x = "", y = "", title = "NUTs") +
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank())

xut <- pervasive %>% filter(class == "XUT")
xut_plt <- ggplot(xut %>% 
                    group_by(source) %>% 
                    arrange(log2FoldChange) %>% 
                    mutate(rn = row_number())) + 
  geom_line(aes(x=log2FoldChange, y=rn/length(unique(rn)), color=source)) +
  theme_classic() + 
  labs(x = "", y = "", title = "XUTs") +
  theme(axis.text.y=element_blank())

cut <- pervasive %>% filter(class == "CUT")
cut_plt <- ggplot(cut %>% 
                    group_by(source) %>% 
                    arrange(log2FoldChange) %>% 
                    mutate(rn = row_number())) + 
  geom_line(aes(x=log2FoldChange, y=rn/length(unique(rn)), color=source)) +
  theme_classic() + 
  labs(x = "", y = "", title = "CUTs") +
  theme(axis.text.x=element_blank())


plt <- ggarrange(cut_plt, nut_plt, sut_plt, xut_plt, common.legend = T, 
          legend = "bottom", ncol = 2, nrow = 2)

# Annotate the figure by adding a common labels
annotate_figure(plt,
                bottom = text_grob("log2FC (mutant vs. control @ 37C)", size = 14),
                left = text_grob("Total Fraction of Transcripts in Class", rot = 90),
)
