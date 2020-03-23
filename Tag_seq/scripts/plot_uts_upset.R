library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(UpSetR)

# read in differential expression results
files <- unlist(snakemake@input)
#files <- c('outputs/deseq2/res_wt.csv', 'outputs/deseq2/res_csl4.csv', 
#           'outputs/deseq2/res_enp1.csv')
log2fc <- files[1:3] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub(".csv", "", source)) %>%
  filter(padj < .05)

# upset plot --------------------------------------------------------------

# filter to pervasive transcripts and make a vector that indicates transcript type

log2fc$pervasive <- log2fc$X1
log2fc$pervasive <- ifelse(grepl(pattern = "ST", x = log2fc$pervasive),
                           "pervasive", log2fc$pervasive)
log2fc$pervasive <- ifelse(grepl(pattern = "[FR]-", x = log2fc$pervasive),
                           "pervasive", log2fc$pervasive)
log2fc$pervasive <- ifelse(grepl(pattern = "NUT", x = log2fc$pervasive),
                           "pervasive", log2fc$pervasive)
log2fc$pervasive <- ifelse(grepl(pattern = "CUT", x = log2fc$pervasive),
                           "pervasive", log2fc$pervasive)

log2fc <- log2fc %>% 
  filter(pervasive == "pervasive") %>% 
  filter(log2FoldChange < 0)

csl4_down <- log2fc %>%
  filter(source == "csl4") %>%
  filter(log2FoldChange < 0) %>%
  select(X1)
csl4_down_vec <- as.character(csl4_down$X1)
names(csl4_down_vec) <- csl4_down$X1

wt_down <- log2fc %>%
  filter(source == "wt") %>%
  filter(log2FoldChange < 0) %>%
  select(X1)
wt_down_vec <- as.character(wt_down$X1)
names(wt_down_vec) <- wt_down$X1

enp1_down <- log2fc %>%
  filter(source == "enp1") %>%
  filter(log2FoldChange < 0) %>%
  select(X1)
enp1_down_vec <- as.character(enp1_down$X1)
names(enp1_down_vec) <- enp1_down$X1

down <- fromList(list(wt = wt_down_vec,
                      enp1 = enp1_down_vec,
                      csl4 = csl4_down_vec))
down$X1 <- unique(c(names(wt_down_vec), names(enp1_down_vec), names(csl4_down_vec)))
down$type <- down$X1

down$type <- ifelse(grepl(pattern = "[FR]-", x = down$type),
                    "XUT", down$type)
down$type <- ifelse(grepl(pattern = "NUT", x = down$type),
                    "NUT", down$type)
down$type <- ifelse(grepl(pattern = "CUT", x = down$type),
                    "CUT", down$type)
down$type <- ifelse(grepl(pattern = "ST", x = down$type),
                    "SUT", down$type)

# color legend:
# NUT = yellow
# CUT = blue/grey
# SUT = grey
# XUT = brown
pdf(snakemake@output[["upset"]], width = 7, height = 5)
UpSetR::upset(down, 
              queries = list(
                list(query = elements, 
                     params = list("type", c("CUT","NUT", "SUT", "XUT")), color = "#446455", active = T),
                list(query = elements, 
                     params = list("type", c("NUT", "SUT", "XUT")), color = "#FDD262", active = T),
                list(query = elements, 
                     params = list("type", c("SUT", "XUT")), color = "#D3DDDC", active = T),
                list(query = elements, 
                     params = list("type", c("XUT")), color = "#C7B19C", active = T)))
dev.off()