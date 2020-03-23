library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(UpSetR)

# read in differential expression results
files <- list.files("outputs/deseq2", full.names = T)
files <- files[grepl(pattern = "sig\\.csv$", x = files)]

log2fc <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("_sig.csv", "", source))

# reset levels so control comes first
# log2fc$source <- factor(log2fc$source, 
#                         levels = c("wt", "csl4", "enp1"))


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
pdf("outputs/figures/uts_upset.pdf", width = 7, height = 5)
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
