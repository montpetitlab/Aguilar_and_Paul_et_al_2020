library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(UpSetR)

# read in differential expression results
# files <- c(snakemake@input[['wt_bm766']], snakemake@input[['csl4_bm5']],
#            snakemake@input[['csl4_bm766']], snakemake@input[['enp1_bm5']],
#            snakemake@input[['enp1_bm766']])
files <- c("outputs/deseq2/res_wt_bm766_v_wt.csv", "outputs/deseq2/res_csl4_bm5_v_wt.csv",
           "outputs/deseq2/res_csl4_bm766_v_wt.csv", "outputs/deseq2/res_enp1_bm5_v_wt.csv",
           "outputs/deseq2/res_enp1_bm766_v_wt.csv")

log2fc <- files[1:5] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("_v_wt.csv", "", source)) %>%
  replace_na(list(log2FoldChange = 0))
  

# plot --------------------------------------------------------------

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
  filter(pervasive == "pervasive")
log2fc$type <- log2fc$X1


log2fc$type <- ifelse(grepl(pattern = "[FR]-", x = log2fc$type),
                      "XUT", log2fc$type)
log2fc$type <- ifelse(grepl(pattern = "NUT", x = log2fc$type),
                      "NUT", log2fc$type)
log2fc$type <- ifelse(grepl(pattern = "CUT", x = log2fc$type),
                      "CUT", log2fc$type)
log2fc$type <- ifelse(grepl(pattern = "ST", x = log2fc$type),
                      "SUT", log2fc$type)

log2fc$change <- ifelse(log2fc$padj < 0.05, "sig", "non-sig")
log2fc$change <- ifelse(log2fc$change == "sig",
                        ifelse(log2fc$log2FoldChange > 0, "up", "down"),
                        "non-sig")
log2fc$change <- ifelse(is.na(log2fc$change), "non-sig", log2fc$change)


# REPRESSED ----------------------------------------------------------------

# log2fc_down <- log2fc %>% 
#   filter(pervasive == "pervasive") %>% 
#   filter(log2FoldChange < 0)
# 
# 
# wt_bm766_down <- log2fc_down %>%
#   filter(source == "wt_bm766") %>%
#   select(X1)
# wt_bm766_down_vec <- as.character(wt_bm766_down$X1)
# names(wt_bm766_down_vec) <- wt_bm766_down$X1
# 
# csl4_bm5_down <- log2fc_down %>%
#   filter(source == "csl4_bm5") %>%
#   select(X1)
# csl4_bm5_down_vec <- as.character(csl4_bm5_down$X1)
# names(csl4_bm5_down_vec) <- csl4_bm5_down$X1
# 
# csl4_bm766_down <- log2fc_down %>%
#   filter(source == "csl4_bm766") %>%
#   select(X1)
# csl4_bm766_down_vec <- as.character(csl4_bm766_down$X1)
# names(csl4_bm766_down_vec) <- csl4_bm766_down$X1
# 
# 
# enp1_bm5_down <- log2fc %>%
#   filter(source == "enp1_bm5") %>%
#   select(X1)
# enp1_bm5_down_vec <- as.character(enp1_bm5_down$X1)
# names(enp1_bm5_down_vec) <- enp1_bm5_down$X1
# 
# enp1_bm766_down <- log2fc %>%
#   filter(source == "enp1_bm766") %>%
#   select(X1)
# enp1_bm766_down_vec <- as.character(enp1_bm766_down$X1)
# names(enp1_bm766_down_vec) <- enp1_bm766_down$X1
# 
# 
# down <- fromList(list(wt_bm766 = wt_bm766_down_vec,
#                       enp1_bm5 = enp1_bm5_down_vec,
#                       enp1_bm766 = enp1_bm766_down_vec,
#                       csl4_bm5 = csl4_bm5_down_vec,
#                       csl4_bm766 = csl4_bm766_down_vec))
# down$X1 <- unique(c(names(wt_bm766_down_vec),
#                     names(enp1_bm5_down_vec),
#                     names(enp1_bm766_down_vec),
#                     names(csl4_bm5_down_vec),
#                     names(csl4_bm766_down_vec)))
# down$type <- down$X1
# 
# 
# down$type <- ifelse(grepl(pattern = "[FR]-", x = down$type),
#                     "XUT", down$type)
# down$type <- ifelse(grepl(pattern = "NUT", x = down$type),
#                     "NUT", down$type)
# down$type <- ifelse(grepl(pattern = "CUT", x = down$type),
#                     "CUT", down$type)
# down$type <- ifelse(grepl(pattern = "ST", x = down$type),
#                     "SUT", down$type)
# 
# # color legend:
# # NUT = yellow
# # CUT = blue/grey
# # SUT = grey
# # XUT = brown
# # pdf(snakemake@output[["upset"]], width = 7, height = 5)
# UpSetR::upset(down, 
#               queries = list(
#                 list(query = elements, 
#                      params = list("type", c("CUT","NUT", "SUT", "XUT")), color = "#446455", active = T),
#                 list(query = elements, 
#                      params = list("type", c("NUT", "SUT", "XUT")), color = "#FDD262", active = T),
#                 list(query = elements, 
#                      params = list("type", c("SUT", "XUT")), color = "#D3DDDC", active = T),
#                 list(query = elements, 
#                      params = list("type", c("XUT")), color = "#C7B19C", active = T)))
# # dev.off()
# 

# correlation matrix ------------------------------------------------------

# log2fc_wide <- log2fc %>%
#   dplyr::select(source, X1, log2FoldChange) %>%
#   filter(grepl(pattern = "CUT", X1)) %>%
#   pivot_wider(id_cols = c(X1),
#               names_from = source, values_from = log2FoldChange) %>%
#   replace_na(list(wt_bm766 = 0, csl4_bm5 = 0, csl4_bm766 = 0, enp1_bm5 = 0, enp1_bm766 = 0))
# 
# 
# correlation <- cor(log2fc_wide %>% 
#                      dplyr::select(wt_bm766, csl4_bm5, csl4_bm766, enp1_bm5, enp1_bm766) %>%
#                      filter(!is.na(wt_bm766)))
# # pdf(snakemake@output[['rna_plt_1']], width = 5, height = 4)
# corrplot(correlation, method = "color", type="upper", addCoef.col = "white")
# # dev.off()
# rcorr_pvalues <- rcorr.adjust(as.matrix(res_wide[, c(4:8)]))
# print(rcorr_pvalues)
