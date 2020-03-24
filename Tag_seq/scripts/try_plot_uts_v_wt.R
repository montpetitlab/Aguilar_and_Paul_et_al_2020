library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# read in differential expression results
files <- c(snakemake@input[['wt_bm766']], snakemake@input[['csl4_bm5']],
           snakemake@input[['csl4_bm766']], snakemake@input[['enp1_bm5']],
           snakemake@input[['enp1_bm766']])
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


# log2fc_tally <- log2fc %>%
#   group_by(source, type, change) %>%
#   tally()
# 
# 
# ggplot(log2fc_tally, aes(x = source, y = n, fill = type)) +
#   geom_col() +
#   theme_minimal() +
#   facet_wrap(~change, scales = "free")



# try scatter plots -------------------------------------------------------
# log2fc_wide <- log2fc %>%
#   dplyr::select(source, X1, log2FoldChange, type) %>%
#   #filter(grepl(pattern = "CUT", X1)) %>%
#   pivot_wider(id_cols = c(X1, type),
#               names_from = source, values_from = log2FoldChange) %>%
#   replace_na(list(wt_bm766 = 0, csl4_bm5 = 0, csl4_bm766 = 0, enp1_bm5 = 0, enp1_bm766 = 0))


# ggplot(log2fc_wide) +
#   geom_point(aes(x = wt_bm766, y = enp1_bm5), color = 'black', alpha = 1/10) +
#   geom_point(aes(x = wt_bm766, y = enp1_bm766), color = 'red', alpha = 1/10) +
#   theme_minimal()
# 
# ggplot(log2fc_wide) +
#   geom_point(aes(x = enp1_bm5, y = enp1_bm766), color = 'black', alpha = 1/10) +
#   theme_minimal() +
#   lims(x = c(-10, 10), y = c(-10, 10))
# 
# 
# ggplot(log2fc %>% filter(source == "enp1_bm5")) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), color = change), alpha = 1/10) +
#   theme_minimal() +
#   lims(x = c(-10, 10), y = c(0, 20))
# 
# ggplot(log2fc %>% filter(source == "enp1_bm766")) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), color = change), alpha = 1/10) +
#   theme_minimal() +
#   lims(x = c(-10, 10), y = c(0, 20))
# 
# 
# ggplot(log2fc %>% 
#          filter(change != "non-sig") %>%
#          filter(source %in% c("enp1_bm5", "enp1_bm766", "wt_bm766")), 
#        aes(x = log2FoldChange, fill = source)) +
#   geom_histogram(alpha = .5) +
#   theme_minimal() +
#   xlim(c(-10, 10))

log2fc$plasmid <- ifelse(grepl("bm5", log2fc$source), "bm5", "bm766")
log2fc$strain <- ifelse(grepl("wt", log2fc$source), "wt", NA)
log2fc$strain <- ifelse(grepl("enp1", log2fc$source), "enp1", log2fc$strain)
log2fc$strain <- ifelse(grepl("csl4", log2fc$source), "csl4", log2fc$strain)

ggplot(log2fc %>% 
         filter(source %in% c("enp1_bm5", "enp1_bm766", "wt_bm766")), 
       aes(x = log2FoldChange, color = strain)) +
  geom_freqpoly(aes(linetype = plasmid)) +
  theme_minimal() +
  #facet_wrap(~source) + 
  xlim(c(-10, 10)) + 
  scale_color_manual(values = c(wt = "black", enp1 = "#129060"))


ggplot(log2fc %>% 
         filter(source %in% c("csl4_bm5", "csl4_bm766", "wt_bm766")), 
       aes(x = log2FoldChange, color = strain)) +
  geom_freqpoly(aes(linetype = plasmid)) +
  theme_minimal() +
  #facet_wrap(~source) + 
  xlim(c(-10, 10)) + 
  scale_color_manual(values = c(wt = "black", csl4 = "#47a4e4"))


ggplot(log2fc, aes(x = log2FoldChange, color = strain)) +
  geom_freqpoly(aes(linetype = plasmid)) +
  theme_minimal() +
  facet_wrap(~type, scales = "free") + 
  xlim(c(-10, 10)) + 
  scale_color_manual(values = c(wt = "black", csl4 = "#47a4e4", enp1 = "#129060"))
