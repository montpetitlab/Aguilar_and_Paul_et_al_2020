library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(corrplot)

# read in differential expression results into single dataframe
files <- list.files("RNA_Seq/outputs/deseq2", ".csv$", full.names = T)
files <- files[!grepl(pattern = "*sig*", x = files)]
res <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("RNA_Seq\\/outputs\\/deseq2\\/", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  select(source, X1, log2FoldChange)

# pivot wide so each experiment has a col of log2fold changes
res_wide <- pivot_wider(data = res, id_cols = X1, names_from = source, values_from = log2FoldChange)
# remove rows with NAs
res_wide <- filter(res_wide, !is.na(res_srm1))
# make gene names rownames
res_wide <- as.data.frame(res_wide)
rownames(res_wide) <- res_wide$X1
res_wide <- res_wide[ , -1]

# correlate
correlation <- cor(res_wide)
pdf(snakemake@output[['corr_plot']], width = 5, height = 4)
corrplot(correlation, method = "color", type="upper", addCoef.col = "white")
dev.off()
