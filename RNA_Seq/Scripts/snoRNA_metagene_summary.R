library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

# generate column names
mat_colnames <- read_tsv(snakemake@input[["mat"]], skip = 2, n_max = 1, col_names = F)
mat_colnames <- as.character(as.data.frame(mat_colnames[1, ]))
mat_colnames <- c(mat_colnames[-1], "tmp")
# get number of bins
bins <- unname(table(mat_colnames)[1])
mat_colnames <- paste0(mat_colnames, "+", 1:bins)

mat <- read_tsv(snakemake@input[["mat"]], skip = 3, col_names = mat_colnames)
# remove last column
mat <- mat[ , -ncol(mat)]
# transform data to long format
mat_long <- pivot_longer(data = mat, cols = `HI.3202.005.Index_20.BM_shift90-NAB2-PrA-CSL4pH-URA-37C-d1+1`:`HI.4535.001.Index_4.CA_WT-2+140`, 
                         names_to = "sample_bin", values_to = "rpkm")

mat_long <- separate(data = mat_long, col = sample_bin, sep = "\\+",
                     into = c("sample", "bin"))

# make a grouping column to combine csl4-ph, enp1, and wt samples together
mat_long$strain <- NA
mat_long$strain <- ifelse(grepl(pattern = "CSL4pH", mat_long$sample), "csl4ph", mat_long$strain)
mat_long$strain <- ifelse(grepl(pattern = "enp1", mat_long$sample), "enp1", mat_long$strain)
mat_long$strain <- ifelse(grepl(pattern = "WT", mat_long$sample), "wt", mat_long$strain)

mat_grouped <- mat_long %>%
  filter(!is.na(strain)) %>%
  group_by(bin, strain) %>%
  summarize(mean_rpkm = mean(rpkm))

mat_grouped$bin <- as.numeric(mat_grouped$bin)

pdf(snakemake@output[['plt']], height = 4, width = 7)
ggplot(mat_grouped, aes(x = bin, y = mean_rpkm, color = strain)) +
  geom_line(size=1.2) + 
  labs(x ="snoRNA genes", y = "RPKM") + 
  geom_vline(aes(xintercept = 16), linetype = "dashed") +
  geom_vline(aes(xintercept = max(bin) - 16), linetype = "dashed") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_manual(values = c(wt = "#000000", csl4ph = "#47a4e4", enp1 = "#129060"))
dev.off()

