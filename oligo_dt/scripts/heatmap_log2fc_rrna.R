library(tximport)
library(dplyr)
library(readr)
library(DESeq2)

setwd("~/github/Paul_et_al_2019/oligo_dt/")
rRNA <- c("ETS1-1", "RDN18-1", "ITS1-1", "RDN58-1", "ITS2-1", 
          "RDN25-1", "ETS2-1", "RDN5-1", "RDN5-2") # "RDN37-1"
# rRNA <- c('ETS2-1::chrXII:451575-451785', 'RDN37-1::chrXII:451575-458432', 
#           'RDN25-1::chrXII:451786-455181', 'ITS2-1::chrXII:455182-455413', 
#           'RDN58-1::chrXII:455414-455571', 'ITS1-1::chrXII:455572-455932', 
#           'RDN18-1::chrXII:455933-457732', 'ETS1-1::chrXII:457733-458432', 
#           'RDN5-1::chrXII:459676-459796', 'RDN5-2::chrXII:468813-468931')



files <- unlist(snakemake@input)
files <- list.files("outputs/deseq2", full.names = T)
res <- files[1:3] %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  select(source, X1, log2FoldChange) %>%
  separate(col = X1, into = c("transcript"), sep = "::", remove = T) %>%
  filter(transcript %in% rRNA)

# pivot wide so each experiment has a col of log2fold changes
res_wide <- pivot_wider(data = res, id_cols = transcript, names_from = source, values_from = log2FoldChange)

res_wide <- res_wide %>%
  tidyr::replace_na(list(csl4 = 0, enp1 = 0, srm1 = 0))

# make gene names rownames
res_wide <- as.data.frame(res_wide)
rownames(res_wide) <- res_wide$transcript
#res_dropped <- rbind(c("RDN18-1", 0, 0, 0), c("ETS2-1", 0, 0, 0))
#rownames(res_dropped) <- c("RDN18-1", "ETS2-1")
#colnames(res_dropped) <- colnames(res_wide)
#res_wide <- rbind(res_wide, res_dropped)
res_wide <- res_wide[order(match(res_wide$transcript, rRNA)), ]
res_wide <- res_wide[ , -1]
res_wide <- as.matrix(sapply(res_wide, as.numeric)) 
rownames(res_wide) <- rRNA


pdf(snakemake@output[["heatmap"]], width = 6, height = 6)
pheatmap(mat                  = res_wide,
         show_rownames        = T,
         show_colnames        = T,
         cluster_rows         = F,
         cluster_cols         = T,
         annotation_names_row = F, 
         display_numbers      = T,
         number_color         = "black",
         color                 = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(21),
         breaks               = seq(-20, 20, by = 2))
dev.off()

