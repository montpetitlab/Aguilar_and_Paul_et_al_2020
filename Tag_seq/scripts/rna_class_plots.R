# RNA class
library(ggplot2)
library(readr)
library(org.Sc.sgd.db)
library(tidyr)
library(purrr)
library(dplyr)

# read in RNA classes and convert to long format
rna_class <- read_csv(snakemake@input[['rna']], skip = 1) %>%
  dplyr::select(I, II, III, IV, V, VI, VII, VIII, IX, X) %>%
  pivot_longer(cols = I:X, names_to = "class", values_to = "gene")

# convert gene names to ORF names
# all pervasive trancripts and genes that are already in ORF format 
# will fail and print and error message
orf <- sapply(rna_class$gene, function(x) try(get(x, org.Sc.sgdCOMMON2ORF)))

# parse the output to a dataframe that maps the two values
lookup <- data.frame(common = rep(NA, length(orf)), orf = rep(NA, length(orf)))
#lookup_common <- vector()
#lookup_orf <- vector()
x <- 0
for(i in 1:length(orf)){
  len <- length(orf[[i]])
  for(j in 1:len){
    x <- length(lookup$common[!is.na(lookup$common)]) + 1 # create a new index for placement into df
    lookup$common[x] <- names(orf[i])
    if(class(orf[[i]]) == "try-error") {
      next
    } else {
      lookup$orf[x] <- orf[[i]][j]
    }
  }
}

# remove CUTs and SUTs from the dataframe
lookup <- lookup %>%
  filter(!grepl(pattern = "CUT", common)) %>%
  filter(!grepl(pattern = "SUT", common)) %>%
  filter(!is.na(common))

# place transcripts that were already annotated as ORFs in ORF col
lookup$orf <- ifelse(is.na(lookup$orf), lookup$common, lookup$orf)

# translate RNA classes from common to ORF
rna_class <- left_join(rna_class, lookup, by = c("gene" = "common"))
# remove NAs (aka SUTs and CUTS)
rna_class <- rna_class %>%
  filter(!is.na(orf))

# read in results of differential expression

files <- c(snakemake@input[['csl4']], snakemake@input[['enp1']],
           snakemake@input[['wt']])

res <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("\\.csv", "", source)) %>%
  filter(padj < .05) %>%
  dplyr::select(source, X1, log2FoldChange)

# separate transcript names from coords
res <- separate(data = res, col = X1, into = c("transcript", "coord"), sep = "::")
res$transcript <- gsub("_mRNA", "", res$transcript)

# filter results to rna classes
res <- res %>%
  filter(transcript %in% rna_class$orf)

# join with class info
res <- left_join(res, rna_class, by = c("transcript" = "orf"))

pdf(snakemake@output[["rna_plt_all"]], height = 3.5, width = 7)
ggplot(res, aes(x = class, y = log2FoldChange)) +
  theme_minimal() +
  geom_boxplot(outlier.size = -1) + 
  facet_wrap(~source) +
  ylim(c(-5, 7))
dev.off()


pdf(snakemake@output[["rna_plt_1_10"]], height = 3.5, width = 7)
ggplot(res %>%
         filter(class %in% c("I", "X")), 
       aes(x = class, y = log2FoldChange)) +
  theme_minimal() +
  geom_boxplot(outlier.size = -1) + 
  facet_wrap(~source) +
  ylim(c(-4, 2.5))
dev.off()