# RNA class
library(ggplot2)
library(readr)
library(dplyr)
library(org.Sc.sgd.db)
library(tidyr)
library(purrr)
library(corrplot)      
library(Hmisc)          
library(RcmdrMisc)

# read in RNA classes and convert to long format
rna_class <- read_csv(snakemake@input[['rna']], skip = 1) %>%
  dplyr::select(I, II, III, IV, V, VI, VII, VIII, IX, X) %>%
  pivot_longer(cols = I:X, names_to = "class", values_to = "gene")

rna_class <- read_csv("inputs/rna_classes.csv", skip = 1) %>%
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

files <- c(snakemake@input[['wt_bm766']], snakemake@input[['csl4_bm5']],
           snakemake@input[['csl4_bm766']], snakemake@input[['enp1_bm5']],
           snakemake@input[['enp1_bm766']])
files <- c("outputs/deseq2/res_wt_bm766_v_wt.csv", "outputs/deseq2/res_csl4_bm5_v_wt.csv",
           "outputs/deseq2/res_csl4_bm766_v_wt.csv", "outputs/deseq2/res_enp1_bm5_v_wt.csv",
           "outputs/deseq2/res_enp1_bm766_v_wt.csv")
res <- files %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs\\/deseq2\\/res_", "", source)) %>%
  mutate(source = gsub("_v_wt.csv", "", source)) %>%
  #filter(padj < .05) %>%
  dplyr::select(source, X1, log2FoldChange, padj) %>%
  replace_na(list(log2FoldChange = 0))


# separate transcript names from coords
res <- separate(data = res, col = X1, into = c("transcript", "coord"), sep = "::")
res$transcript <- gsub("_mRNA", "", res$transcript)

# filter results to rna classes
res <- res %>%
  filter(transcript %in% rna_class$orf)

# join with class info
res <- left_join(res, rna_class, by = c("transcript" = "orf"))

res$sig <- ifelse(res$padj < .05, "sig", "nonsig")


# plot all ----------------------------------------------------------------

pdf(snakemake@output[["rna_plt_all"]], height = 7, width = 7)
ggplot(res, aes(x = class, y = log2FoldChange)) +
  theme_minimal() +
  geom_boxplot(outlier.size = -1) +
  facet_wrap(~source) +
  ylim(c(-7, 7))
dev.off()

# plot accumulation: I -------------------------------------------------------

resI <- res %>% 
  filter(class == "I") %>%
  group_by(source, class) %>% 
  arrange(log2FoldChange) %>% 
  mutate(rn = row_number()) %>%
  mutate(prop = rn/length(unique(rn)))


resI$plasmid <- ifelse(grepl("bm5", resI$source), "bm5", "bm766")
resI$strain <- ifelse(grepl("wt", resI$source), "wt", NA)
resI$strain <- ifelse(grepl("enp1", resI$source), "enp1", resI$strain)
resI$strain <- ifelse(grepl("csl4", resI$source), "csl4", resI$strain)


pdf(snakemake@output[["rna_I_accum"]], height = 6, width = 7)
ggplot(resI) + 
  geom_line(aes(x=log2FoldChange, y= prop, color=strain, linetype = plasmid), size = 1.1, alpha = .7) +
  theme_minimal() + 
  facet_wrap(~class, scales = "free") + 
  theme(legend.position = "bottom", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Total Fraction of Transcripts in Class",
       x = "log2FC (mutant vs. control with bm5)") +
  scale_color_manual(values = c(csl4 = "#47a4e4", wt = "black", 
                                enp1 = "#129060"))

dev.off()


# plot accumulation: I-III -------------------------------------------------------

resIII <- res %>% 
  filter(class %in% c("I", "II", "III")) %>%
  group_by(source) %>% 
  arrange(log2FoldChange) %>% 
  mutate(rn = row_number()) %>%
  mutate(prop = rn/length(unique(rn)))


resIII$plasmid <- ifelse(grepl("bm5", resIII$source), "bm5", "bm766")
resIII$strain <- ifelse(grepl("wt", resIII$source), "wt", NA)
resIII$strain <- ifelse(grepl("enp1", resIII$source), "enp1", resIII$strain)
resIII$strain <- ifelse(grepl("csl4", resIII$source), "csl4", resIII$strain)


pdf(snakemake@output[["rna_III_accum"]], height = 6, width = 7)
ggplot(resIII) + 
  geom_line(aes(x=log2FoldChange, y= prop, color=strain, linetype = plasmid), size = 1.1, alpha = .7) +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Total Fraction of Transcripts in Class",
       x = "log2FC (mutant vs. control with bm5)") +
  scale_color_manual(values = c(csl4 = "#47a4e4", wt = "black", 
                                enp1 = "#129060"))

dev.off()

# plot accumulation: X -------------------------------------------------------

resX <- res %>% 
  filter(class == "X") %>%
  group_by(source, class) %>% 
  arrange(log2FoldChange) %>% 
  mutate(rn = row_number()) %>%
  mutate(prop = rn/length(unique(rn)))


resX$plasmid <- ifelse(grepl("bm5", resX$source), "bm5", "bm766")
resX$strain <- ifelse(grepl("wt", resX$source), "wt", NA)
resX$strain <- ifelse(grepl("enp1", resX$source), "enp1", resX$strain)
resX$strain <- ifelse(grepl("csl4", resX$source), "csl4", resX$strain)


pdf(snakemake@output[["rna_X_accum"]], height = 6, width = 7)
ggplot(resX) + 
  geom_line(aes(x=log2FoldChange, y= prop, color=strain, linetype = plasmid), size = 1.1, alpha = .7) +
  theme_minimal() + 
  facet_wrap(~class, scales = "free") + 
  theme(legend.position = "bottom", 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(y = "Total Fraction of Transcripts in Class",
       x = "log2FC (mutant vs. control with bm5)") +
  scale_color_manual(values = c(csl4 = "#47a4e4", wt = "black", 
                                enp1 = "#129060"))
dev.off()


# correlation coefficient -------------------------------------------------

res_wide <- res %>%
  filter(class == "I") %>%
  dplyr::select(-sig, -padj, -class) %>%
  pivot_wider(id_cols = c(gene, transcript, coord),
              names_from = source, values_from = log2FoldChange)


correlation <- cor(res_wide %>% 
                     dplyr::select(wt_bm766, csl4_bm5, csl4_bm766, enp1_bm5, enp1_bm766) %>%
                     filter(!is.na(wt_bm766)))
pdf(snakemake@output[['corr']], width = 5, height = 4)
corrplot(correlation, method = "color", type="upper", addCoef.col = "white")
dev.off()
rcorr_pvalues <- rcorr.adjust(as.matrix(res_wide[, c(4:8)]))
print(rcorr_pvalues)