# percent of mapped reads that fall into each category
library(readr)
library(tximport)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(dplyr)

## Read in metadata and counts
samples <- read_csv(snakemake@input[['samples']])
counts <- tximport(files = samples$filename, type = "salmon", txOut = TRUE)

## extract normalized counts
ddsTxi <- DESeqDataSetFromTximport(counts,
                                   colData = samples,
                                   design = ~ group)

dds <- DESeq(ddsTxi)
dds <- estimateSizeFactors(dds)
dds_counts <- counts(dds, normalized=TRUE)
colnames(dds_counts) <- colData(dds)$sample
dds_counts <- as.data.frame(dds_counts)
dds_counts$transcript <- rownames(dds_counts)

## define transcript class based on transcript name
dds_counts <- separate(data = dds_counts, col = transcript, 
                       into = c("class", "coord"), sep = "::", remove = F)
## edit class
## CUT == CUT
## NUT == NUT
## ST  == SUT 

## *_mRNA == mRNA
## ..-... == XUT

dds_counts$class <- ifelse(grepl(pattern = "mRNA", x = dds_counts$class),
                                "mRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "NME1", x = dds_counts$class),
                           "snoRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "ST", x = dds_counts$class),
                           "pervasive", dds_counts$class)
snRNA <- c("LSR1", "snR14", "snR7-L", "snR7-S", "snR6", "snR19")
dds_counts$class <- ifelse(grepl(pattern = "sn", x = dds_counts$class),
                           ifelse(!dds_counts$class %in% snRNA,
                           "snoRNA", "snRNA"), dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "LSR1", x = dds_counts$class),
                           "snRNA", dds_counts$class)
rRNA <- c("ETS2-1", "RDN37-1", "RDN25-1", "ITS2-1", "RDN58-1", "ITS1-1", 
          "RDN18-1", "ETS1-1", "RDN5-1", "RDN5-2")
dds_counts$class <- ifelse(dds_counts$class %in% rRNA, "rRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "[FR]-", x = dds_counts$class),
                           "pervasive", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "NUT", x = dds_counts$class),
                           "pervasive", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "CUT", x = dds_counts$class),
                           "pervasive", dds_counts$class)

# transform to long 
dds_counts_long <- dds_counts %>%
  dplyr::select(-transcript, -coord) %>%
  pivot_longer(names_to = "sample", values_to = "norm_count", -class)

# collapse to group
dds_counts_long <- left_join(dds_counts_long, samples)

# count the number of reads quantified per class
dds_counts_long_summarized <- dds_counts_long %>%
  group_by(group, class) %>%
  summarise(sum_class = sum(norm_count))

# count the total number of reads quantified
dds_counts_sample <- dds_counts_long %>%
  group_by(group) %>%
  summarize(sum_sample = sum(norm_count))

# join counts together
dds_counts_summarized <- left_join(dds_counts_long_summarized, dds_counts_sample)

# calculate percent mapping
dds_counts_summarized$percent_mapping <- dds_counts_summarized$sum_class / dds_counts_summarized$sum_sample


pdf(snakemake@output[["perc_map_ribo"]], height = 5, width = 7)
ggplot(dds_counts_summarized, aes(x = group, y = percent_mapping, fill = class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "fraction of reads mapped (ribo -)", x = "") +
  coord_flip() +
  scale_fill_manual(values = c(mRNA = "#999999", pervasive = "#377EB8", 
                               snoRNA = "#FF8C00", snRNA = "#000000",
                               rRNA = "#4B0082")) 
dev.off()

pdf(snakemake@output[["perc_map_ribo_zoom"]], height = 5, width = 7)
ggplot(dds_counts_summarized, aes(x = group, y = percent_mapping, fill = class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "fraction of reads mapped (ribo -)", x = "") +
  coord_flip(ylim = c(0, .26)) +
  scale_fill_manual(values = c(mRNA = "#999999", pervasive = "#377EB8", 
                               snoRNA = "#FF8C00", snRNA = "#000000",
                               rRNA = "#4B0082")) +
  theme(legend.position="none")
dev.off()

  
