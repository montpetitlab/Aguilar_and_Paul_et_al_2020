# percent of mapped reads that fall into each category
setwd("~/github/Paul_et_al_2019/RNA_Seq/")
library(readr)
library(tximport)
library(DESeq2)
library(tidyr)
library(ggplot2)

## Read in metadata and counts
# samples <- read_csv(snakemake@input[['samples']])
samples <- read_csv("inputs/ribo-samples.csv")
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
## sn* == snNRA
## *_mRNA == mRNA
## ..-... == XUT

dds_counts$class <- ifelse(grepl(pattern = "mRNA", x = dds_counts$class),
                                "mRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "NME1", x = dds_counts$class),
                           "mRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "ST", x = dds_counts$class),
                           "pervasive", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "sn", x = dds_counts$class),
                           "snRNA", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "[FR]-", x = dds_counts$class),
                           "pervasive", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "NUT", x = dds_counts$class),
                           "pervasive", dds_counts$class)
dds_counts$class <- ifelse(grepl(pattern = "CUT", x = dds_counts$class),
                           "pervasive", dds_counts$class)

# transform to long 
dds_counts_long <- dds_counts %>%
  select(-transcript, -coord) %>%
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


ggplot(dds_counts_summarized, aes(x = group, y = percent_mapping, fill = class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(y = "fraction of reads mapped (ribo -)", x = "") +
  coord_flip() +
  scale_fill_manual(values = c(mRNA = "#999999", pervasive = "#377EB8", snRNA = "#000000")) 

  
