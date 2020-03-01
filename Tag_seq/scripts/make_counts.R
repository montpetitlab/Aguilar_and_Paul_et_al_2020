library(tximport)
library(readr)

tx2gene <- read_tsv(snakemake@input[['tx2gene']])
samples <- read_csv(snakemake@input[['samples']])
counts <- tximport(files = samples$filename, type = "salmon", tx2gene = tx2gene)
# 3' TagSeq transcripts all originate from the 3' end of the transcript,
# and so the transcript length information imported with salmon is 
# irrelevant for count normalization. 
counts <- counts$counts
colnames(counts) <- samples$sample
write.csv(counts, snakemake@output[['counts']], 
          quote = F, row.names = T)
