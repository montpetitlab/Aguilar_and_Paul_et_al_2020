library(DESeq2)
library(tximport)
library(readr)

samples <- read_csv(snakemake@input[['samples']])
counts <- tximport(files = samples$filename, type = "salmon", txOut = TRUE)
## Read in metadata and counts

## Run differential expression 
dds <- DESeqDataSetFromTximport(counts,
                                colData = samples,
                                design = ~ group)
ds <- DESeq(dds, test="Wald")
# resultsNames(ds)

## Extract and write results
# designating mutant first means that genes with a small p value and 
# positive log FC are more highly expressed in mutant than in wt.

res_csl4 <- results(ds, contrast = c("group", "csl4ph", "ctrl"))
write.csv(res_csl4, snakemake@output[['csl4']], quote = F)

res_enp1 <- results(ds, contrast = c("group", "enp1", "ctrl"))
write.csv(res_enp1, snakemake@output[['enp1']], quote = F)

res_srm1 <- results(ds, contrast = c("group", "srm1", "ctrl"))
write.csv(res_srm1, snakemake@output[['srm1']], quote = F)
