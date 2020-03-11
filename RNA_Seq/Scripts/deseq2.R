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

res_csl4ph <- results(ds, contrast = c("group", "csl4ph", "wt"))
write.csv(res_csl4ph, snakemake@output[['csl4ph']], quote = F)
res_csl4ph_sig <- subset(res_csl4ph, padj < 0.01)
write.csv(res_csl4ph_sig, snakemake@output[['csl4ph_sig']], quote = F)

res_dis3 <- results(ds, contrast = c("group", "dis3", "wt"))
write.csv(res_dis3, snakemake@output[['dis3']], quote = F)
res_dis3_sig <- subset(res_dis3, padj < 0.01)
write.csv(res_dis3_sig, snakemake@output[['dis3_sig']], quote = F)

res_enp1 <- results(ds, contrast = c("group", "enp1", "wt"))
write.csv(res_enp1, snakemake@output[['enp1']], quote = F)
res_enp1_sig <- subset(res_enp1, padj < 0.01)
write.csv(res_enp1_sig, snakemake@output[['enp1_sig']], quote = F)

res_rrp6 <- results(ds, contrast = c("group", "rrp6", "wt"))
write.csv(res_rrp6, snakemake@output[['rrp6']], quote = F)
res_rrp6_sig <- subset(res_rrp6, padj < 0.01)
write.csv(res_rrp6_sig, snakemake@output[['rrp6_sig']], quote = F)

res_srm1 <- results(ds, contrast = c("group", "srm1", "wt"))
write.csv(res_srm1, snakemake@output[['srm1']], quote = F)
res_srm1_sig <- subset(res_srm1, padj < 0.01)
write.csv(res_srm1_sig, snakemake@output[['srm1_sig']], quote = F)
