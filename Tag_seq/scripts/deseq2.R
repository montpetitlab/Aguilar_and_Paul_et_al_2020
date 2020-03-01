library(DESeq2)

## Read in metadata and counts
info <- read.csv(snakemake@input[['samples']])
counts <- read.csv(snakemake@input[['counts']], row.names = 1)
counts <- as.data.frame(apply(counts, 1:2, round))

## Run differential expression 
dds <- DESeqDataSetFromMatrix(counts,
                              colData = info,
                              design = ~ group)
ds <- DESeq(dds, test="Wald")
# resultsNames(ds)

## Extract and write results
# extract significant genes for each PAB1 vs. empty.
# listing PAB1 (BM766) first means that genes with a small p value and 
# positive log FC are more highly expressed in PAB1 than in empty.
res_wt <- results(ds, contrast = c("group", "BM_WT_pBM766", "BM_WT_pBM5"))
write.csv(res_wt, snakemake@output[['wt']], quote = F)
res_wt_sig <- subset(res_wt, padj < 0.01)
write.csv(res_wt_sig, snakemake@output[['wt_sig']], quote = F)

res_csl4 <- results(ds, contrast = c("group", "BM_csl4ph_pBM766", "BM_csl4ph_pBM5"))
write.csv(res_csl4, snakemake@output[['csl4']], quote = F)
res_csl4_sig <- subset(res_csl4, padj < 0.01)
write.csv(res_csl4_sig, snakemake@output[['csl4_sig']], quote = F)

res_enp1 <- results(ds, contrast = c("group", "BM_enp1_1_pBM766", "BM_enp1_1_pBM5"))
write.csv(res_enp1, snakemake@output[['enp1']], quote = F)
res_enp1_sig <- subset(res_enp1, padj < 0.01)
write.csv(res_enp1_sig, snakemake@output[['enp1_sig']], quote = F)
