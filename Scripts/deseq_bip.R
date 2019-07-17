options( warn = -1 )
library('DESeq2')
library('data.table')
library('stats')
library('ggplot2')
library('tidyverse')
setwd("~/Desktop/Paul_et_al_2019/")
wt_vs_mutants_d <- fread("./results/output_feature_count.csv", header = T, sep = ",")
#---------------------------------- DESEQ2 RB ---------------------------------------------------

#Choosing columns for RB data
sample_name = c('GeneID', 'Rrp6_1', 'Rrp6_2', 'Dis3_1', 'Dis3_2', 'Csl4_90_1_rb', 'Csl4_90_2_rb', 'Csl4_90_3_rb',
                'Enp1_90_1_rb', 'Enp1_90_2_rb','Srm1_90_1_rb', 'Srm1_90_2_rb', 'wt_90_1', 'wt_90_2' )
wt_vs_mutants_d <- wt_vs_mutants_d[!duplicated(wt_vs_mutants_d$GeneID), ]
colnames(wt_vs_mutants_d) <- sample_name

sample_name = c('Rrp6_1', 'Rrp6_2', 'Dis3_1', 'Dis3_2', 'Csl4_90_1_rb', 'Csl4_90_2_rb', 'Csl4_90_3_rb',
                'Enp1_90_1_rb', 'Enp1_90_2_rb','Srm1_90_1_rb', 'Srm1_90_2_rb', 'wt_90_1', 'wt_90_2' )


#Selecting data
expr_rb = wt_vs_mutants_d[, sample_name, with=FALSE]


#Setting up colData
col_rb = data.frame(sample = sample_name,
                    treatment = c('rrp6_90','rrp6_90', 'dis3_90', 'dis3_90', 'csl4_90','csl4_90','csl4_90', 'enp1_90', 'enp1_90',
                                  'srm1_90', 'srm1_90','Ctrl', 'Ctrl'))

col_rb$treatment <- factor(col_rb$treatment, levels=c('rrp6_90', 'dis3_90', 'csl4_90', 'enp1_90', 'srm1_90', 'Ctrl'))
                                                       
                                                    

dds_rb <- DESeqDataSetFromMatrix(countData = data.frame(expr_rb),
                                 colData = col_rb,
                                 design = ~ treatment)
rownames(dds_rb) <- wt_vs_mutants_d$GeneID

#Applying the model
dds_rb <- estimateSizeFactors(dds_rb)
dds_rb <- estimateDispersions(dds_rb)
#res_lrt <- DESeq(dds_rb, test='LRT', reduced = ~1)
res_wald <- DESeq(dds_rb, test='Wald')

#getting results ...
res_rrp6_90 <- results(res_wald, pAdjustMethod = "BH", contrast=c('treatment', 'rrp6_90', 'Ctrl'))
res_dis3_90 <- results(res_wald, pAdjustMethod = "BH", contrast=c('treatment', 'dis3_90', 'Ctrl'))
res_csl4_90 <- results(res_wald, pAdjustMethod = "BH", contrast=c('treatment', 'csl4_90', 'Ctrl'))
res_enp1_90 <- results(res_wald, pAdjustMethod = "BH", contrast=c('treatment', 'enp1_90', 'Ctrl'))
res_srm1_90 <- results(res_wald, pAdjustMethod = "BH", contrast=c('treatment', 'srm1_90', 'Ctrl'))



#saving results for RB data
tmp <- cbind(as.data.frame(res_rrp6_90), as.data.frame(res_dis3_90), as.data.frame(res_csl4_90), as.data.frame(res_enp1_90), as.data.frame(res_srm1_90))
log2fc <- tmp[ , grepl( "log2FoldChange" , names( tmp) ) ]
p_v <- tmp[ , grepl( "pvalue" , names( tmp ) ) ]
p_v <- p_v[apply(p_v[, -1], MARGIN = 1, function(x) any(x < 0.01)), ]
new_log2fc <- log2fc[rownames(p_v),]
colnames(new_log2fc) = c('rrp6', 'dis3', 'csl4', 'enp1', 'srm1')
new_log2fc<- new_log2fc %>% rownames_to_column("gene_id")
write.table(new_log2fc, file="./results/log2FC_all_rb.tsv", quote = F, sep='\t', row.names = TRUE, col.names = TRUE)

