options( warn = -1 )
library('edgeR')
library('data.table')
library('stats')
library('ggplot2')
library(dplyr)
library(tidyverse)
setwd("~/Paul_et_al_2019/Tag_seq/")
dr <- list.files("salmon_out", pattern = "*count.sf*", recursive = T, full.names = T)
salmon_out <- read.csv(dr[1], sep = "\t")
sample_name <- sub("_001.trim.trimbb", "", strsplit(dr[1], split = "/")[[1]][2])
colnames(salmon_out) <- c("Name", sample_name)
for (i in 2:length(dr)){
  temp <- read.csv(dr[i], sep = "\t")
  sample_name <- sub("_001.trim.trimbb", "", strsplit(dr[i], split = "/")[[1]][2])
  colnames(temp) <- c("Name", sample_name)
  salmon_out <- merge(salmon_out, temp, by = "Name", all = T)
  }
write.csv(salmon_out, file = "salmon_output.csv")
counts <- readDGE(dr)
group <- c("csl4ph_pBM5", "csl4ph_pBM5","csl4ph_pBM5","csl4ph_pBM766", "csl4ph_pBM766",
           "csl4ph_pBM766", "enp1_1_pBM5", "enp1_1_pBM5", "enp1_1_pBM5",  
           "enp1_1_pBM766", "enp1_1_pBM766", "enp1_1_pBM766", "WT_pBM5",      
           "WT_pBM5", "WT_pBM5", "WT_pBM766", "WT_pBM766", "WT_pBM766")
labels <- group
y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) # eliminate compositional biases between libraries
plotMD(cpm(y, log=TRUE), column=1) # check that majority of transcripts are around 0
plotMDS(y, labels = labels)

design <- model.matrix(~ 0 + group) # construct model matrix
colnames(design) <- levels(as.factor(group))

y <- estimateDisp(y, design, robust = TRUE) # estimate dispersion
plotBCV(y) # plots the biological coefficient of variation for each gene
fit <- glmQLFit(y, design, robust=TRUE) # estimate values of the GLM coefficients for each gene
difex <- function(fit, contrast, outfile){
  qlf <- glmQLFTest(fit, contrast = contrast)
  tr <- glmTreat(fit, contrast=contrast, lfc=log2(1.2)) # narrow down the list of DE genes
  print(summary(decideTests(tr)))
  tags <- topTags(tr, n = 12379)
  tags <- tags$table
  tags$mrna <- row.names(tags)
  write.table(tags, file = paste0("Tag_seq/output_edgeR/", outfile, ".tsv"), 
              sep = "\t", row.names = F, quote = F)
}

WT_pBM5vsWT_pBM766 = makeContrasts(WT_pBM5 - WT_pBM766, levels = design)
difex(fit = fit, contrast = WT_pBM5vsWT_pBM766, outfile = "WT_pBM5vsWT_pBM766")

WT_pBM5vscsl4ph_pBM5 = makeContrasts(csl4ph_pBM5 - WT_pBM5, levels = design)
difex(fit = fit, contrast = WT_pBM5vscsl4ph_pBM5, outfile ='WT_pBM5vscsl4ph_pBM5')

WT_pBM5vscsl4ph_pBM766 = makeContrasts(csl4ph_pBM766 - WT_pBM5, levels = design)
difex(fit = fit, contrast = WT_pBM5vscsl4ph_pBM766, outfile ='WT_pBM5vscsl4ph_pBM766')

WT_pBM5vsenp1_1_pBM5 = makeContrasts(enp1_1_pBM5 - WT_pBM5, levels = design)
difex(fit = fit, contrast = WT_pBM5vsenp1_1_pBM5, outfile = 'WT_pBM5vsenp1_1_pBM5')

WT_pBM5vsenp1_1_pBM766 = makeContrasts(enp1_1_pBM766 - WT_pBM5, levels = design)
difex(fit = fit, contrast = WT_pBM5vsenp1_1_pBM766, outfile = 'WT_pBM5vsenp1_1_pBM766')

enp1_1_pBM5vsenp1_1_pBM766 = makeContrasts(enp1_1_pBM766 - enp1_1_pBM5, levels = design)
difex(fit = fit, contrast = enp1_1_pBM5vsenp1_1_pBM766, outfile = 'enp1_1_pBM5vsenp1_1_pBM766')

csl4ph_pBM5vscsl4ph_pBM766 = makeContrasts(csl4ph_pBM766 - csl4ph_pBM5, levels = design)
difex(fit = fit, contrast = csl4ph_pBM5vscsl4ph_pBM766, outfile ='csl4ph_pBM5vscsl4ph_pBM766')
