setwd("~/DATA2/hisat2_out/")
library("Rsubread")
myFiles <- list.files(pattern="*bam")
myFiles_dT <- myFiles[grep("mR|-m[123]", myFiles)]
myFiles_rd <- myFiles[grep("dpl|-d[123]|CA_WT-[12].bam|CA_Rrp6|CA_Dis3", myFiles)]
wt_vs_wt37c <- featureCounts(myFiles_rd,
                             
                             # annotation
                             annot.inbuilt=NULL,
                             annot.ext="/home/biplab/project_exosome/input_files/gtf_files/Saccharomyces_cerevisiae.R64-1-1.85.gtf",
                             isGTFAnnotationFile=TRUE,
                             GTF.featureType="exon",
                             GTF.attrType="gene_id",
                             chrAliases=NULL,
                             
                             # level of summarization
                             useMetaFeatures=FALSE,
                             
                             # overlap between reads and features
                             allowMultiOverlap=TRUE,
                             minOverlap=1,
                             fracOverlap=0,
                             fracOverlapFeature=0,
                             largestOverlap=FALSE,
                             readExtension5=0,
                             readExtension3=0,
                             read2pos=NULL,
                             
                             # multi-mapping reads
                             countMultiMappingReads=FALSE,
                             
                             # fractional counting
                             fraction=FALSE,
                             
                             # long reads
                             isLongRead=FALSE,
                             
                             # read filtering
                             minMQS=0,
                             splitOnly=FALSE,
                             nonSplitOnly=FALSE,
                             primaryOnly=FALSE,
                             ignoreDup=FALSE,
                             
                             # strandness
                             strandSpecific=2,
                             
                             # exon-exon junctions
                             juncCounts=FALSE,
                             genome=NULL,
                             
                             # parameters specific to paired end reads
                             isPairedEnd=TRUE,
                             requireBothEndsMapped=TRUE,
                             checkFragLength=TRUE,
                             minFragLength=20,
                             maxFragLength=600,
                             countChimericFragments=FALSE,	
                             autosort=FALSE,
                             
                             # number of CPU threads
                             nthreads=4,
                             
                             # read group
                             byReadGroup=FALSE,
                             
                             # miscellaneous
                             maxMOp=10,
                             reportReads=NULL,
                             tmpDir=".",
                             verbose=FALSE)

wt_mutants <- data.frame(wt_vs_wt37c$annotation,wt_vs_wt37c$counts,stringsAsFactors=FALSE)
write.table(wt_mutants, file = "~/project_exosome/input_files/expression_data/fc_mRNA_snRNA_rd.csv", sep = "\t", quote = FALSE, col.names = TRUE)

wt_vs_wt37c <- featureCounts(myFiles_dT,
                             
                             # annotation
                             annot.inbuilt=NULL,
                             annot.ext="/home/biplab/project_exosome/input_files/gtf_files/Saccharomyces_cerevisiae.R64-1-1.85.gtf",
                             isGTFAnnotationFile=TRUE,
                             GTF.featureType="exon",
                             GTF.attrType="gene_id",
                             chrAliases=NULL,
                             
                             # level of summarization
                             useMetaFeatures=FALSE,
                             
                             # overlap between reads and features
                             allowMultiOverlap=TRUE,
                             minOverlap=1,
                             fracOverlap=0,
                             fracOverlapFeature=0,
                             largestOverlap=FALSE,
                             readExtension5=0,
                             readExtension3=0,
                             read2pos=NULL,
                             
                             # multi-mapping reads
                             countMultiMappingReads=FALSE,
                             
                             # fractional counting
                             fraction=FALSE,
                             
                             # long reads
                             isLongRead=FALSE,
                             
                             # read filtering
                             minMQS=0,
                             splitOnly=FALSE,
                             nonSplitOnly=FALSE,
                             primaryOnly=FALSE,
                             ignoreDup=FALSE,
                             
                             # strandness
                             strandSpecific=2,
                             
                             # exon-exon junctions
                             juncCounts=FALSE,
                             genome=NULL,
                             
                             # parameters specific to paired end reads
                             isPairedEnd=TRUE,
                             requireBothEndsMapped=TRUE,
                             checkFragLength=TRUE,
                             minFragLength=20,
                             maxFragLength=600,
                             countChimericFragments=FALSE,	
                             autosort=FALSE,
                             
                             # number of CPU threads
                             nthreads=4,
                             
                             # read group
                             byReadGroup=FALSE,
                             
                             # miscellaneous
                             maxMOp=10,
                             reportReads=NULL,
                             tmpDir=".",
                             verbose=FALSE)

wt_mutants <- data.frame(wt_vs_wt37c$annotation,wt_vs_wt37c$counts,stringsAsFactors=FALSE)
write.table(wt_mutants, file = "~/project_exosome/input_files/expression_data/fc_mRNA_snRNA_dT.csv", sep = "\t", quote = FALSE, col.names = TRUE)




