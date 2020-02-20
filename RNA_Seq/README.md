# Paul_et_al_2019
This github repository contain code for analysis of RNA-Seq and Tag-Seq data which is published as Paul et. al., 2019.

## RNA-Seq analysis:
### Alignment to the the genome by using hisat2 software, below the code for alignment
```
for f in `ls *.fastq | sed 's/_R[12].fastq//g' | sort -u`

do

hisat2 -p 4 --fr -x ~/DATA1/annotation/R64-1-1/Sequence/WholeGenomeFasta/genome --known-splicesite-infile ~/DATA1/annotation/R64-1-1/Annotation/Genes/yeast_splice_sites.txt -1 ${f}_R1.fastq -2 ${f}_R2.fastq | samtools view -Sb > ~/DATA2/hisat2_out/${f}.bam

done
```
### Counting reads by featurecount R package by running R scripts file 
### Differential expression analysis by using the R package DESeq2 by running R scripts file
### Log2FC data from DESeq2 were visualized by running R script for figures
### snoRNA analysis were done by using Deeptools in command lines. Scripts for generating snoRNA figure
convert bam file to bw file:
```
bamCoverage -b reads.bam -o coverage.bw
```
Computmatrix for plotting:
```
computeMatrix scale-regions -S <list of bigwig files> -R <bed file> --beforeRegionStartLength 200 --afterRegionStartLength 200 --skipZeros -o matrix.mat.gz
```
Code for metageneplot:

code for heatmap:
