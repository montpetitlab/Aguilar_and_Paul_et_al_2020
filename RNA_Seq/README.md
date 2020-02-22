# Paul_et_al_2019
This github repository contain code for analysis of RNA-Seq and Tag-Seq data which is published as Paul et. al., 2019.

In order to run these scripts in your computer, you need to have following software installed in your computer.
1. Samtools
2. Hisat2
3. Deeptools
4. R
5. R package(DESeq2, Rsubread, ggplot2, dplyr)

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
bamCoverage -b <file.bam> -o <file.bw>
```
bw files that are generated in previous step are used as input for Computmatrix command of deeptool, which devide every genes into equal number of bins and give counts of reads for each bins in a matrix form. Shells commands for this is given below:
```
computeMatrix scale-regions -S <list_of_bigwig_files.bw> -R <bed_file.bed> --beforeRegionStartLength 200 --afterRegionStartLength 200 --skipZeros -o matrix.mat.gz
```
The matrix.mat.gz file that is generated from computeMatrix command is directly used as input to generate metageneplot. Deeptools command for matageneplots are given below:
```
plotProfile -m matrix.mat.gz --perGroup -out ExampleProfile2.png
```
The matrix.mat.gz file that is generated from computeMatrix command is also use directly as input to generate heatmap. Deeptool command for matageneplots are given below:
```
plotHeatmap -m matrix.mat.gz -out ExampleHeatmap1.png
```
