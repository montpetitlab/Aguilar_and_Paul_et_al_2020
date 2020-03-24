# Paul_et_al_2019
This github repository contain code for analysis of Tag-Seq data which is published as Paul et. al., 2019.

In order to run these scripts in your computer, you need to have following software installed in your computer.
1. Trimmomatic
2. BBTools
3. Salmon
4. Deeptools
5. R
6. R package(DESeq2, Rsubread, ggplot2, dplyr)

## Tag-Seq analysis:
### Pre-processing of fastq files
fastq files are pre-processed by trimmomatic and bbtools using following shell scripts:
```
for i in *fastq; do
j=${i%.*}
echo $j
trimmomatic SE $i $j.trim.fastq LEADING:12

done
```

```
for i in *trim.fastq; do
j=${i%.*}
~/sw/bbmap/bbduk.sh in=$i out=$j.trimbb.fastq \
ref=~/DATA2/test_tagseq/data/polyA.fa,/home/biplab/DATA2/test_tagseq/data/truseq.fa k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20

done
```
### Convert gtf file to bed file
```
awk '{if($3=="transcript")print $1 "\t" $4 "\t" $5 "\t" $10}' new_gtf_file.gtf | sed 's/"//g' | sed 's/;//g'
```
### Generating multifasta file for all transcripts, which includes mRNA, snoRNA and pervasive transcripts.
```
bedtools getfasta -fi genome.fa -bed all_transcript1.bed -name > all_transcript1.fa
```
### Counting read for features were calculated by Salmon, which use fastq file as input and output count files
```
for i in *trimbb.fastq; do
j=${i%.*}
salmon quant -i ~/DATA2/Tag_Seq_Analysis/all_transcript1 -l A -r $i --writeUnmappedNames -o ~/DATA2/Tag_Seq_Analysis/pervasive_salmon_outy/$j

done
```

### Differential expression analysis by using the R package edgeR by running edgeR.R scripts file
### Log2FC data from edgeR were visualized by running R script for figures
### snoRNA analysis were done by using Deeptools in command lines. Scripts for generating snoRNA figure
Aligning reads to genome:
```
for f in `ls *.fastq | sed 's/_R[12].fastq//g' | sort -u`

do

hisat2 -p 4 --fr -x ~/DATA1/annotation/R64-1-1/Sequence/WholeGenomeFasta/genome --known-splicesite-infile ~/DATA1/annotation/R64-1-1/Annotation/Genes/yeast_splice_sites.txt -1 ${f}_R1.fastq -2 ${f}_R2.fastq | samtools view -Sb > ~/DATA2/hisat2_out/${f}.bam

done
```

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
