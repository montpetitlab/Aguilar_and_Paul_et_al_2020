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
## Taq-Seq analysis:
