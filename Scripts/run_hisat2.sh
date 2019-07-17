for f in `ls *.fastq | sed 's/_R[12].fastq//g' | sort -u`

do

hisat2 -p 4 --fr -x ~/DATA1/annotation/R64-1-1/Sequence/WholeGenomeFasta/genome --known-splicesite-infile ~/DATA1/annotation/R64-1-1/Annotation/Genes/yeast_splice_sites.txt -1 ${f}_R1.fastq -2 ${f}_R2.fastq | samtools view -Sb > ~/DATA2/hisat2_out/${f}.bam

done
