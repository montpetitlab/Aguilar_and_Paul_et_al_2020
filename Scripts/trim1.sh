for i in *fastq; do
j=${i%.*}
echo $j
trimmomatic SE $i $j.trim.fastq LEADING:12

done
