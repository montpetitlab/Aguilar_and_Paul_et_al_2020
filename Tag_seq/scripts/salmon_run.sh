for i in *trimbb.fastq; do
j=${i%.*}
salmon quant -i ~/DATA2/Tag_Seq_Analysis/all_transcript1 -l A -r $i --writeUnmappedNames -o ~/DATA2/Tag_Seq_Analysis/pervasive_salmon_outy/$j
#htseq-count -f bam -s no $i /home/biplab/gtf_files/Saccharomyces_cerevisiae.R64-1-1.85.gtf > $j.htseq_out
#samtools view -Sb  $i > $j.bam 
#echo $j

done
