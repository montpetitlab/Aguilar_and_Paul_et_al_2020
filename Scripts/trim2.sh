for i in *trim.fastq; do
j=${i%.*}
~/sw/bbmap/bbduk.sh in=$i out=$j.trimbb.fastq \
ref=~/DATA2/test_tagseq/data/polyA.fa,/home/biplab/DATA2/test_tagseq/data/truseq.fa k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
#samtools view -Sb  $i > $j.bam 
#echo $j

done





#~/sw/bbmap/bbduk.sh in=BM_WT_pBM766_1_S32_L002_R1_001.trimm.fastq.gz out=BM_WT_pBM766_1_S32_L002_R1_001.trimbb.fastq.gz \
 # ref=~/DATA2/test_tagseq/data/polyA.fa,/home/biplab/DATA2/test_tagseq/data/truseq.fa k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
