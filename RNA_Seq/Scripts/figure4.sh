bamCoverage -bs 1 -b HI.3202.005.Index_3.BM_Ctrl-NAB2-PrA-CSL4-URA-d2.sorted.bam -o ctrl.bw
bamCoverage -bs 1 -b HI.3202.005.Index_20.BM_shift90-NAB2-PrA-CSL4pH-URA-37C-d1.sorted.bam -o csl4_ph.bw
bamCoverage -bs 1 -b HI.3758.002.Index_14.BM-MO_dpl-enp1-37C-90-4-1.sorted.bam -o enp1.bw

computeMatrix scale-regions -R ~/transcriptomics/snoRNA_BED6.bed -S ctrl.bw csl4_ph.bw enp1.bw -b 200 -a 200 --skipZeros -o matrix.mat.gz

plotHeatmap -m matrix.mat.gz -out ExampleHeatmap.png
