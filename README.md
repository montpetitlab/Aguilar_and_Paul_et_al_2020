# Aguilar and Paul et al 2020
This github repository contain code for analysis of ribo-purified, oligo-dT enriched, and 3' tag-seq rna-sequencing data from Aguilar and Paul et. al., 2020.

Each folder contains a snakefile that coordinates preprocessing, analysis, and visualization of each dataset. 

To run each pipeline, install miniconda. Then run:

```
conda create -n ap2020 snakemake-minimal
conda activate ap2020
```

Then navigate into a directory and run the snakefile:

```
snakemake --use-conda
```
