#!/bin/bash
#### Description: Prepocesses raw sequencing data for   
####              1. Ribo-seq
####                  a) adapter trimming
####                  b) ncRNA contamination
####                  c) Read length based on phasing
####              2. RNA-seq  
####                  a) adapter trimming
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.

mkdir differential_translation && dt_dir=$_
conda activate microprotein_process_env
# Make a reference to count with Plastid
./scripts/differential_translation/make_plastid_reference.sh

mkdir -p $dt_dir/quantitation/transcript_cds_rpkm
./scripts/differential_translation/calculate_rpkm.sh

mkdir $dt_dir/quantitation/gene_cds_counts
./scripts/differential_translation/calculate_gene_cds_counts.sh

Rscript ./scripts/run_deseq2.R