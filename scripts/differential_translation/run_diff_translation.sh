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

mkdir -p $dt_dir/quantitation/cds_
./scripts/differential_translation/calculate_rpkm.sh

# mkdir $dt_dir/cds_counts
./scripts/differential_translation/count_cds_reads.sh

Rscript ./scripts/run_deseq2.R