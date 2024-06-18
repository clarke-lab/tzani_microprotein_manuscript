#!/bin/bash
#### Description: Identify canonical and novel ORFs where gene expression and/or 
#### translation changes upon temperature reduction  
####              1. Make a plastid reference for counting
####              2. Count against the reference 
####              3. Run DESeq2 and the deltaTE method
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.

mkdir differential_translation 

# 1. Make a reference to count with Plastid
./scripts/differential_translation/make_plastid_reference.sh

# 2. Count RPFs and RNA against the reference 
./scripts/differential_translation/count_cds_reads.sh

# 3. Run the DESeq with the delta TE method
source activate microprotein_r_env
Rscript ./scripts/run_deseq2.R
conda deactivate