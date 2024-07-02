#!/bin/bash

#### Description: Identify novel ORFs for Chinese hamster 
####    1. Run ORF-RATER
####    2. Filter ORF-RATER output
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie


# 1. Run ORF-RATER 
./scripts/orf_discovery/find_orfs.sh

# 2. Preprocess the data
source activate microprotein_r_env
Rscript ./scripts/orf_discovery/filter_orfrater.R
conda deactivate