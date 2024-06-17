#!/bin/bash

#### Description: Run the pipeline preprocess the NGS data for the study 
####    1. Download the data from NCBI SRA
####    2. Process the data
####    3. Determine the P-site for Ribo-seq data


#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie


# 1. Download the data 
./scripts/ngs/get_ngs_data.sh

# 2. Preprocess the data
./scripts/ngs/get_ngs_data.sh

# 3. Determine the P-site
./scripts/ngs/calc_psite.sh

# 4. Make coverage plots shown in the manuscript
./scripts/ngs/make_coverage_tracks.sh

