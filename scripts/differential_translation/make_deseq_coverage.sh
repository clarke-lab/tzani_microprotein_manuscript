#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source activate microprotein_process_env

mkdir differential_translation/coverage && cov_dir=$_

deseq_out_dir=differential_translation/deseq_output

rnaseq=sequencing/rnaseq_se/mapped/individual
riboseq=sequencing/riboseq_chx/mapped/individual

sample_names=("nts_r1" "nts_r2" "nts_r3" "nts_r4" "ts_r1" "ts_r2" "ts_r3" "ts_r4")

readarray -t rna_scale < <(cut -f 2 $deseq_out_dir/rna_scale_factors.txt)
readarray -t ribo_scale < <(cut -f 2 $deseq_out_dir/rpf_scale_factors.txt)

for i in "${!sample_names[@]}"; do

bamCoverage -b $rnaseq/${sample_names[i]}.bam \
     -o $cov_dir/"${sample_names[i]}"_rnaseq_scaled.bw \
     --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None  \
     --skipNonCoveredRegions --scaleFactor ${rna_scale[i]}

bamCoverage -b $riboseq/${sample_names[i]}.bam \
     -o $cov_dir/"${sample_names[i]}"_riboseq_scaled.bw \
     --outFileFormat bigwig -bs 1 -p 70 --normalizeUsing None  \
     --skipNonCoveredRegions --scaleFactor ${ribo_scale[i]}

done

conda deactivate