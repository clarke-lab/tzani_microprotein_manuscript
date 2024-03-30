#!/bin/bash

# Description: Download the RIBO-seq and RNA-seq data 
#              from ENA using ffq 
#
# Credit for ffq:
# Gálvez-Merchán, Á., et al. (2022). Metadata retrieval from sequence databases with ffq. 
# bioRxiv 2022.05.18.492548.
# 
#
# Author: Clarke Lab. NIBRT 

# 1. Total RNASeq (Single-end)
mkdir -p $PWD/sequencing/rnaseq_se/raw_data && rna_dir=$_

rna_samples=("SRR16796972" "SRR16796971" "SRR16796960" "SRR16796949" 
             "SRR16796946" "SRR16796945" "SRR16796944" "SRR16796943")

for sample in ${rna_samples[@]}; do
    (cd $rna_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

echo "RNA-seq downloaded"

# 2. CHX riboseq
mkdir -p sequencing/riboseq_chx/raw_data && chx_dir=$_

chx_samples=("SRR16796955" "SRR16796954" "SRR16796953" "SRR16796952" 
             "SRR16796951" "SRR16796950" "SRR16796948" "SRR16796947")

for sample in ${chx_samples[@]}; do
    (cd $chx_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done    

echo "CHX Ribo-seq downloaded"

# 3. Harr riboseq
mkdir -p $PWD/sequencing/riboseq_harr/raw_data && harr_dir=$_

harr_samples=("SRR16796942" "SRR16796941" "SRR16796970" "SRR16796969" 
              "SRR16796968" "SRR16796967" "SRR16796966" "SRR16796965")

for sample in ${harr_samples[@]}; do
    (cd $harr_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

echo "Harr Ribo-seq downloaded"

# 4. No drug riboseq
mkdir -p $PWD/sequencing/riboseq_nd/raw_data && nd_dir=$_

nd_samples=("SRR16796964" "SRR16796963" "SRR16796962" "SRR16796961" 
            "SRR16796959" "SRR16796958" "SRR16796957" "SRR16796956")

for sample in ${nd_samples[@]}; do
    (cd $nd_dir && ffq --ftp $sample | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O -s)
done

echo "ND Ribo-seq downloaded"
echo "Done"
# end