#!/bin/bash
#### Description: Gene level differential translation counting   
####              1. Count the CHX Ribo-seq RPFs in individual 
####                 samples for the TS and NTS groups
####              2. Count the reads for individual 
####                 samples for the TS and NTS groups
####              3. Combine for import into R    
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

dt_dir=differential_translation
plastid_ref=$dt_dir/plastid_reference
chx_dir=sequencing/riboseq_chx/mapped/individual
rna_dir=sequencing/riboseq_se/mapped/individual

mkdir $dt_dir/cds_counts && count_dir=$_
# count the CHX riboseq data

for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files $chx_dir/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  --offset 12 \
  --min_length 28 \
  --max_length 31 \
  $plastid_ref/annotation_gene.positions \
  $count_dir/riboseq_$i
done

# count the RNASeq data
for i in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
do
cs count \
  --count_files $rna_dir/$i".bam" \
  --countfile_format BAM \
  --fiveprime \
  $plastid_ref/annotation_gene.positions \
  $count_dir/rnaseq_$i
done

# merge the RNASeq and Ribo-seq counts into a single table
# python scripts/differential_translation/combine_gene_cds_counts.py 