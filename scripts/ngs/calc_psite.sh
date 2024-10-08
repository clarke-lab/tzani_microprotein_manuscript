#!/bin/bash
#### Description: calculate offsets for Ribo-seq   
####              data using Plastid              
####              the RNA-seq data is also analyse for comparison 
####              this analysis is carried out for 
####              1. individual 
####              2. merged samples   
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source activate microprotein_process_env
mkdir sequencing/plastid_analysis && plastid_dir=$_

mkdir $plastid_dir/annotation
# make a Plastid annotation for the known CDS start positions
# across the genome
metagene generate $plastid_dir/annotation/cgr_orfs \
--landmark cds_start \
--annotation_files reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf

# 1. individual sample psite offsets
for seqtype in riboseq_chx riboseq_harr riboseq_nd
  do
  
  mkdir -p $plastid_dir/individual_files/$seqtype/offset
  mkdir $plastid_dir/individual_files/$seqtype/periodicity

  while read -ra a ;
    do
      psite $plastid_dir/annotation/cgr_orfs_rois.txt  \
      $plastid_dir/individual_files/$seqtype/offset/${a[1]} \
      --min_length 28 \
      --max_length 31 \
      --require_upstream \
      --count_files sequencing/$seqtype/mapped/individual/${a[1]}".bam"

      phase_by_size $plastid_dir/annotation/cgr_orfs_rois.txt \
      $plastid_dir/individual_files/$seqtype/periodicity/${a[1]} \
      --count_files sequencing/$seqtype/mapped/individual/${a[1]}".bam" \
      --fiveprime  \
      --offset 12 \
      --codon_buffer 5 \
      --min_length 28 \
      --max_length 31

    done < data/seq_files/"$seqtype".txt
done

mkdir $plastid_dir/merged
# 2. merged samples psite offset
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
  do
  
  psite $plastid_dir/annotation/cgr_orfs_rois.txt \
  $plastid_dir/merged/$seqtype \
  --min_length 28 \
  --max_length 31 \
  --require_upstream \
  --min_counts 10 \
  --count_files sequencing/$seqtype/mapped/merged/$seqtype".bam"

done

# Determine the Proportion of reads found to be in frame at the 12 nt offset
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
do
  phase_by_size $plastid_dir/annotation/cgr_orfs_rois.txt \
  $plastid_dir/merged/$seqtype \
  --count_files sequencing/$seqtype/mapped/merged/$seqtype".bam" \
  --fiveprime \
  --offset 12 \
  --codon_buffer 5 \
  --min_length 28 \
  --max_length 31
done
conda deactivate