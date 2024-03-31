#!/bin/bash
#### Description: ORF-RATER analysis run via docker   
####                           
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# set the working directory
work_dir=/microprotein_analysis

mkdir $work_dir/orf_identification/orfrater && out_dir=$_

threads=60

prune_transcripts.py \
--inbed $work_dir/orf_identification/cgr.orfrater.annotation.reference.bed \
--summarytable $out_dir/tid_removal_summary.txt \
-p $threads \
--minlen 28 \
--maxlen 31 \
$work_dir/reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
$work_dir/sequencing/riboseq_chx/mapped/merged/riboseq_chx.bam \
$work_dir/sequencing/riboseq_nd/mapped/merged/riboseq_nd.bam \
--pseudogenes $work_dir/orf_identification/pseudogene_transcript_ids.txt \
--outbed $out_dir/transcripts.bed \
-v \
--force > $out_dir/1prune.log

make_tfams.py --force \
--inbed $out_dir/transcripts.bed \
--tfamstem $out_dir/tfams >   \
$out_dir/2make_tfams.log

# # echo "identifying orf candidates"
FILE=$out_dir/orf.h5
if test -f "$FILE"; then
    rm $FILE # need to remove this as problem with indexing if exists
fi    

find_orfs_and_types.py \
$work_dir/reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
--codons NTG \
-p $threads \
--force \
-v \
--tfamstem $out_dir/tfams \
--inbed $out_dir/transcripts.bed \
--orfstore $out_dir/orf.h5 > \
$out_dir/3find_ORF.log

mkdir $work_dir/orfrater_analysis/chx
psite_trimmed.py \
$work_dir/sequencing/riboseq_chx/mapped/merged/riboseq_chx.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $out_dir/chx \
--tallyfile tallies.txt \
-v \
--cdsbed $work_dir/orf_identification/cgr.orfrater.annotation.reference.bed \
-p $threads \
--force > $out_dir/psite_chx.log

mkdir $work_dir/orfrater_analysis/nd
psite_trimmed.py \
$work_dir/sequencing/riboseq_nd/mapped/merged/riboseq_nd.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $out_dir/nd \
--tallyfile tallies.txt \
--cdsbed $work_dir/orf_identification/cgr.orfrater.annotation.reference.bed \
-p $threads \
-v \
--force > $out_dir/psite_nd.log

mkdir $out_dir/harr
psite_trimmed.py \
$work_dir/sequencing/riboseq_harr/mapped/merged/riboseq_harr.bam \
--minrdlen 28 \
--maxrdlen 31 \
--subdir $out_dir/harr \
--tallyfile tallies.txt \
--cdsbed $work_dir/orf_identification/cgr.orfrater.annotation.reference.bed \
-p $threads \
-v \
--force > $out_dir/psite_har.log

regress_orfs.py \
$work_dir/sequencing/riboseq_harr/mapped/merged/riboseq_harr.bam \
--startonly \
--subdir $out_dir/harr \
--orfstore $out_dir/orf.h5 \
--inbed $out_dir/transcripts.bed \
-p 50 \
-v \
--startcount 1 \
--force > $out_dir/regress_start.log

regress_orfs.py \
$work_dir/sequencing/riboseq_chx/mapped/merged/riboseq_chx.bam \
--subdir $out_dir/chx \
--orfstore $out_dir/orf.h5 \
--inbed $out_dir/transcripts.bed \
--restrictbystarts $out_dir/harr \
--startcount 1 \
-p $threads \
-v \
--force > $out_dir/chx_regress_stop.log

regress_orfs.py \
$work_dir/sequencing/riboseq_nd/mapped/merged/riboseq_nd.bam \
--subdir $out_dir/nd \
--orfstore $out_dir/orf.h5 \
--inbed $out_dir/transcripts.bed \
--restrictbystarts $out_dir/harr \
--startcount 1 \
-p $threads \
--force > $out_dir/nodrug_regress_stop.log

rate_regression_output.py \
$out_dir/harr \
$out_dir/chx \
$out_dir/nd \
--orfstore $out_dir/orf.h5 \
--ratingsfile $out_dir/orfratings.h5 \
-p $threads \
-v \
--CSV $out_dir/rate_regression.csv \
--force > $out_dir/rate.regression.log

 make_orf_bed.py \
 --inbed $out_dir/transcripts.bed \
 --ratingsfile $out_dir/orfratings.h5 \
 --minlen 5 \
 --force \
 --outbed $out_dir/orfrater_predictions.reference.bed