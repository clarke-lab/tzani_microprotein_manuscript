#!/bin/bash

# Description: Prepocesses raw sequencing data for   
#              1. Ribo-seq
#                  a) adapter trimming
#                  b) ncRNA contamination
#                  c) Read length based on phasing
#              2. RNA-seq  
#                  a) adapter trimming
# 
# Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# num_cores=$(getconf_NPROCESSORS_ONLN)
# set -x

mkdir reference_genome
cat data/reference_genome_files.txt | parallel -j 4 wget -P reference_genome {}
gunzip reference_genome/*.gz

# build the index
if ! [ -d reference_genome/star_index ]; then
mkdir reference_genome/star_index && star_index=$_
STAR --runThreadN 32 \
     --runMode genomeGenerate \
     --sjdbOverhang 31 \
     --genomeChrBinNbits 16 \
     --genomeDir $star_index \
     --genomeFastaFiles reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna \
     --sjdbGTFfile reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf
else 
  star_index=reference_genome/star_index
fi

# Ribo-seq
# build indexes for contamination filtering
srna_fasta=data/small_rna_fasta

mkdir reference_genome/srna_bt_index && srna_index=$_

# a) rRNA
if ! [ -f $srna_index/rRNA.1.ebwt ]; then
    bowtie-build $srna_fasta/rna_central_v22_rRNA.fasta $srna_index/rRNA
fi

# b) tRNA
if ! [ -f $srna_index/tRNA.1.ebwt ]; then
    bowtie-build $srna_fasta/rna_central_v22_tRNA.fasta $srna_index/tRNA
fi

# c) snoRNA
if ! [ -f $srna_index/snoRNA.1.ebwt ]; then
    bowtie-build $srna_fasta/rna_central_v22_snoRNA.fasta $srna_index/snoRNA
fi

# trim and filter contamination from each sample 
# for each type of Ribo-seq data
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
   mkdir -p sequencing/$seqtype/preprocessed/trimmed
   mkdir  sequencing/$seqtype/preprocessed/rRNA_filter
   mkdir  sequencing/$seqtype/preprocessed/tRNA_filter
   mkdir  sequencing/$seqtype/preprocessed/snoRNA_filter
   mkdir sequencing/$seqtype/preprocessed/complete
   mkdir -p sequencing/$seqtype/mapped/individual

while read -ra a ;
  do
    # run cutadapt; note chx data has adapter removed by sequencing provider / -m flag used to remove 
    # short reads post trimming for HARR/ND

    if ! [ -f sequencing/$seqtype/preprocessed/trimmed/${a[0]} ]; then   
      if [ $seqtype == 'riboseq_chx' ]
      then
        cutadapt  --report=full -a AGATCGGAAGAGCACACGTCT -j 50 -m 20 \
        -o sequencing/$seqtype/preprocessed/trimmed/${a[0]} sequencing/$seqtype/raw_data/${a[0]}
      else
        cutadapt  --discard-untrimmed -m 20 --report=full -a AGATCGGAAGAGCACACGTCT -j 50 --minimum-length 1 \
        -o sequencing/$seqtype/preprocessed/trimmed/${a[0]} sequencing/$seqtype/raw_data/${a[0]}
      fi
    fi
  
# Perform action here
  
    # a) filter rRNA
    if ! [ -f sequencing/$seqtype/preprocessed/rRNA_filter/${a[1]}_unaligned.fq ]; then
        gunzip -c sequencing/$seqtype/preprocessed/trimmed/${a[0]} | \
        bowtie -v 2 -p 50 -l 20 -norc \
        $srna_index/rRNA \
        -q - \
        sequencing/$seqtype/preprocessed/rRNA_filter/${a[1]}_aligned.fq \
        --un sequencing/$seqtype/preprocessed/rRNA_filter/${a[1]}_unaligned.fq
    fi

    # b) filter snoRNA
    if ! [ -f sequencing/$seqtype/preprocessed/snoRNA_filter/${a[1]}_unaligned.fq ]; then
        bowtie -v 2 -p 50 -l 20 -norc \
        $srna_index/snoRNA \
        -q sequencing/$seqtype/preprocessed/rRNA_filter/${a[1]}_unaligned.fq \
        sequencing/$seqtype/preprocessed/snoRNA_filter/${a[1]}_aligned.fq \
        --un sequencing/$seqtype/preprocessed/snoRNA_filter/${a[1]}_unaligned.fq
    fi

    # c) filter tRNA
    if ! [ -f sequencing/$seqtype/preprocessed/tRNA_filter/${a[1]}_unaligned.fq ]; then
       bowtie -v 2 -p 50 -l 20 -norc \
        $srna_index/tRNA \
        -q sequencing/$seqtype/preprocessed/snoRNA_filter/${a[1]}_unaligned.fq \
        sequencing/$seqtype/preprocessed/tRNA_filter/${a[1]}_aligned.fq \
        --un sequencing/$seqtype/preprocessed/tRNA_filter/${a[1]}_unaligned.fq
    fi

    # d) remove unecessary files
    # rm sequencing/$seqtype/preprocessed/*_filter/${a[1]}.*.Aligned.out.sam
    # rm sequencing/$seqtype/preprocessed/*_filter/${a[1]}.*.Log.progress.out
    # rm sequencing/$seqtype/preprocessed/*_filter/${a[1]}.*.SJ.out.tab
  

    # e) filter the read lengths based on phasing qc (28nt to 31nt)
    if ! [ -f sequencing/$seqtype/preprocessed/complete/${a[1]}.fastq ]; then
      awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 28 && length(seq) <= 31) {print header, seq, qheader, qseq}}' \
      < sequencing/$seqtype/preprocessed/tRNA_filter/${a[1]}_unaligned.fq > sequencing/$seqtype/preprocessed/complete/${a[1]}.fastq
    fi
    # map the retained RPFs to the PICR reference genome

  echo "mapping to genome"

   STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMismatchNmax 2 \
    --genomeDir $star_index \
    --readFilesIn sequencing/$seqtype/preprocessed/complete/${a[1]}.fastq \
    --outFileNamePrefix sequencing/$seqtype/mapped/individual/${a[1]} \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv sequencing/$seqtype/mapped/individual/${a[1]}"Aligned.sortedByCoord.out.bam" \
    sequencing/$seqtype/mapped/individual/${a[1]}".bam"

    # index bam
    samtools index sequencing/$seqtype/mapped/individual/${a[1]}".bam"

  done < data/seq_files/"$seqtype".txt
done

# merge and map the files for each treatment for ORF-RATER
for seqtype in riboseq_harr riboseq_nd riboseq_chx
  do
    mkdir sequencing/$seqtype/mapped/merged

    # merge the RPF for each Ribo-seq type 
    cat sequencing/$seqtype/preprocessed/complete/*.fastq > \
    sequencing/$seqtype/mapped/merged/$seqtype.fastq

    STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --outFilterMismatchNmax 2 \
    --seedSearchStartLmaxOverLread .5 \
    --genomeDir $star_index \
    --readFilesIn sequencing/$seqtype/mapped/merged/$seqtype.fastq \
    --outFileNamePrefix sequencing/$seqtype/mapped/merged/$seqtype \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    # rename
    mv sequencing/$seqtype/mapped/merged/$seqtype"Aligned.sortedByCoord.out.bam" \
    sequencing/$seqtype/mapped/merged/$seqtype".bam"

    # index bam
    samtools index sequencing/$seqtype/mapped/merged/$seqtype.bam
done

# 2. RNASeq data 
mkdir -p sequencing/rnaseq_se/preprocessed/complete
mkdir -p sequencing/rnaseq_se/mapped/individual

# individual samples
while read -ra a ;
 do
     # run cutadapt
     cutadapt  -q 30 -m 20 --report=full \
     -a AGATCGGAAGAGCACACGTCT -j 32 \
     --minimum-length 1 \
     -o sequencing/rnaseq_se/preprocessed/complete/${a[0]} sequencing/rnaseq_se/raw_data/${a[0]}

    STAR \
       --outSAMtype BAM SortedByCoordinate \
       --runThreadN 16 \
       --outFilterMismatchNmax 2 \
       --seedSearchStartLmaxOverLread .5 \
       --genomeDir $star_index \
       --readFilesIn  sequencing/rnaseq_se/preprocessed/complete/${a[0]} \
       --outFileNamePrefix sequencing/rnaseq_se/mapped/individual/${a[1]} \
       --outFilterMultimapNmax 1 \
       --outFilterMatchNmin 16 \
       --alignEndsType EndToEnd \
       --readFilesCommand zcat \
       --outMultimapperOrder Random \
       --outSAMattributes All

     mv sequencing/rnaseq_se/mapped/individual/${a[1]}Aligned.sortedByCoord.out.bam sequencing/rnaseq_se/mapped/individual/${a[1]}.bam
     samtools index sequencing/rnaseq_se/mapped/individual/${a[1]}.bam

done < data/seq_files/rnaseq_se.txt

# map the merged RNA-seq reads
mkdir -p sequencing/rnaseq_se/mapped/merged 

zcat sequencing/rnaseq_se/preprocessed/complete/* > \
sequencing/rnaseq_se/mapped/merged/rnaseq_se.fastq 

STAR \
    --outFilterType BySJout \
    --runThreadN 16 \
    --seedSearchStartLmaxOverLread .5 \
    --outFilterMismatchNmax 2 \
    --genomeDir $star_index \
    --readFilesIn sequencing/rnaseq_se/mapped/merged/rnaseq_se.fastq \
    --outFileNamePrefix sequencing/rnaseq_se/mapped/merged/rnaseq_se \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All

    mv sequencing/rnaseq_se/mapped/merged/rnaseq_seAligned.sortedByCoord.out.bam \
    sequencing/rnaseq_se/mapped/merged/rnaseq_se.bam

    # index bam
    samtools index sequencing/rnaseq_se/mapped/merged/rnaseq_se.bam
# end

# count the reads and determine read length distribution 

echo "Determing counts and read distributions"
./scripts/ngs/calc_read_stats.sh