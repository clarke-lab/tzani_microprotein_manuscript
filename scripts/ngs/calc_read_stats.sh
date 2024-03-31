#!/bin/bash
#### Description: Count the reads from each stage of preprocessing for the Ribo-seq data   
####              1. Raw reads
####              2. Trimmed
####              3. rRNA, tRNA, snoRNA
####              3. Final RPFs 28-31nt   
#### 
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# 1. count the raw reads
mkdir -p sequencing/stats/read_counts && count_dir=$_
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for f in sequencing/$seqtype/raw_data/*fastq.gz; 
        do echo -ne "$f\t"  && echo $(($(zcat $f | echo $((`wc -l`/4))))); 
    done > $count_dir/$seqtype.raw.counts && sed  -i '1i file\traw_read_number' $count_dir/$seqtype.raw.counts
done

# 2. count the reads removed by trimming 
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for f in sequencing/$seqtype/preprocessed/trimmed/*fastq.gz; 
        do echo -ne "$f\t"  && echo $(($(zcat $f | echo $((`wc -l`/4))))); 
    done > $count_dir/$seqtype.trimmed.counts && sed  -i '1i file\traw_read_number' $count_dir/$seqtype.trimmed.counts
done

# 3. count the reads remaining after filtering non-coding RNA
for seqtype in riboseq_chx riboseq_harr riboseq_nd
do
    for step in rRNA_filter snoRNA_filter tRNA_filter
    do
        for f in sequencing/"$seqtype"/preprocessed/"$step"/*_unaligned.fq; 
            do echo -ne "$f\t"  && echo $(($(cat $f | echo $((`wc -l`/4))))); 
        done > $count_dir/$seqtype.$step.counts && sed  -i '1i file\traw_read_number' $count_dir/$seqtype.$step.counts
    done
done

# 4. count the reads removed after preprocessing
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
do
    for f in sequencing/$seqtype/preprocessed/complete/*fastq*;
        do 
        if [ $seqtype = "rnaseq_se" ]; then
            echo -ne "$f\t"  && echo $(($(zcat $f | echo $((`wc -l`/4))))); 
        else
            echo -ne "$f\t"  && echo $(($(cat $f | echo $((`wc -l`/4))))); 
        fi

    done > $count_dir/$seqtype.final.counts && sed  -i '1i file\traw_read_number' $count_dir/$seqtype.final.counts
done

mkdir sequencing/stats/length_distribution && dist_dir=$_

# 2. count read length distributions
for seqtype in riboseq_chx riboseq_harr riboseq_nd rnaseq_se
do  
    if [ $seqtype = "rnaseq_se" ]; then
                cat sequencing/$seqtype/mapped/merged/rnaseq_se.fastq | \
                awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $dist_dir/$seqtype.txt
    else
        for sample in nts_r1 nts_r2 nts_r3 nts_r4 ts_r1 ts_r2 ts_r3 ts_r4
        do
                cat sequencing/$seqtype/preprocessed/tRNA_filter/"$sample"_unaligned.fq | \
                awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $dist_dir/$seqtype.$sample.txt
        done
    fi
done
