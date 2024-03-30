#!/bin/bash
#### Description: Make coverage tracks for the manuscript figures
####             1. Merged Data for Figure 2B
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

### 1 merged transcriptome tracks for new annotations 

# # 1. Merged Data for Figure 2B
# ## Demonstrate the coverage at transcript-level
# ## a) CHX full coverage
# ## b) CHX p-site offset
# ## c) CHX a-site offset
# ## d) HARR-ND p-site

mkdir -p sequencing/coverage_tracks/merged && merged_dir=$_

for seqtype in riboseq_harr riboseq_chx riboseq_nd
do
  if ! [ -f  sequencing/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam.bai" ]; then
    ## sort and index transcriptome aligned BAMs
    samtools sort \
    sequencing/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.bam" \
    -o sequencing/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"

    samtools index \
    sequencing/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam"
  fi

  make_wiggle \
     --count_files sequencing/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
     --countfile_format BAM \
     --fiveprime \
     --offset 12 \
     --output_format bedgraph \
     -o $merged_dir/$seqtype.psite.transcriptome.nonorm.bedgraph

done

for seqtype in riboseq_chx riboseq_harr riboseq_nd riboseq_chx
 do
     make_wiggle \
     --count_files data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
     --countfile_format BAM \
     --fiveprime \
     --offset 12 \
     --output_format bedgraph \
     -o $merged_dir/$seqtype.psite.transcriptome.nonorm.bedgraph 
    
done

# full coverage CHX track
make_wiggle \
     --count_files sequencing/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sorted.bam \
     --countfile_format BAM \
     --fiveprime \
     --offset 0 \
     --output_format bedgraph \
     -o $merged_dir/riboseq_chx.full.transcriptome.nonorm.bedgraph 

bamCoverage -b sequencing/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.sorted.bam \
-o $merged_dir/chx.fullcov.transcriptome.wig \
--outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing CPM --smoothLength 25

# extract selected transcripts plotting

extract_transcript_coverage() {
    local transcript="$1"
    local in_file="$2"   
    local out_dir="$3"
    local out_file="$4"
 
    if [ ! -f "$in_file" ]; then
        echo "Error: File '$in_file' not found."
        return 1
    fi

    grep "$transcript" "$in_file" > $out_dir/"$transcript"_"$out_file"
}

chx_transcriptome_p_fw="$merged_dir"/riboseq_chx.psite.transcriptome.nonorm.bedgraph_fw.wig
chx_transcriptome_p_rc=$merged_dir/riboseq_chx.psite.transcriptome.nonorm.bedgraph_rc.wig

chx_transcriptome_f_fw=$merged_dir/chx.fullcov.transcriptome.wig
harr_transcriptome_p_fw=$merged_dir/riboseq_harr.psite.transcriptome.nonorm.bedgraph_fw.wig
harr_transcriptome_p_rc=$merged_dir/riboseq_harr.psite.transcriptome.nonorm.bedgraph_rc.wig
nd_transcriptome_p_fw=$merged_dir/riboseq_nd.psite.transcriptome.nonorm.bedgraph_fw.wig
nd_transcriptome_p_rc=$merged_dir/riboseq_nd.psite.transcriptome.nonorm.bedgraph_rc.wig

# Fig2b
mkdir $merged_dir/fig2b && fig2b=$_

extract_transcript_coverage 'XM_027423276.2' "$chx_transcriptome_p_fw" $fig2b 'chx_p_transcriptome.wig'
extract_transcript_coverage 'XM_027423276.2' $chx_transcriptome_f_fw $fig2b 'chx_f_transcriptome.wig'
extract_transcript_coverage 'XM_027423276.2' $harr_transcriptome_p_fw $fig2b 'harr_p_transcriptome.wig'
extract_transcript_coverage 'XM_027423276.2' $nd_transcriptome_p_fw $fig2b 'nd_p_transcriptome.wig'

# Fig3d
mkdir $merged_dir/fig3d && fig3d=$_

extract_transcript_coverage 'XM_027392226.1' $chx_transcriptome_p_fw $fig3d chx_p_transcriptome.wig
extract_transcript_coverage 'XM_027392226.1' $chx_transcriptome_f_fw $fig3d chx_f_transcriptome.wig
extract_transcript_coverage 'XM_027392226.1' $harr_transcriptome_p_fw $fig3d harr_p_transcriptome.wig
extract_transcript_coverage 'XM_027392226.1' $nd_transcriptome_p_fw $fig3d nd_p_transcriptome.wig






bamCoverage -b data/riboseq_chx/mapped/merged/riboseq_chx.bam -o $merged_dir/chx.fullcov.genome.bedgraph \
  --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing CPM --smoothLength 25

for seqtype in riboseq_harr riboseq_nd riboseq_harr riboseq_chx
 do
     make_wiggle \
     --count_files data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
     --countfile_format BAM \
     --fiveprime \
     --offset 0 \
     --output_format bedgraph \
     -o $merged_dir/$seqtype.full.transcriptome.nonorm.bedgraph 
    
done

for seqtype in riboseq_chx
do
    make_wiggle \
    --count_files data/$seqtype/mapped/merged/$seqtype.bam \
    --countfile_format BAM \
    --fiveprime \
    --offset 0 \
    --output_format bedgraph \
    -o $merged_dir/$seqtype.full.genome.nonorm.bedgraph 
    
done

# for seqtype in riboseq_harr riboseq_nd riboseq_harr riboseq_chx
# do
#     make_wiggle \
#     --count_files data/$seqtype/mapped/merged/$seqtype"Aligned.toTranscriptome.out.sorted.bam" \
#     --countfile_format BAM \
#     --fiveprime \
#     --offset 0 \
#     --output_format bedgraph \
#     --normalize \
#     -o $merged_dir/$seqtype.fullcov.transcriptome.bedgraph \
    
# done

conda deactivate


bamCoverage -b data/riboseq_chx/mapped/merged/riboseq_chx.bam  -o results/alignment_tracks/merged/riboseq_chx.fullcov.transcriptome.bedgraph --outFileFormat bedgraph -bs 1 -p 70 --normalizeUsing BPM 




    merged_dir=results/alignment_tracks/merged
    make_wiggle \
        --count_files data/riboseq_chx/mapped/merged/riboseq_chxAligned.toTranscriptome.out.bam \
        --countfile_format BAM \
        --fiveprime \
        --offset 0 \
        --output_format bedgraph \
        -o $merged_dir/chx.full.transcriptome.nonorm.bedgraph 