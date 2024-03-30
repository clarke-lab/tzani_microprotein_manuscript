#!/usr/bin/env Rscript --vanilla

# Author: Colin Clarke
# Date:   2024-02-06
# Purpose: This script returns high confidence novel CGR ORFs
# from the ORFRATER output

# load libraries
suppressMessages(library(ORFik))
suppressMessages(library(tidyverse))

# load custom functions
source(paste0(getwd(),"/scripts/orf_discovery/orf_filtering_functions.R"))

# make a folder for outputs
out_dir=paste0(getwd(), "/orf_identification/orf_filtered/")

if (!dir.exists(out_dir)) {
    system(paste0("mkdir ", getwd(), "/orf_identification/orf_filtered"), 
    intern = T, ignore.stderr = T, ignore.stdout = T)
}

############################################################
# 1. Genome Annotation Data
############################################################

ref_genome_dir <- paste0(getwd(),"/reference_genome/")

# make a txdb for the NCBI reference annotation
if (!file.exists(paste0(out_dir,"PICRH_ORFik.gtf.db"))) {

    # GTF file compatible with ORFik
    system(paste0(
    "sed 's/XM_/XM./g; s/XR_/XR./g; s/NR_/NR./g; s/NM_/NM./g' ",
    ref_genome_dir, "GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf | ",
    "grep -v unassigned_transcript ",
    "> ", out_dir, "PICRH_ORFik.gtf"), 
    intern = TRUE)

    gtf <- paste0(out_dir, "PICRH_ORFik.gtf")
    fasta <- paste0(ref_genome_dir,"GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna")
    organism <- "Cricetulus griseus"

    print("Constructing NCBI reference TXDB")
    makeTxdbFromGenome(gtf,
                   fasta,
                   organism,
                   optimize = TRUE,
                   pseudo_5UTRS_if_needed = 100)

}

# load reference txdb
txdb <- loadTxdb(paste0(out_dir, "PICRH_ORFik.gtf.db"),
                 chrStyle = NULL)


############################################################
# 2. import the Ribo-seq alignments
############################################################

# p-site offsets for each read length [determined via plastid]
shifts <- data.frame(cbind(fraction = c(28, 29, 30, 31), 
                           offsets_start = c(-12, -12, -12, -12)))

# 2 genome aligned

data_dir <- paste0(getwd(),"/sequencing/")

if (!file.exists(paste0(out_dir,"genome_footprints.rds"))) {

    # p-site offsets for each read length [determined via plastid]
    shifts <- data.frame(cbind(fraction = c(28, 29, 30, 31), 
    offsets_start = c(-12, -12, -12, -12)))

    bam_file_chx_genome <- paste0(data_dir,"riboseq_chx/mapped/merged/riboseq_chx.bam")
    genome_footprints_chx <- readBam(bam_file_chx_genome)
    shiftedFootprints_chx <- shiftFootprints(genome_footprints_chx, shifts)

    bam_file_harr_genome <- paste0(data_dir,"riboseq_harr/mapped/merged/riboseq_harr.bam")
    genome_footprints_harr <- readBam(bam_file_harr_genome)
    shiftedFootprints_harr <- shiftFootprints(genome_footprints_harr, shifts)

    bam_file_nd_genome <- paste0(data_dir,"riboseq_nd/mapped/merged/riboseq_nd.bam")
    genome_footprints_nd <- readBam(bam_file_nd_genome)
    shiftedFootprints_nd <- shiftFootprints(genome_footprints_nd, shifts)

    bam_file_rna <- paste0(data_dir,"rnaseq_se/mapped/merged/rnaseq_se.bam")
    genome_rna <- readBam(bam_file_rna)

    genome_footprints <- list(
        shiftedFootprints_chx = shiftedFootprints_chx,
        shiftedFootprints_harr = shiftedFootprints_harr,
        shiftedFootprints_nd = shiftedFootprints_nd,
        genome_rna = genome_rna
    )

    saveRDS(genome_footprints, 
    file = paste0(out_dir,"genome_footprints.rds"))

} else {
    
    genome_footprints <- readRDS(paste0(out_dir,"genome_footprints.rds"))

    shiftedFootprints_chx = genome_footprints$shiftedFootprints_chx
    shiftedFootprints_harr = genome_footprints$shiftedFootprints_harr
    shiftedFootprints_nd = genome_footprints$shiftedFootprints_nd
    genome_rna = genome_footprints$genome_rna
}

############################################################
# 3. Initial filter for ORF-RATER outputs
############################################################

orfrater_dir=paste0(getwd(), "/orf_identification/orfrater/")

# 3a. Filter by ORF-Rate score and protein length
orftable <- read.csv(file = paste0(orfrater_dir,"rate_regression.csv"), header = T) %>%
  filter(orfrating >= 0.5 & AAlen >= 5)

paste0(dim(orftable)[1], " ORFs with an orfrating >= 0.5 and longer than 5aa")

# 3b. Remove unwanted ORF-types
orftable <- orftable %>%
  filter(!orftype %in% c("Ciso", "Giso", "NCiso", "Niso", 
  "new_iso", "internal", "LOOF", "truncation"))

# 3c. import NCBI annotation table 
suppressWarnings(
  suppressMessages(
    CriGri_PICRH_1_0_annotation <- read_delim(
      paste0(ref_genome_dir,"GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt"),
      "\t",escape_double = FALSE, trim_ws = TRUE) %>%
      filter(str_detect(product_accession, "XM_|NM_|NR_|XR_")) %>%
      mutate(tid = product_accession)
  )
)

# 3d. add NCBI transcript annotation, select and rename columns
orftable_annotated <- orftable %>%
  left_join(CriGri_PICRH_1_0_annotation, by = "tid") %>%
  dplyr::select(
    `ORF-RATER name` = orfname,
    `ORF type` = orftype,
    `ORF-RATER score` = orfrating,
    `Associated Gene symbol` = symbol,
    `Associated Gene Name` = name,
    `Transcript family` = tfam,
    `Transcript ID` = tid,
    `Start codon` = codon,
    `Length (AAs)` = AAlen,
    `Start Annotated?` = annot_start,
    `Stop Annotated` = annot_stop,
    `Transcript start position` = tcoord,
    `Transcript stop position` = tstop,
    `Chromosome` = chrom,
    `Strand` = strand.x,
    `Genomic start position` = gcoord,
    `Genomic stop position` = gstop
  )

# 3e. Rename ORF types
orftable_annotated <- orftable_annotated %>%
  mutate(`ORF type` = case_when(
    `ORF type` == "Xiso" ~ "Isoform",
    `ORF type` == "Siso" ~ "Isoform",
    `ORF type` == "annotated" ~ "Annotated",
    `ORF type` == "new" ~ "New",
    `ORF type` == "upstream" ~ "Upstream",
    `ORF type` == "downstream" ~ "Downstream",
    `ORF type` == "start_overlap" ~ "Start overlap",
    `ORF type` == "stop_overlap" ~ "Stop overlap",
    `ORF type` == "extension" ~ "Extension",
  ))


paste0(dim(orftable_annotated)[1], " remain following removal of truncated, internal, LOOF and low confidence isoforms")
print("ORFs type by start codon")
print(table(orftable_annotated$`ORF type`, orftable_annotated$`Start codon`))

saveRDS(orftable_annotated, 
    file = paste0(out_dir,"orfs_remaining_post_first_filter.rds"))

############################################################
# 4. Rharr-Rnd TIS enrichment filter
############################################################

# 4a. import Ribo-seq counts aligned against transcriptome
if (!file.exists(paste0(out_dir,"transcriptome_footprints.rds"))) {

    bam_file_harr <- paste0(data_dir,
    "riboseq_harr/mapped/merged/riboseq_harrAligned.toTranscriptome.out.bam")
    footprints_harr_transcriptome <- readBam(bam_file_harr)

    bam_file_nd <- paste0(data_dir,
    "riboseq_nd/mapped/merged/riboseq_ndAligned.toTranscriptome.out.bam")
    footprints_nd_transcriptome <- readBam(bam_file_nd)

    transcriptome_footprints <- list(
        footprints_harr_transcriptome = footprints_harr_transcriptome,
        footprints_nd_transcriptome = footprints_nd_transcriptome
    )

    saveRDS(transcriptome_footprints, 
    file = paste0(out_dir,"transcriptome_footprints.rds"))

} else {
    
    transcriptome_footprints <- readRDS(paste0(out_dir,"transcriptome_footprints.rds"))

    footprints_harr_transcriptome = transcriptome_footprints$footprints_harr_transcriptome
    footprints_nd_transcriptome = transcriptome_footprints$footprints_nd_transcriptome 

}

# 4b.Reads per transcript in each transcriptome aligned BAM
harr_transcript_counts <- count_transcript_alignments(footprints_harr_transcriptome)
nd_transcript_counts <- count_transcript_alignments(footprints_nd_transcriptome)

# 4c. merge the transcript count matrices 
transcript_counts <- harr_transcript_counts %>%
  dplyr::rename(nharr = counts) %>%
  left_join(nd_transcript_counts, by = "transcript") %>%
  dplyr::rename(nnd = counts) %>%
  dplyr::rename("Transcript ID" = transcript)

# 4c. count the reads overlapping the TIS in Harr and ND Riboseq
# Make a Grange object for ORF TISs
orf_start_codons <- GRanges(
  seqnames = orftable_annotated$`Transcript ID`,
  ranges = IRanges(
    start = orftable_annotated$`Transcript start position`-2,
    width = 5,
    strand = "*"),
  names = orftable_annotated$`ORF-RATER name`
)

# count reads overlapping the ORF start codons 
orf_tis_counts_harr <- countOverlaps(orf_start_codons, 
footprints_harr_transcriptome)

orf_tis_counts_nd <- countOverlaps(orf_start_codons, 
footprints_nd_transcriptome)

# merge TIS overlap counts, counts on full transcript, 
# calculate Rharr-Rnd
# retain ORFs above minium threshold
orftable_tis_filtered <- bind_cols(
  orftable_annotated,
  data.frame(
    xharr = orf_tis_counts_harr,
    xnd = orf_tis_counts_nd)) %>%
    filter(xharr >= 10) %>%
  left_join(transcript_counts, by = "Transcript ID") %>%
  mutate(
    Rharr = (xharr / nharr) * 10,
    Rnd = (xnd / nnd) * 10,
    Rharr_min_Rnd = Rharr - Rnd
  ) %>%
  mutate(tis_enriched = case_when(
    Rharr_min_Rnd >= 0.01 ~ "TRUE",
    TRUE ~ "FALSE"
  )) %>%
  filter(tis_enriched == "TRUE") %>%
  dplyr::select(-c(tis_enriched, Rharr, Rnd))

paste0(dim(orftable_tis_filtered)[1], " remain following TIS enirchment filtering")
print("ORFs types by start codon")
table(orftable_tis_filtered$`ORF type`, orftable_tis_filtered$`Start codon`)

saveRDS(orftable_tis_filtered, 
    file = paste0(out_dir,"orfs_remaining_post_tis_filter.rds"))

############################################################
# 5. Prepare for ORFik and FLOSS analysis
############################################################

# 5a. create an ORFIK compatiable ID
orftable_tis_filtered <- orftable_tis_filtered %>%
  arrange(`Transcript ID`) %>%
  mutate(orfik_id = ave(`Transcript ID`,
                        `Transcript ID`,
                        FUN = function(x) {
                          ifelse(duplicated(x),
                                 paste0(x, "_", seq_along(x)),
                                 paste0(x, "_1")
                          )
                        }
  )) %>%
  mutate(orfik_id = sub("_", ".", orfik_id))

# 5b. Make ORF-RATER txdb

orfrater_dir=paste0(getwd(), "/orf_identification/orfrater/")

system(paste0("cp ", orfrater_dir, "orfrater_annotation.gtf ", 
out_dir))

if (!file.exists(paste0(out_dir,"orfrater_annotation.gtf.db"))) {

    gtf_orf <- paste0(out_dir, "orfrater_annotation.gtf")
    fasta <- paste0(ref_genome_dir,"GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna")
    organism <- "Cricetulus griseus"

    print("Constructing ORF-RATER txdb")
    makeTxdbFromGenome(gtf_orf,
                    fasta,
                    organism,
                    optimize = TRUE,
                    pseudo_5UTRS_if_needed = 100)
    
}

txdb_orf <- loadTxdb(paste0(out_dir,"orfrater_annotation.gtf.db"),
                     chrStyle = NULL)

# extract ORFRATER CDS
cds_novel_orfs <- loadRegion(txdb_orf, part = "cds")

# reorder and add orfik_id to ORFRATER CDS Granges
order_indices <- match(orftable_tis_filtered$`ORF-RATER name`, names(cds_novel_orfs))
sum(orftable_tis_filtered$`ORF-RATER name` == names(cds_novel_orfs[order_indices])) / length(orftable_tis_filtered$`ORF-RATER name`)
cds_novel_orfs <- cds_novel_orfs[order_indices]
names(cds_novel_orfs) <- orftable_tis_filtered$orfik_id

# 5c. Define ORF group CDS Granges for analysis
# Annotated
reference_orfs <- orftable_tis_filtered %>%
  filter(`ORF type` == "Annotated")

reference_cds_orfs <- cds_novel_orfs[names(cds_novel_orfs) %in% reference_orfs$orfik_id]

# Upstream
upstream_novel_orfs <- orftable_tis_filtered %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap")

upstream_novel_orfs_cds <- cds_novel_orfs[names(cds_novel_orfs) %in% upstream_novel_orfs$orfik_id]

# Other
other_novel_orfs <- orftable_tis_filtered %>%
  filter(`ORF type` != "Annotated" & `ORF type` != "Upstream" & `ORF type` != "Start overlap")

other_novel_orfs_cds <- cds_novel_orfs[names(cds_novel_orfs) %in% other_novel_orfs$orfik_id]

############################################################
# 6. ORFIK analysis using CHX RiboSeq
############################################################

# load the reference genome fasta
cgr_fasta <- FaFile(paste0(ref_genome_dir,"GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"))

# orfrate identified NCBI annotated cds
orfik_chx_reference_orfs <- computeFeatures(
  grl = reference_cds_orfs,
  RFP = shiftedFootprints_chx,
  RNA = genome_rna,
  Gtf = txdb,
  faFile = cgr_fasta,
  sequenceFeatures = T,
  riboStart = 28,
  riboStop = 31,
  uorfFeatures = F
)

reference_orfs <- bind_cols(reference_orfs, orfik_chx_reference_orfs)

# upstream cds
orfik_chx_upstream_orfs <- computeFeatures(
  grl = upstream_novel_orfs_cds,
  RFP = shiftedFootprints_chx,
  RNA = genome_rna,
  Gtf = txdb,
  faFile = cgr_fasta,
  sequenceFeatures = T,
  riboStart = 28,
  riboStop = 31,
  uorfFeatures = T
)

upstream_novel_orfs <- bind_cols(upstream_novel_orfs, orfik_chx_upstream_orfs)

# other cds
orfik_chx_other_orfs <- computeFeatures(
  grl = other_novel_orfs_cds,
  RFP = shiftedFootprints_chx,
  RNA = genome_rna,
  Gtf = txdb,
  faFile = cgr_fasta,
  sequenceFeatures = T,
  riboStart = 28,
  riboStop = 31,
  uorfFeatures = F
)

other_novel_orfs <- bind_cols(other_novel_orfs, orfik_chx_other_orfs)

orftable_tis_filtered <- bind_rows(
  reference_orfs,
  upstream_novel_orfs,
  other_novel_orfs
)


############################################################
# 5. FLOSS classification and ORF filtering
############################################################

# all annotated NCBI CDS
cds <- loadRegion(txdb, part = "cds")
ncbi_annotated_orfs <- computeFeatures(
  grl = cds,
  RFP = shiftedFootprints_chx,
  RNA = genome_rna,
  Gtf = txdb,
  faFile = cgr_fasta,
  sequenceFeatures = T,
  riboStart = 28,
  riboStop = 31,
  uorfFeatures = F
)

# determine the cutoffs for FLOSS based on ref CDS with > 10 counts
classifier_ref <- flossCutoffs(
  ncbi_annotated_orfs$countRFP[ncbi_annotated_orfs$countRFP > 10],
  ncbi_annotated_orfs$floss[ncbi_annotated_orfs$countRFP > 10]
)

# apply cutoff to classify
orftable_tis_filtered$floss_classification <- as.character(
  classifier_ref$classify(orftable_tis_filtered$countRFP, 
                                  orftable_tis_filtered$floss)
)

# add classifications to orftable
orftable_tis_filtered <- orftable_tis_filtered %>%
  mutate(floss_classification = case_when(
    floss_classification == "Good" ~ "Good",
    floss_classification == "Extreme" ~ "Extreme",
    is.na(floss_classification) ~ "ok"
  ))

# remove ORFs with extreme floss classification
orftable_final <- orftable_tis_filtered %>%
  filter(floss_classification != "Extreme")

# save the annotation
saveRDS(orftable_final, 
file = paste0(out_dir,"final_orfs.rds"))

paste0(dim(orftable_final)[1], " remain following FLOSS filter")
print("ORFs type by start codon")
print(table(orftable_final$`ORF type`, orftable_final$`Start codon`))

############################################################
# 7. Make the protein fasta for proteomics 
############################################################

microprotein_annotations <- orftable_final %>%
  filter(`ORF type` != "Annotated") %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap" | `ORF type` == "New") %>%
  filter(`Length (AAs)` < 100)

print("microprotein type by start codon")
table(microprotein_annotations$`ORF type`, microprotein_annotations$`Start codon`)

cds_novel_orfs <- loadRegion(txdb_orf, part = "cds")

# create the fasta file
fasta_out_file=paste0(out_dir,"microproteins.fasta")
orf_to_fasta(microprotein_annotations, cds_novel_orfs, fasta_out_file, cgr_fasta)

############################################################
# 8. output list to be used for DE analysis
############################################################

selected_new_for_ts_de <- orftable_final %>%
  filter(`ORF type` == "New") %>%
  filter(str_detect(`ORF-RATER name`, "NR|XR")) %>%
  group_by(`Transcript ID`) %>%
  top_n(`Length (AAs)`, n = 1)

write(selected_new_for_ts_de$`ORF-RATER name`, 
file = paste0(out_dir,"lncrna_novel_orfs.txt"))
