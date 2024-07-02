## --------------------------------------------------------------------------------------------------------------------
package_list <- c(
  "tidyverse", "ggpubr", "ggpp", "Biostrings", "viridis", 
  "ggforce", "grid", "GenomicFeatures", "ORFik", "cubar"
)

suppressMessages(lapply(package_list, require, character.only = TRUE))

source(file = "scripts/manuscript_prep/utility_functions.R")
source(file = "scripts/orf_discovery/orf_filtering_functions.R")


## ----make_results_dir------------------------------------------------------------------------------------------------
results_dir <- "manuscript/section_2.3/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ----load_data-------------------------------------------------------------------------------------------------------
orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")

# Genome FASTA
cgr_fasta <- FaFile(paste0("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"))

# granges for ORF-RATER identified CDS
txdb_orf <- loadTxdb("orf_identification/orf_filtered/orfrater_annotation.gtf.db",
                     chrStyle = NULL)

cds_granges <- loadRegion(txdb_orf, part = "cds")

txdb_orf <- loadTxdb("orf_identification/orf_filtered/orfrater_annotation.gtf.db",
                     chrStyle = NULL)

cds_novel_orfs <- loadRegion(txdb_orf, part = "cds")

# codon information
ctab <- get_codon_table(gcid = '1')


## ----sorf_summary----------------------------------------------------------------------------------------------------
paste0("___________________________________________________________________________________")
paste0("                   Number and % of smORFs in the uORF, ouORFs and New              ")
paste0("___________________________________________________________________________________")
orftable %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap" | `ORF type` == "New") %>%
  mutate(category = case_when(
    `Length (AAs)` < 100 ~ "microprotein",
    TRUE ~ "not_microprotein"
  )) %>%
  group_by(`ORF type`, category) %>%
  summarize(count_smorf = n()) %>%
  mutate(percentage_smorf = count_smorf / sum(count_smorf) * 100) %>%
  filter(category == "microprotein")


paste0("___________________________________________________________________________________")
paste0("                                Median smORF length.                               ")
paste0("___________________________________________________________________________________")
orftable %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap" | `ORF type` == "New") %>%
  mutate(category = case_when(
    `Length (AAs)` < 100 ~ "microprotein",
    TRUE ~ "not_microprotein"
  )) %>%
  filter(category == "microprotein") %>%
  group_by(`ORF type`) %>%
  summarise(median_length = round(mean(`Length (AAs)`)))

paste0("___________________________________________________________________________________")
paste0("              uORF and ouORF smORF category by start codon                         ")
paste0("___________________________________________________________________________________")
orftable %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap") %>%
  mutate(category = case_when(
    `Length (AAs)` < 100 ~ "microprotein",
    TRUE ~ "not_microprotein"
  )) %>%
  group_by(`ORF type`, category, `Start codon`) %>%
  filter(category == "microprotein") %>%
  summarize(count_tis = n()) %>%
  mutate(percentage_tis = round(count_tis / sum(count_tis) * 100, 1))

paste0("___________________________________________________________________________________")
paste0("                   number of smorfs in each category                               ")
paste0("___________________________________________________________________________________")
orftable %>%
  filter(`ORF type` == "New") %>%
  filter(str_detect(`ORF-RATER name`, "XR_|NR_")) %>%
  filter(`Length (AAs)` < 100) %>% 
  group_by(`ORF type`, `Start codon`) %>%
  summarize(count_tis = n()) %>%
  mutate(percentage_tis = round(count_tis / sum(count_tis) * 100, 1))


## ----fig3a-----------------------------------------------------------------------------------------------------------
XM_027392226.1 <- orftable %>%
  filter(`Transcript family` == "XM_027392226.1")

novel_CDS_start <- XM_027392226.1[XM_027392226.1$"ORF type" == "Upstream","Transcript start position"]
novel_CDS_stop <-XM_027392226.1[XM_027392226.1$"ORF type" == "Upstream","Transcript stop position"]

# account for the different encoding of the CDS on NCBI
canonical_CDS_start <- 167
canonical_CDS_stop <- 671

# load the p site offset coverage
chx_data_psite <- read_delim("sequencing/coverage_tracks/merged/fig3d/XM_027392226.1_chx_p.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  mutate(norm=(coverage/sum(coverage)),
         coverage_type = "chx_psite") %>%
  mutate(frame = "NA") %>% 
  mutate(frame = ifelse(position >= novel_CDS_start & position <= novel_CDS_stop, (position - novel_CDS_start) %% 3,
    ifelse(position >= canonical_CDS_start & position <= canonical_CDS_stop, (position - canonical_CDS_start) %% 3, NA)
  )) %>%
  mutate(scaled_coverage = scale_between_0_1(coverage))

# load the full coverage
chx_data_full <- read_delim("sequencing/coverage_tracks/merged/fig3d/XM_027392226.1_chx_f.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  mutate(norm=coverage/sum(coverage), 
         coverage_type = "chx_full")  %>%
  mutate(scaled_coverage = scale_between_0_1(coverage)) %>%
  mutate(frame = NA) 

bind_rows(chx_data_full,chx_data_psite) %>%
  ggplot(aes(x = position, y = scaled_coverage, fill = frame, alpha = as.character(coverage_type))) +
  geom_hline(yintercept = 0, linewidth=0.1) +
   annotate("rect", xmax = c(novel_CDS_start, canonical_CDS_start), 
            xmin=c(novel_CDS_stop,canonical_CDS_stop), ymin=c(1.01, -0.01), 
            ymax=c(1.1,-0.1), alpha = .7, 
            fill = c("#999933", "#CC79A7")) +
  annotate("text", x=c(87, 430), y=c(1.0595, -0.0495), label = c("34aa uORF", "Canonical Ddit3 ORF"), size=3, color="black") +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_viridis(option = "C") +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_x_continuous(limits = c(0,850),expand = c(0, 0), breaks=seq(0, 850, by = 100)) +
  # scale_y_continuous(expand = expansion(mult = c(0.1), add = c(0, 0))) +
 labs(
    x = "Transcript coordinates (nts)",
    y = "Ribo-seq coverage (Scaled BPM)",
    title = "DNA damage inducible transcript 3",
    subtitle = "XM_027392226.1"
  ) +
  theme_bw() +
  theme(validate = FALSE, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "CHX", size = 2)

ggsave(filename = paste(results_dir, "Figure 3a.png", sep = ""), 
       width = 6, height = 3, device = "png", dpi = 2000)


## ----fig3b-----------------------------------------------------------------------------------------------------------
XM_027397764.2 <- orftable %>%
  filter(`Transcript ID` == "XM_027397764.2")

novel_CDS_start <- XM_027397764.2[XM_027397764.2$"ORF type" == "Start overlap","Transcript start position"]
novel_CDS_stop <- XM_027397764.2[XM_027397764.2$"ORF type" == "Start overlap","Transcript stop position"]

canonical_CDS_start <- 202
canonical_CDS_stop <- 745


# load the p site offset coverage
chx_data_psite <- read_delim("sequencing/coverage_tracks/merged/fig3e/XM_027397764.2_chx_p.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  mutate(norm=(coverage/sum(coverage)),
         coverage_type = "chx_psite") %>%
  mutate(frame = "NA") %>% 
  mutate(frame = ifelse(position >= novel_CDS_start & position <= novel_CDS_stop, (position - novel_CDS_start) %% 3,
    ifelse(position >= canonical_CDS_start & position <= canonical_CDS_stop, (position - canonical_CDS_start) %% 3, NA)
  )) %>%
  mutate(scaled_coverage = scale_between_0_1(coverage))

# need to make a new row, because deeptools bins adjoining nt positions with an identifical value
# makes the region look aritifcally uncovered when plotting with ggplot
missing_row = data.frame(transcript = "XM_027418155.2",
                         position = 351,
                         coverage = 7.63380)

# load the full coverage
chx_data_full <- read_delim("sequencing/coverage_tracks/merged/fig3e/XM_027397764.2_chx_f.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  bind_rows(missing_row) %>%
  mutate(norm=coverage/sum(coverage), 
         coverage_type = "chx_full")  %>%
  mutate(scaled_coverage = scale_between_0_1(coverage)) %>%
  mutate(frame = NA) 


chx_data_full %>%
  filter(position > 345 & position < 400)

bind_rows(chx_data_full,chx_data_psite) %>%
  ggplot(aes(x = position, y = scaled_coverage, fill = frame, alpha = as.character(coverage_type))) +
  geom_hline(yintercept = 0, linewidth=0.1) +
   annotate("rect", xmax = c(novel_CDS_start, canonical_CDS_start), 
            xmin=c(novel_CDS_stop,canonical_CDS_stop), 
            ymin=c(1.01, -0.01), ymax=c(1.1,-0.1), 
            alpha = .7, 
            fill = c("#332288", "#CC79A7")) +
  annotate("text",x=c(210, 450), y=c(1.0595, -0.0495), label = c("68aa start overlap ORF", "Canonical Rab31 ORF"), size=2.75, color="black") +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_viridis(option = "C") +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_x_continuous(limits = c(0,820),expand = c(0, 0), breaks=seq(0, 850, by = 100)) +
  # scale_y_continuous(expand = expansion(mult = c(0.1), add = c(0, 0))) +
 labs(
    x = "Transcript coordinates (nts)",
    y = "Ribo-seq coverage (Scaled BPM)",
    title = "RAB31, member RAS oncogene family",
    subtitle = "XM_027397764.2"
  ) +
  theme_bw() +
  theme(validate = FALSE, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "CHX", size = 2)

ggsave(filename = paste(results_dir, "Figure 3b.png", sep = ""), 
       width = 6, height = 3, device = "png", dpi = 2000)


## ----fig3c-----------------------------------------------------------------------------------------------------------
XR_004770827.1 <- orftable %>%
  filter(`Transcript ID` == "XR_004770832.1" )

novel_CDS_start <- XR_004770827.1[XR_004770827.1$"ORF type" == "New","Transcript start position"]
novel_CDS_stop <- XR_004770827.1[XR_004770827.1$"ORF type" == "New","Transcript stop position"]


# load the p site offset coverage
chx_data_psite <- read_delim("sequencing/coverage_tracks/merged/fig3f/XR_004770832.1_chx_p.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  mutate(norm=(coverage/sum(coverage)),
         coverage_type = "chx_psite") %>%
  mutate(frame = "NA") %>% 
  mutate(frame = ifelse(position >= novel_CDS_start & position <= novel_CDS_stop, (position - novel_CDS_start) %% 3,NA)) %>%
  mutate(scaled_coverage = scale_between_0_1(coverage))

missing_row = data.frame(transcript = "XR_004770832.1",
                         position = 324,
                         coverage = 0.153785000)

# load the full coverage
chx_data_full <- read_delim("sequencing/coverage_tracks/merged/fig3f/XR_004770832.1_chx_f.wig", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,
                         show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) %>%
  bind_rows(missing_row) %>%
  mutate(norm=coverage/sum(coverage), 
         coverage_type = "chx_full")  %>%
  mutate(scaled_coverage = scale_between_0_1(coverage)) %>%
  mutate(frame = NA)


         
bind_rows(chx_data_full,chx_data_psite) %>%
  ggplot(aes(x = position, y = scaled_coverage, fill = frame, alpha = as.character(coverage_type))) +
  geom_hline(yintercept = 0, linewidth=0.1) +
   annotate("rect", xmin = c(novel_CDS_start, novel_CDS_start), 
            xmax=c(novel_CDS_stop,novel_CDS_stop), ymin=c(1.01, -0.01), 
            ymax=c(1.1,-0.1), alpha = .7, 
            fill = c("white", "#35968CFF")) +
  annotate("text", x=c(87, 250), y=c(1.0595, -0.0495), label = c("", "66aa novel ORF"), size=3, color="black") +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_viridis(option = "C") +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_x_continuous(limits = c(0,850),expand = c(0, 0), breaks=seq(0, 850, by = 100)) +
  # scale_y_continuous(expand = expansion(mult = c(0.1), add = c(0, 0))) +
 labs(
    x = "Transcript coordinates (nts)",
    y = "Ribo-seq coverage (Scaled BPM)",
    title = "LOC103164303 uncharacterized transcript",
    subtitle = "XR_004770827.1"
  ) +
  theme_bw() +
  theme(validate = FALSE, legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "CHX", size = 2)

ggsave(filename = paste(results_dir, "Figure 3c.png", sep = ""), 
       width = 6, height = 3, device = "png", dpi = 2000)


## ----fig3d-----------------------------------------------------------------------------------------------------------
safe_colorblind_palette <- c(
   "#CC79A7", "#88CCEE", "#888888",  "#332288",
            "#117733", "#AA4499", "#035968", "#999933")

orftable %>%
  mutate(`ORF type` = recode_factor(`ORF type`,
    `Annotated` = "Annotated",
    `Extension` = "Extension",
    `Isoform` = "Isoform",
    `Start overlap` = "Start overlap",
    `Stop overlap` = "Stop overlap",
    `Downstream` = "Downstream",
    `New` = "New",
    `Upstream` = "Upstream",
    .default = `ORF type`
  )) %>%
  ggplot(aes(y = `Length (AAs)`, x = `ORF type`, fill = `ORF type`)) +
  geom_boxplot(position = "dodge", size = 0.2, outlier.size = 0.5) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_bw() +
  scale_fill_manual(values = safe_colorblind_palette) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_zoom2(ylim = c(0, 100)) +
  labs(x = "")

ggsave(filename = paste(results_dir, "Figure 3d.png", sep = ""), 
       width = 6, height = 4, device = "png", dpi = 2000)


## ----supp_figure_4---------------------------------------------------------------------------------------------------
uorf_length_plot <- orftable %>%
  filter(`ORF type` == "Upstream" & `Length (AAs)` < 100) %>%
  ggplot(aes(x = sort(`Length (AAs)`))) +
  geom_histogram(bins = 50, fill = "#999933") +
  labs(y = "# uORFs", x = "uORF Length (aa)") +
  geom_text_npc(aes(npcx = 0.9, npcy = 0.9, label = paste0("Mean length:\n", 
                                                           round(mean(`Length (AAs)`)), " aa")), size = 3) +
  theme_bw()

ouorf_length_plot <- orftable %>%
  filter(`ORF type` == "Start overlap" & `Length (AAs)` < 100) %>%
  ggplot(aes(x = sort(`Length (AAs)`))) +
  geom_histogram(bins = 50, fill = "#332288") +
  labs(y = "# ouORFs", x = "ouORF Length (aa)") +
  theme_bw() +
  geom_text_npc(aes(npcx = 0.9, npcy = 0.9, label = paste0("Mean length:\n", 
                                                           round(mean(`Length (AAs)`)), " aa")), size = 3)

new_orf_length_plot <- orftable %>%
  filter(`ORF type` == "New" & `Length (AAs)` < 100) %>%
  ggplot(aes(x = sort(`Length (AAs)`))) +
  geom_histogram(bins = 50, fill = "#35968CFF") +
  labs(y = "# New ORFs", x = "ORF Length (aa)") +
  geom_text_npc(aes(npcx = 0.9, npcy = 0.9, label = paste0("Mean length:\n", 
                                                           round(mean(`Length (AAs)`)), " aa")), size = 3) +
  theme_bw()

ggarrange(uorf_length_plot,
  ggplot() +
    theme_void(),
  ouorf_length_plot, ggplot() +
    theme_void(),
  new_orf_length_plot,
  ncol = 5, widths = c(1, 0.05, 1, 0.05, 1)
)

ggsave(filename = paste(results_dir, "Supplementary Figure 6.png", sep = ""),
       dpi = 2000, device = "png", width = 9, height = 3)


## ----aa_usage--------------------------------------------------------------------------------------------------------
# upstream microproteins
upstream_sorf <- orftable %>%
  filter(`Length (AAs)` < 100) %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap")

fasta_out_file="manuscript/section_2.3/upstream_sorfs.fasta"
orf_to_fasta(upstream_sorf, cds_granges, fasta_out_file, cgr_fasta)

# lncRNA encoded
lncrna_sorf <- orftable %>%
  filter(`Length (AAs)` < 100) %>%
  filter(`ORF type` == "New") %>%
  filter(str_detect(`ORF-RATER name`, "XR|NR"))

fasta_out_file="manuscript/section_2.3/lncrna_sorfs.fasta"
orf_to_fasta(lncrna_sorf, cds_granges, fasta_out_file, cgr_fasta)

# calculate expected amino acid frequency

# frequency of nucleotides in the PICR genome
key_mapping <- c("A" = 0.2671548, "T" = 0.2217720, "C" = 0.2538705, "G" = 0.2572025)

# translate
frequency_df <- data.frame(matrix(NA, nrow = 64, ncol = 3))
colnames(frequency_df) <- c("Amino Acid","Codon","Expected Frequency")

for (i in 1:64) {
  frequency_df[i,1] <- GENETIC_CODE[[names(GENETIC_CODE)[i]]]
  frequency_df[i,2] <- names(GENETIC_CODE)[i]
  codon_string <- unlist(str_split(names(GENETIC_CODE)[i],""))
  frequency_df[i,3] <- key_mapping[codon_string[1]] * key_mapping[codon_string[2]] *  key_mapping[codon_string[3]]
}

frequency_df <- frequency_df %>%
  filter(`Amino Acid` != "*") %>% # no stop
  group_by(`Amino Acid`) %>%
  summarize(`Expected Frequency` = sum(`Expected Frequency`))

# NCBI annotated proteome
fasta <- readAAStringSet("reference_proteome/GCF_003668045.3_CriGri-PICRH-1.0_protein.faa")
annotated_protein_coding <-colSums(alphabetFrequency(fasta)[, AA_STANDARD]) / (sum(colSums(alphabetFrequency(fasta)[, AA_STANDARD])))
annotated_protein_coding <- data.frame(annotated_protein_coding) %>%
as_tibble(rownames = "Amino Acid")

fasta <- readAAStringSet("manuscript/section_2.3/upstream_sorfs.fasta")
short_upstream_ORF <- colSums(alphabetFrequency(fasta)[, AA_STANDARD]) / (sum(colSums(alphabetFrequency(fasta)[, AA_STANDARD])))
short_upstream_ORF<- data.frame(short_upstream_ORF) %>%
as_tibble(rownames = "Amino Acid")

fasta <- readAAStringSet("manuscript/section_2.3/lncrna_sorfs.fasta")
short_ncRNA_ORF <- colSums(alphabetFrequency(fasta)[, AA_STANDARD]) / (sum(colSums(alphabetFrequency(fasta)[, AA_STANDARD])))
short_ncRNA_ORF<- data.frame(short_ncRNA_ORF) %>%
as_tibble(rownames = "Amino Acid")

frequency_df <- frequency_df %>%
  left_join(annotated_protein_coding, by = "Amino Acid")  %>%
  left_join(short_upstream_ORF, by = "Amino Acid")  %>%
  left_join(short_ncRNA_ORF, by = "Amino Acid") 

frequency_df <- frequency_df %>%
  dplyr::rename("Expected" = `Expected Frequency`,
                "Annotated\n(>= 100aa)" = annotated_protein_coding,
                "New\n(< 100aa)" = short_ncRNA_ORF,
                "Upstream &\nStart overlap\n(< 100 aa)" = short_upstream_ORF)


## ----supp_figure_5---------------------------------------------------------------------------------------------------
# display all AAs
frequency_df %>%
  pivot_longer(names_to = "type", values_to = "Freq",  `Expected`:`New\n(< 100aa)`) %>%
  # mutate(type = c("Expected","Annotated >= 100aa", "uORFs < 100aa", "New < 100aa")) %>%
  # pivot_longer(names_to = "AA", values_to = "Freq", cols = -contains("type")) %>%
  mutate(aa_name = AMINO_ACID_CODE[`Amino Acid`]) %>%
  ggplot(aes(x = aa_name, y = Freq, fill = type)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  theme_bw() +
  labs(x = "", y = "Frequency", fill = "") +
  scale_fill_manual(values = c("#CC79A7","#79A7CC", "#35968CFF", "#999933")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult=c(0, 0.05))) +
  theme(legend.position = "top")

ggsave(filename = paste(results_dir, "Supplementary Figure 7.png", sep = ""), 
       width = 7, height = 3.5, device = "png", dpi = 2000)


## ----fig3e-----------------------------------------------------------------------------------------------------------
frequency_df %>%
  pivot_longer(names_to = "type", values_to = "Freq",  `Expected`:`New\n(< 100aa)`) %>%
  # mutate(type = c("Expected","Annotated >= 100aa", "uORFs < 100aa", "New < 100aa")) %>%
  # pivot_longer(names_to = "AA", values_to = "Freq", cols = -contains("type")) %>%
  mutate(aa_name = AMINO_ACID_CODE[`Amino Acid`]) %>%
    filter(aa_name != "His" & aa_name != "Leu" & aa_name != "Ser" & aa_name != "Val" & aa_name != "Thr" & aa_name != "Phe") %>%
  ggplot(aes(x = aa_name, y = Freq, fill = type)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  theme_bw() +
  labs(x = "", y = "Frequency", fill = "") +
  scale_fill_manual(values = c("#CC79A7","#79A7CC", "#35968CFF", "#999933")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult=c(0, 0.05))) +
  theme(legend.position = "top", legend.text = element_text(size=7))

ggsave(filename = paste(results_dir, "Figure 3e.png", sep = ""), 
       width = 6, height = 3.5, device = "png", dpi = 2000)


## ----fig3f-----------------------------------------------------------------------------------------------------------
#calculate RSCU for highly translated

canonical <- orftable %>%
  filter(`Length (AAs)` >= 100) %>%
  filter(`ORF type` == "Annotated") 

# Import the RiboSeq annotated CDS
canonical_orf_granges <- cds_granges[names(cds_granges) %in% canonical$`ORF-RATER name`]

canonical_orf_cds <- extractTranscriptSeqs(cgr_fasta, 
                                               canonical_orf_granges)

# highly expressed genes from the RNA-seq data are used as a reference set
high_translation_set <- orftable %>%
  filter(`ORF type` == "Annotated") %>%
  filter(`Length (AAs)` >= 100) %>%
 slice_max(fpkmRFP,n = 500)

# determine the couunts 
canonical_cf <- count_codons(canonical_orf_cds)

high_translation_set <- high_translation_set[high_translation_set$`ORF-RATER name` %in% rownames(canonical_cf), ]

rscu_heg <- est_rscu(canonical_cf[high_translation_set$`ORF-RATER name`,], 
                     codon_table = ctab)

# ncbi annotated ORF codon frequency
annotated_cds<-readDNAStringSet("reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_cds_from_genomic.fna")
annotated_cf <- count_codons(annotated_cds)

# upstream ORF codon frequency
upstream_sorf <- orftable %>%
  filter(`Length (AAs)` < 100) %>%
  filter(`ORF type` == "Upstream" | `ORF type` == "Start overlap")

upstream_orf_granges <- cds_granges[names(cds_granges) %in% upstream_sorf$`ORF-RATER name`]
upstream_orf_cds <- extractTranscriptSeqs(cgr_fasta, upstream_orf_granges)
upstream_orf_cf <- count_codons(upstream_orf_cds)

# new ORF codon frequency
lncrna_sorf <- orftable %>%
  filter(`Length (AAs)` < 100) %>%
  filter(`ORF type` == "New") %>%
  filter(str_detect(`ORF-RATER name`, "XR|NR")) # only use ORFs on ncRNA

new_orf_granges <- cds_granges[names(cds_granges) %in% lncrna_sorf$`ORF-RATER name`]
new_orf_cds <- extractTranscriptSeqs(cgr_fasta, new_orf_granges)
new_orf_cf <- count_codons(new_orf_cds)

# codon adaptation index
# determine the CAI for the three of sets
combined_cai <- bind_rows(
  data.frame(cai=get_cai(annotated_cf, rscu = rscu_heg)) %>%
            mutate(class="Annotated\n(>= 100 aa)"),
            data.frame(cai=get_cai(upstream_orf_cf, rscu = rscu_heg))%>%
    mutate(class="Upstream &\nStart overlap\n(< 100 aa)"),
  data.frame(cai=get_cai(new_orf_cf, rscu = rscu_heg)) %>%
  mutate(class="New\n(< 100 aa)")
  )

# Kolmogorov–Smirnov
ks_test_can_uORF <- ks.test(combined_cai$cai[combined_cai$class == "Annotated\n(>= 100 aa)"], combined_cai$cai[combined_cai$class == "Upstream &\nStart overlap\n(< 100 aa)"])
ks_test_can_ncORF <- ks.test(combined_cai$cai[combined_cai$class == "Annotated\n(>= 100 aa)"], combined_cai$cai[combined_cai$class == "New\n(< 100 aa)"])
ks_test_ncORF_uORF <- ks.test(combined_cai$cai[combined_cai$class == "New\n(< 100 aa)"], combined_cai$cai[combined_cai$class == "Upstream &\nStart overlap\n(< 100 aa)"])

stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
    "Annotated\n(>= 100 aa)",     "Upstream &\nStart overlap\n(< 100 aa)", "< 2.2e-16",
    "Annotated\n(>= 100 aa)",     "New\n(< 100 aa)", "< 2.2e-16",
      "New\n(< 100 aa)",     "Upstream &\nStart overlap\n(< 100 aa)", "0.42"
  )

ggboxplot(combined_cai, x = "class", y = "cai", fill="class") +
    stat_pvalue_manual(
    stat.test,
    size=3,
    y.position = 1.15, step.increase = 0.1,
    label = "p.adj"
    ) +
  scale_fill_manual(values=c("#CC79A7", "#999933", "#035968")) +
  labs(x = "", y="Codon Adaptation Index (CAI)", fill = "ORF type") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = paste(results_dir, "Figure 3f.png", sep = ""), 
       dpi = 2000, device = "png", height =3.5, width = 4)


## --------------------------------------------------------------------------------------------------------------------
print("section 2.3 results complete")

