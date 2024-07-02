## ----setup-----------------------------------------------------------------------------------------
package_list <- c("tidyverse", "writexl", "ggpubr", "viridis", "patchwork")

lapply(package_list, require, character.only = TRUE)

options(dplyr.summarise.inform = FALSE)

# source("scripts/manuscript/utility_functions.R")


## ----make_results_dir------------------------------------------------------------------------------
results_dir <- "manuscript/section_2.1/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ---- cell_density_info----------------------------------------------------------------------------
# import the cell density data for both the CHX and HARR riboseq data
cell_density_data <- read_delim("data/cell_density/ts_cell_density_ngs.txt", "\t",
  escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# calculate the difference in cell density at 72hrs post seeding for the TS and NTS cell lines
cd_summary <- cell_density_data %>%
  group_by(Experiment, condition) %>%
  summarise(mean_density = mean(cell_density)) 

paste0("___________________________________________________________________________________ ")
paste0("                         Initiation experiment.                                    ")
paste0("___________________________________________________________________________________")

initiation_ts_delta <- cd_summary %>%
  filter(Experiment == "Initiation") %>%
  summarise(`CD diff` = 100 - round((mean_density[condition == "TS"] / mean_density[condition == "NTS"]) * 100))

paste0("___________________________________________________________________________________ ")
paste0("                         Elongation experiment.                                    ")
paste0("___________________________________________________________________________________")
paste0("Approx. ", initiation_ts_delta$`CD diff`, "% cell density reduction in CD in TS samples")

elongation_ts_delta <- cd_summary %>%
  filter(Experiment == "Elongation") %>%
  summarise(`CD diff` = 100 - round((mean_density[condition == "TS"] / mean_density[condition == "NTS"]) * 100))

paste0("Approx. ", elongation_ts_delta$`CD diff`, "% cell density reduction in CD in TS samples")


## ---- supp_data_1----------------------------------------------------------------------------------
supp_data_1 <- cell_density_data %>%
  pivot_wider(names_from= c(condition),values_from=cell_density) %>%
  rename(`NTS cell density (cells/ml)`=NTS, `TS cell density (cells/ml)`=TS)

fn <- paste(results_dir, "Supplementary Data 1a.xlsx", sep = "")

suppressMessages(if (file.exists(fn)) {
  file.remove(fn)})

write_xlsx(list(`A`=supp_data_1),
path = fn, format_headers = TRUE)


## ----supp_figure_1---------------------------------------------------------------------------------
ggboxplot(cell_density_data, x = "condition", y = "cell_density", 
                               color = "black", add = "jitter", fill = "condition") +
  facet_wrap(~Experiment, ncol = 2) +
  theme_bw() +
  stat_compare_means(method = "t.test", label.y = 2100000, label.x = 1.65) +
  labs(x = "", y = "Cell Density (cells/ml)", fill = "") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9")) +
  theme(strip.background = element_rect(fill = "white"))

ggsave(filename = paste(results_dir, "Supplementary Figure 1.png", sep = ""), 
       width = 6, height = 3, device = "png", dpi = 2000)


## ----import_read_counts----------------------------------------------------------------------------
combined_counts <- read_delim("sequencing/stats/read_counts/combined_counts.txt", 
    delim = "\t", escape_double = FALSE, show_col_types = FALSE,trim_ws = TRUE) %>%
  filter(file !="file") %>%
  mutate(raw_read_number=as.numeric(raw_read_number))

sample_file_list <- fs::dir_ls("data/seq_files/", regexp = "\\.txt$")

sample_ids <- map_dfr(sample_file_list, read_tsv, col_names=F, show_col_types = FALSE) %>%
  mutate(sample = gsub(".fastq.gz","",X1))

combined_counts <- combined_counts %>%
  mutate(seq_type=case_when(
    str_detect(file, "riboseq_chx") ~ "Cycloheximide",
    str_detect(file, "riboseq_harr") ~ "Harringtonine",
    str_detect(file, "riboseq_nd") ~ "No drug", 
    str_detect(file, "rnaseq_se") ~ "RNAseq",
  )) 

combined_counts <- combined_counts %>%
  mutate(stage=case_when(
    str_detect(file, "raw_data") ~ "raw",
    str_detect(file, "trimmed") ~ "trimmed",
    str_detect(file, "rRNA") ~ "rRNA_retained", 
    str_detect(file, "snoRNA") ~ "snoRNA_retained",
    str_detect(file, "tRNA") ~ "tRNA_retained",
    str_detect(file, "complete") ~ "passed",
  )) %>%
  mutate(sample=gsub("_unaligned.fq", "",gsub("\\.fastq(?:\\.gz)?$", "", gsub(".*/", "", file))))

combined_counts <- combined_counts %>%
left_join(sample_ids,by = "sample") %>%
  mutate(sample=case_when(
    str_detect(sample,"SRR")  ~ X2,
    !str_detect(sample,"SRR") ~ sample
  )) %>%
  select(-X1, -X2)


## ----transform_read_counts-------------------------------------------------------------------------
riboseq_counts <- combined_counts %>%
  filter(seq_type != "RNAseq") %>%
  pivot_wider(id_cols = c(sample, seq_type), 
                      names_from = stage, 
                      values_from = raw_read_number) %>%
  mutate(
  Raw=raw,
  Cutadapt = raw - trimmed,
  rRNA = trimmed - rRNA_retained,
  snoRNA = rRNA_retained - snoRNA_retained,
  tRNA = snoRNA_retained - tRNA_retained,
  `<28nt or >32nt` = tRNA_retained - passed,
  Retained = passed
  ) %>%
  dplyr::select(sample, seq_type, c("Raw","Cutadapt", "rRNA", "snoRNA", "tRNA", 
                                  "<28nt or >32nt", "Retained")) %>%
  mutate(sample = toupper(sample))


## ----riboseq_read_info-----------------------------------------------------------------------------
total_read_counts<- riboseq_counts %>%
  rowwise() %>%
  mutate(raw=sum(c_across(Cutadapt:Retained))) %>% 
  group_by(seq_type) %>%
  summarise(total_sequenced=round(mean(raw)/1e6))

paste0("___________________________________________________________________________________")
paste0("                         Ribo-seq reads per sample                                          ")
paste0("___________________________________________________________________________________")
paste0("Average of ", total_read_counts[1,2], " million raw reads sequenced for each ", total_read_counts[1,1], " Ribo-seq sample.")
paste0("Average of ", total_read_counts[2,2], " million raw reads sequenced for each ", total_read_counts[2,1], " Ribo-seq sample.")
paste0("Average of ", total_read_counts[3,2], " million raw reads sequenced for each ", total_read_counts[3,1], " Ribo-seq sample.")

ncrna_removed <-riboseq_counts %>%
  rowwise() %>%
  mutate(ncRNA_removed = sum(rRNA,snoRNA,tRNA)) %>%
    rowwise() %>%
  mutate(percent_nc_removed = ncRNA_removed/Raw) %>%
  group_by(seq_type) %>%
  summarise(round(mean(percent_nc_removed)*100))

paste0("___________________________________________________________________________________")
paste0("                      ncRNA mapped Ribo-seq reads removed                          ")
paste0("___________________________________________________________________________________")
paste0(ncrna_removed[1,2], " % of ", ncrna_removed[1,1], " Ribo-seq reads originating from non-coding RNA were removed")
paste0(ncrna_removed[2,2], " % of ", ncrna_removed[2,1], " Ribo-seq reads originating from non-coding RNA were removed")
paste0(ncrna_removed[3,2], " % of ", ncrna_removed[3,1], " Ribo-seq reads originating from non-coding RNA were removed")

retained_total <- riboseq_counts %>%
  group_by(seq_type) %>%
  summarise(floor(sum(Retained)/1e6))

paste0("___________________________________________________________________________________")
paste0("                    Ribo-seq reads reatained post-processing                       ")
paste0("___________________________________________________________________________________")
paste0("Total of ", retained_total[1,2], " million preprocessed ", retained_total[1,1], " RPFs for ORF-RATER")
paste0("Total of ", retained_total[2,2], " million preprocessed ", retained_total[2,1], " RPFs for ORF-RATER")
paste0("Total of ", retained_total[3,2], " million preprocessed ", retained_total[3,1], " RPFs for ORF-RATER")


## ----supp_figure_2---------------------------------------------------------------------------------
stack_order <- c(
  "Cutadapt", "rRNA", "snoRNA", "tRNA", "<28nt or >32nt", "Retained"
)

riboseq_counts %>%
  pivot_longer(cols = c("Cutadapt", "rRNA", "snoRNA", "tRNA", "<28nt or >32nt", "Retained"), 
               values_to = "count", names_to = "stage") %>%
  ggplot(aes(x = sample, y = count, fill = factor(stage, levels = stack_order))) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  # scale_fill_brewer(palette = "Paired") +
  labs(x = "", y = "Proportion of reads", fill = "") +
  facet_wrap(~seq_type, nrow = 3) +
  scale_fill_viridis(discrete = T, option="H", direction = 1) +
  theme(strip.background = element_rect(fill = "white"))

ggsave(filename = paste(results_dir, "Supplementary Figure 2.png", sep = ""), 
       width = 7, height = 8, device = "png", dpi = 2000)


## ----rnaseq_seq_info-------------------------------------------------------------------------------
rnaseq_counts <- combined_counts %>%
  filter(seq_type == "RNAseq") %>%
  mutate(sample = toupper(sample)) %>%
  pivot_wider(id_cols = c(sample, seq_type), 
                      names_from = stage, 
                      values_from = raw_read_number)

paste0("___________________________________________________________________________________")
paste0("                        RNA-seq reads per sample                                   ")
paste0("___________________________________________________________________________________")
paste0("Average of ", floor(mean(rnaseq_counts$passed)/1e6), " million reads passed preprocesing for each RNA-seq sample.")


## ----supp_data_2-----------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 2.xlsx", sep = "")

suppressMessages(if (file.exists(fn)) {
  file.remove(fn)})

riboseq_output_counts <- riboseq_counts %>% 
  dplyr::rename(`cutadapt removed` = "Cutadapt",
         `rRNA mapped` = "rRNA",
         `snoRNA mapped` = "snoRNA",
         `tRNA mapped` = "tRNA")

rnaseq_output_counts <- rnaseq_counts %>%
  select(sample, -seq_type, raw, passed) %>%
  dplyr::rename("preprocessed" = passed, Raw=raw)

write_xlsx(list(a = riboseq_output_counts %>% filter(seq_type=="Cycloheximide") %>% 
                  select(-seq_type), 
                b = riboseq_output_counts %>% filter(seq_type=="Harringtonine") %>% 
                  select(-seq_type), 
                c = riboseq_output_counts %>% filter(seq_type=="No drug") %>% 
                  select(-seq_type), 
                d = rnaseq_output_counts),
  path = fn,
  format_headers = TRUE)


## ----fig_1d----------------------------------------------------------------------------------------
length_distribution_files <- fs::dir_ls("sequencing/stats/length_distribution/", 
                        regexp = "\\.txt$")

length_distirbution_data <- length_distribution_files %>%
  map_dfr(read_table2, .id = "source", col_names = FALSE, col_types = cols()) %>%
  mutate(ribotype=case_when(
    str_detect(source,"riboseq_chx") ~ "riboseq_chx",
    str_detect(source,"riboseq_harr") ~ "riboseq_harr",
    str_detect(source,"riboseq_nd") ~ "riboseq_nd",
    str_detect(source,"rnaseq_se") ~ "RNA-seq")) %>%
  mutate(sample=case_when(
    str_detect(source,"nts_r1") ~ "NTS_R1",
     str_detect(source,"nts_r2") ~ "NTS_R2",
     str_detect(source,"nts_r3") ~ "NTS_R3",
     str_detect(source,"nts_r4") ~ "NTS_R4",
       str_detect(source,"ts_r1") ~ "TS_R1",
     str_detect(source,"ts_r2") ~ "TS_R2",
     str_detect(source,"ts_r3") ~ "TS_R3",
     str_detect(source,"ts_r4") ~ "TS_R4",
    str_detect(source,"merged") ~ "rnaseq")) %>%
  mutate(readlength=X2,
         count=X1) %>%
  dplyr::select(sample,ribotype,readlength, count)

labels <- c(riboseq_chx = "Ribo-seq (CHX)", riboseq_harr = "Ribo-seq (HARR)", riboseq_nd = "Ribo-seq (ND)", `RNA-seq` = "RNA-seq")

dist_plot_data <- length_distirbution_data %>%
  group_by(ribotype,readlength) %>%
  summarise(count=sum(count))

read_lengths_hist_plot <- dist_plot_data %>%
  filter(readlength >= 25 & readlength <= 34) %>%
  mutate(ribotype = factor(ribotype, levels = c("riboseq_harr", "riboseq_chx", "riboseq_nd", "RNA-seq"))) %>%
  group_by(ribotype) %>%
  mutate(total = count) %>%
  ggplot(aes(x = readlength, y = total, fill = ifelse(readlength > 27 & readlength < 32 & ribotype != "RNA-seq", "RPF", "Non-RPF"))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ribotype, ncol = 2, labeller = labeller(ribotype = labels), scales = "free_x") +
  labs(y = "Number of reads", x = "Read length (nts)", fill = "") +
  scale_y_continuous(breaks = c(20000000, 40000000, 60000000), labels = c("20 million", "40 million", "60 million")) +
  scale_x_continuous(breaks = 25:34) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), legend.position = "bottom",axis.text.x= element_text(size=7)) +
  scale_fill_manual(values = c("#999999", "#0C7BDC"))


## ----phasing_analysis------------------------------------------------------------------------------
phasing_path <- "sequencing/plastid_analysis/merged/"

riboseq_chx_phasing <- read_delim(paste0(phasing_path, "riboseq_chx_phasing.txt"),
  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, show_col_types = FALSE) %>%
  summarise(
    `0` = mean(phase0),
    `1` = mean(phase1),
    `2` = mean(phase2)
  ) %>%
  mutate(ribotype = "Ribo-seq (CHX)")

riboseq_harr_phasing <- read_delim(paste0(phasing_path, "riboseq_harr_phasing.txt"),
  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, show_col_types = FALSE) %>%
  summarise(
    `0` = mean(phase0),
    `1` = mean(phase1),
    `2` = mean(phase2)
  ) %>%
  mutate(ribotype = "Ribo-seq (HARR)")

riboseq_nd_phasing <- read_delim(paste0(phasing_path, "riboseq_nd_phasing.txt"),
  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, show_col_types = FALSE) %>%
  summarise( 
    `0` = mean(phase0),
    `1` = mean(phase1),
    `2` = mean(phase2)
  ) %>%
  mutate(ribotype = "Ribo-seq (ND)")

rnaseq_phasing <- read_delim(paste0(phasing_path, "rnaseq_se_phasing.txt"),
  "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, show_col_types = FALSE) %>%
  summarise(
    `0` = mean(phase0),
    `1` = mean(phase1),
    `2` = mean(phase2)
  ) %>%
  mutate(ribotype = "RNA-seq")

paste0("___________________________________________________________________________________")
paste0("                        Cycloheximide Ribo-seq 28-31nt Phasing                     ")
paste0("___________________________________________________________________________________")
paste0("Approx. ", round(riboseq_chx_phasing$`0`*100), "% of reads are in phase")

paste0("__________________________________________________________________________________ ")
paste0("                        Harringtonine Ribo-seq 28-31nt Phasing                     ")
paste0("___________________________________________________________________________________")
paste0("Approx. ", round(riboseq_harr_phasing$`0`*100), "% of reads are in phase")

paste0("__________________________________________________________________________________ ")
paste0("                           No drug Ribo-seq 28-31nt Phasing phasing                ")
paste0("___________________________________________________________________________________")
paste0("Approx. ", round(riboseq_nd_phasing$`0`*100), "% of reads are in phase")

phasing_plot <- bind_rows(riboseq_chx_phasing, riboseq_harr_phasing, riboseq_nd_phasing, rnaseq_phasing) %>%
  pivot_longer(cols = c(`0`, `1`, `2`), values_to = "Proportion") %>%
  mutate(ribotype = factor(ribotype, levels = c("Ribo-seq (HARR)", "Ribo-seq (CHX)", "Ribo-seq (ND)", "RNA-seq"))) %>%
  ggplot(aes(x = name, y = Proportion, fill = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ribotype, ncol = 2) +
  labs(x = "Frame", y = "Proportion of RPFs", fill = "Frame") +
  scale_fill_viridis(discrete = T, option = "C") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), legend.position = "bottom")


## ----metagene_plot---------------------------------------------------------------------------------
chx_metagene <- read_tsv("orf_identification/orfrater/chx/metagene.txt", show_col_types = FALSE) %>%
  mutate(ribotype = "chx")

harr_metagene <- read_tsv("orf_identification/orfrater/harr/metagene.txt", show_col_types = FALSE) %>%
  mutate(ribotype = "harr")

nd_metagene <- read_tsv("orf_identification/orfrater/nd/metagene.txt", show_col_types = FALSE) %>%
  mutate(ribotype = "nd")

metagene <- bind_rows(chx_metagene, harr_metagene, nd_metagene)

metagene <- metagene %>%
  mutate(Means = rowMeans(.[, 3:6]))

metagene <- metagene %>%
  mutate(Frame = case_when(
    position %% 3 == 0 ~ "0",
    position %% 3 == 1 ~ "1",
    position %% 3 == 2 ~ "2"
  ))

metagene$region <- factor(metagene$region, # Reordering group factor levels
  levels = c("START", "CDS", "STOP")
)


metagene <- metagene %>%
  filter(region == "START") %>%
  filter(position <= 50) %>%
  mutate(ribotype = case_when(
    ribotype == "chx" ~ "Ribo-seq (CHX)",
    ribotype == "harr" ~ "Ribo-seq (HARR)",
    ribotype == "nd" ~ "Ribo-seq (ND)"
  ))
  
 metagene <- metagene %>% 
  mutate(translation_type = case_when(
    ribotype == "Ribo-seq (CHX)" ~ "Ribo-seq (CHX & ND)",
       ribotype == "Ribo-seq (ND)" ~ "Ribo-seq (CHX & ND)",
    ribotype == "Ribo-seq (HARR)" ~ "Ribo-seq (HARR)"
  ))

labels <- c(chx = "Cycloheximide", harr = "Harringtonine", nd = "No drug")

metagene_plot <- metagene %>%
    group_by(translation_type, position) %>%
  summarise(mean_at_position=mean(Means)) %>%
  mutate(Frame = case_when(
    position %% 3 == 0 ~ "0",
    position %% 3 == 1 ~ "1",
    position %% 3 == 2 ~ "2"
  )) %>%
    mutate(translation_type = factor(translation_type, levels = c("Ribo-seq (HARR)", "Ribo-seq (CHX & ND)"))) %>%
    ggplot(aes(x = position, y = mean_at_position,  fill=Frame)) +
    geom_bar(stat="identity") +
  facet_wrap(~translation_type) +
  labs(y = "Averge RPF density", color = "Frame", x = "Transcript position") +
  theme_bw() +
  scale_fill_viridis(discrete = T, option = "C") +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  theme(strip.background = element_rect(fill = "white"), legend.position = "bottom")
  
read_lengths_hist_plot + 
  plot_spacer() + 
  phasing_plot + 
  plot_spacer() + 
  metagene_plot + 
  plot_layout(widths = c(5, 0.25 ,5, 0.25, 10))

  ggsave(filename = paste(results_dir, "Figure 1def.png", sep = ""), 
         width = 15, height = 4, device = "png", dpi = 2000)

metagene %>%   
  mutate(ribotype = factor(ribotype, levels = c("Ribo-seq (HARR)", "Ribo-seq (CHX)", "Ribo-seq (ND)"))) %>%
  ggplot(aes(x = position, y = Means, fill = Frame), alpha = ribotype) +
  geom_bar(stat="identity") +
  labs(y = "Averge RPF density", color = "", x = "Transcript position") +
  facet_wrap(~ribotype, ncol=1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_viridis(discrete = T, option = "C") +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  theme(strip.background = element_rect(fill = "white"), legend.position = "bottom")


ggsave(filename = paste(results_dir, "Supplementary Figure 3.png", sep = ""), 
       width = 6, height = 8, device = "png", dpi = 2000)


## --------------------------------------------------------------------------------------------------
print("section 2.1 results complete")

