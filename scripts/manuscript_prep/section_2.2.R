## ------------------------------------------------------------------------------------------------------------
source(file = "scripts/manuscript_prep/utility_functions.R")


## ----setup---------------------------------------------------------------------------------------------------
package_list <- c("tidyverse","gridExtra","viridis", "ggforce", 
                  "ggpp", "cowplot", "writexl", "patchwork")

lapply(package_list, require, character.only = TRUE)

options(dplyr.summarise.inform = FALSE)


## ----make_results_dir----------------------------------------------------------------------------------------
results_dir <- "manuscript/section_2.2/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ------------------------------------------------------------------------------------------------------------
paste0("___________________________________________________________________________________")
paste0("                 Pseduogenes  Transcript removal prior to ORF-RATER.               ")
paste0("___________________________________________________________________________________")
pseudogene_gene_ids <- read_csv("orf_identification/pseudogene_gene_ids.txt", 
                                col_names = F, show_col_types = F)

paste0(length(unique(pseudogene_gene_ids$X1)), " pseudogenes removed")

paste0("___________________________________________________________________________________")
paste0("                 Transcripts removed by ORF-RATER.               ")
paste0("___________________________________________________________________________________")
tid_removal_summary <- read_delim("orf_identification/orfrater/tid_removal_summary.txt", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE, show_col_types = F) %>%
  filter(!is.na(dropped))

table(tid_removal_summary$dropped)


## ----initial_orfrater_result---------------------------------------------------------------------------------
orfrater_dir=paste0(getwd(), "/orf_identification/orfrater/")

orftable_raw <- read.csv(file = paste0(orfrater_dir,"rate_regression.csv"), 
                         header = T) %>%
  filter(orfrating >= 0.5 & AAlen >= 5)

paste0("___________________________________________________________________________________")
paste0("                              Raw ORFs identified                                  ")
paste0("___________________________________________________________________________________")
paste0(dim(orftable_raw)[1], " ORFs with an orfrating >= 0.5 and longer than 5aa found by ORFRATER")

paste0("___________________________________________________________________________________")
paste0("                              Removed ORF classes                                  ")
paste0("___________________________________________________________________________________")
print(orftable_raw %>%
  filter(orftype %in% c("Ciso", "Giso", "NCiso", "Niso", 
  "new_iso", "internal", "LOOF", "truncation")) %>%
  group_by(orftype) %>%
  summarise(removed=n()) %>%
    arrange(-removed))



paste0("___________________________________________________________________________________")
paste0("                              TIS filter removed                                   ")
paste0("___________________________________________________________________________________")
orfs_post_filter1 <- readRDS(file = "orf_identification/orf_filtered/orfs_remaining_post_first_filter.rds")
orfs_post_tis_filter <- readRDS(file = "orf_identification/orf_filtered/orfs_remaining_post_tis_filter.rds")
paste0(dim(orfs_post_filter1)[1] - dim(orfs_post_tis_filter)[1], " Removed by the TIS filter")

paste0("___________________________________________________________________________________")
paste0("                              Overlap filter removed                                   ")
paste0("___________________________________________________________________________________")
orfs_post_overlap_filter <- readRDS(file = "orf_identification/orf_filtered/orfs_remaining_post_overlap_filter.rds")
paste0(dim(orfs_post_tis_filter)[1] - dim(orfs_post_overlap_filter)[1] , " Removed by the overlap filter")


## ----orf_rater_summary---------------------------------------------------------------------------------------
orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")

paste0("___________________________________________________________________________________")
paste0("                              FLOSS filter removed                                   ")
paste0("___________________________________________________________________________________")
paste0(dim(orfs_post_overlap_filter)[1] - dim(orftable)[1], " Removed by the FLOSS filter")

paste0("___________________________________________________________________________________")
paste0("                              Total                                                ")
paste0("___________________________________________________________________________________")
paste0(dim(orftable)[1], " ORFs confidently identified")

paste0("___________________________________________________________________________________")
paste0("                              All ORF start codon                                  ")
paste0("___________________________________________________________________________________")
table(orftable$`ORF type`, orftable$`Start codon`)

paste0("___________________________________________________________________________________")
paste0("                              Novel                                                ")
paste0("___________________________________________________________________________________")
novel_orfs <- orftable %>%
  filter(`ORF type` != "Annotated") %>%
  summarise(total=length(`ORF-RATER name`))

paste0(novel_orfs$total, " Previously unnannotated Chinese hamster ORFs identified")

paste0("___________________________________________________________________________________")
paste0("                              Novel ORF start codons                               ")
paste0("___________________________________________________________________________________")
percent_start_codon <- orftable %>%
  filter(`ORF type` != "Annotated") %>%
  group_by(`Start codon`) %>%
  summarise(
    count = n(),
    percentage = (n() / nrow(.)) * 100
  ) %>%
  arrange(-percentage)
paste0(floor(percent_start_codon[1,2]), " or ", round(percent_start_codon[1,3],2), "% of novel ORFs identified had an ",percent_start_codon[1,1], " start codon" )
paste0(floor(percent_start_codon[2,2]), " or ", round(percent_start_codon[2,3],2), "% of novel ORFs identified had a ",percent_start_codon[2,1], " start codon" )
paste0(floor(percent_start_codon[3,2]), " or ", round(percent_start_codon[3,3],2), "% of novel ORFs identified had a ",percent_start_codon[3,1], " start codon" )
paste0(floor(percent_start_codon[4,2]), " or ", round(percent_start_codon[4,3],2), "% of novel ORFs identified had a ",percent_start_codon[4,1], " start codon" )

paste0("___________________________________________________________________________________")
paste0("                              Novel ORF types.                                     ")
paste0("___________________________________________________________________________________")
percent_orftype <- orftable %>%
  filter(`ORF type` != "Annotated") %>%
  group_by(`ORF type`) %>%
  summarise(
    count = n(),
    percentage = (n() / nrow(.)) * 100
  ) %>%
  arrange(-percentage)

paste0(floor(percent_orftype[1,2]), " or ", round(percent_orftype[1,3],2), "% of novel ORFs were ",percent_orftype[1,1], " ORFs" )
paste0(floor(percent_orftype[2,2]), " or ", round(percent_orftype[2,3],2), "% of novel ORFs were ",percent_orftype[2,1], " ORFs" )
paste0(floor(percent_orftype[3,2]), " or ", round(percent_orftype[3,3],2), "% of novel ORFs were ",percent_orftype[3,1], " ORFs" )
paste0(floor(percent_orftype[4,2]), " or ", round(percent_orftype[4,3],2), "% of novel ORFs were ",percent_orftype[4,1], " ORFs" )
paste0(floor(percent_orftype[5,2]), " or ", round(percent_orftype[5,3],2), "% of novel ORFs were ",percent_orftype[5,1], " ORFs" )
paste0(floor(percent_orftype[6,2]), " or ", round(percent_orftype[6,3],2), "% of novel ORFs were ",percent_orftype[6,1], " ORFs" )
paste0(floor(percent_orftype[7,2]), " or ", round(percent_orftype[7,3],2), "% of novel ORFs were ",percent_orftype[7,1], " ORFs" )


## ------------------------------------------------------------------------------------------------------------
orfrater_out <- orftable_raw  %>%
filter(orftype %in% c("extension", "upstream", "start_overlap", "new")) %>%
  dplyr::rename(`ORF type` = orftype, 
                `Start codon` = codon) %>%
  mutate(`ORF type` = case_when(
    `ORF type` == "upstream" ~ "Upstream",
    `ORF type` == "start_overlap" ~ "Start overlap",
    `ORF type` == "extension" ~ "Extension",
    `ORF type` == "new" ~ "New",
  ))

orfrater_out <- orfrater_out %>%
  mutate(stage="1. ORF-RATER output") %>%
  dplyr::select(stage,`ORF type`,`Start codon`) %>%
  mutate(`Start codon` = case_when(
   `Start codon` == "ATG" ~ "AUG", 
   `Start codon` == "CTG" ~ "CUG", 
   `Start codon` == "GTG" ~ "GUG",
   `Start codon` == "TTG" ~ "UUG")) %>%
  dplyr::select(stage,`ORF type`,`Start codon`)

tzani_orf_data <- orftable %>%
  filter(`ORF type` %in% c("Extension", "Upstream", "Start overlap", "New"))

tzani <- orftable %>%
  mutate(stage="4. FLOSS filter") %>%
  filter(`ORF type` %in% c("Extension", "Upstream", "Start overlap", "New")) %>%
  dplyr::select(stage,`ORF type`,`Start codon`) %>%
  mutate(`Start codon` = case_when(
   `Start codon` == "ATG" ~ "AUG", 
   `Start codon` == "CTG" ~ "CUG", 
   `Start codon` == "GTG" ~ "GUG",
   `Start codon` == "TTG" ~ "UUG")) %>%
  dplyr::select(stage,`ORF type`,`Start codon`)

tzani_tis_filter <- orfs_post_overlap_filter %>%
  filter(`ORF type` %in% c("Extension", "Upstream", "Start overlap", "New")) %>%
  mutate(stage="2. TIS filter") %>%
  dplyr::select(stage,`ORF type`,`Start codon`) %>%
  mutate(`Start codon` = case_when(
   `Start codon` == "ATG" ~ "AUG", 
   `Start codon` == "CTG" ~ "CUG", 
   `Start codon` == "GTG" ~ "GUG",
   `Start codon` == "TTG" ~ "UUG")) %>%
  dplyr::select(stage,`ORF type`,`Start codon`)

tzani_overlap_filter <- orfs_post_overlap_filter %>%
  filter(`ORF type` %in% c("Extension", "Upstream", "Start overlap", "New")) %>%
  mutate(stage="3. Overlap filter") %>%
  dplyr::select(stage,`ORF type`,`Start codon`) %>%
  mutate(`Start codon` = case_when(
   `Start codon` == "ATG" ~ "AUG", 
   `Start codon` == "CTG" ~ "CUG", 
   `Start codon` == "GTG" ~ "GUG",
   `Start codon` == "TTG" ~ "UUG")) %>%
  dplyr::select(stage,`ORF type`,`Start codon`)

data_for_filter_plot <- bind_rows( orfrater_out, tzani_tis_filter, tzani_overlap_filter, tzani)

data_for_filter_plot <- data_for_filter_plot %>% 
  group_by(stage, `ORF type`, `Start codon`) %>%
  count() %>%
  group_by(stage, `ORF type`) %>%
  mutate(percentage = n / sum(n) * 100) 

data_for_filter_plot %>%
  ggplot(aes(x = `ORF type`, y = percentage,  fill = `Start codon`)) +
  geom_bar(stat="identity") +
  scale_fill_viridis(discrete = T, option = "D") +
  facet_wrap(~factor(stage, 
                     levels=c("1. ORF-RATER output", "2. TIS filter", "3. Overlap filter", "4. FLOSS filter")), 
             nrow = 1) +
  labs(x="", y= "% of ORFs") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult=c(0, 0.05))) +
  geom_text(aes(label = if_else(percentage >= 5, paste0(sprintf("%1.0f", percentage), "%"), "")),
            position=position_stack(vjust=0.5),
            colour=ifelse(data_for_filter_plot$`Start codon`=="UUG", 'black', 'white'), 
            size =3) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") 

 ggsave(filename = paste(results_dir, "Supplementary Figure 4.png", sep = ""), 
        device = "png", dpi = 2000, width = 8, height = 3.75)


## ----eval=FALSE, include=FALSE-------------------------------------------------------------------------------
## uORF_dump_uORFdb <- read_delim("~/Downloads/uORF_dump_uORFdb.tsv",
##     delim = "\t", escape_double = FALSE,
##     trim_ws = TRUE)
## 
## saveRDS(uORF_dump_uORFdb, "uORFdb_dump")


## ------------------------------------------------------------------------------------------------------------


uORF_dump_uORFdb <- readRDS("uORFdb_dump")

species_list <-c(
  "Mus musculus", "Sus scrofa", "Macaca mulatta", "Pongo abelii", "Gallus gallus",
  "Pan troglodytes", "Homo sapiens", "Rattus norvegicus", "Bos taurus",  
  "Xenopus laevis", "Drosophila melanogaster", "Gallus gallus"
)

start_codons <- c("ATG", "CTG", "GTG", "TTG")

orfs <- c("overlapping", "N-terminal extension", "non-overlapping")

selected_uORFdb <- uORF_dump_uORFdb %>% 
  filter(Taxon %in% species_list) %>%
  filter(uORFstartCodon %in% start_codons) %>%
  filter(uORFtype %in% orfs) 

uORFdb <- selected_uORFdb %>% 
  mutate(`ORF type` = case_when(
  uORFtype == "non-overlapping" ~ "Upstream", 
  uORFtype == "overlapping" ~ "Start overlap", 
  uORFtype == "N-terminal extension" ~ "Extension")) %>%
  mutate(`Start codon` = case_when(
   uORFstartCodon == "ATG" ~ "AUG", 
   uORFstartCodon == "CTG" ~ "CUG", 
   uORFstartCodon == "GTG" ~ "GUG",
   uORFstartCodon == "TTG" ~ "UUG")) %>%
  dplyr::select(Taxon,`ORF type`,`Start codon`)

tzani <-orftable %>%
  mutate(Taxon="Cricetulus griseus") %>%
  filter(`ORF type` %in% c("Extension", "Upstream", "Start overlap")) %>%
  dplyr::select(Taxon,`ORF type`,`Start codon`) %>%
  mutate(`Start codon` = case_when(
   `Start codon` == "ATG" ~ "AUG", 
   `Start codon` == "CTG" ~ "CUG", 
   `Start codon` == "GTG" ~ "GUG",
   `Start codon` == "TTG" ~ "UUG")) %>%
  dplyr::select(Taxon,`ORF type`,`Start codon`)

data_for_species_plot <- bind_rows(uORFdb, tzani)

data_for_species_plot <- data_for_species_plot %>%
  group_by(Taxon, `ORF type`, `Start codon`) %>%
  count() %>%
  group_by(Taxon, `ORF type`) %>%
  mutate(percentage = n / sum(n) * 100) 
  # filter(percentage >2)

data_for_species_plot %>%
  ggplot(aes(x = `ORF type`, y = percentage,  fill = `Start codon`)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = if_else(percentage >= 5, paste0(sprintf("%1.0f", percentage), "%"), "")),
            position=position_stack(vjust=0.5),
             colour=ifelse(data_for_species_plot$`Start codon`=="UUG", 'black', 'white'), 
            size =3) +
  scale_fill_viridis(discrete = T, option = "D") +
  facet_wrap(~factor(Taxon, levels=c("Cricetulus griseus","Mus musculus","Rattus norvegicus","Homo sapiens",
                                     "Sus scrofa", "Macaca mulatta", "Pongo abelii", "Gallus gallus",
                                     "Xenopus laevis", "Drosophila melanogaster", "Pan troglodytes", "Bos taurus"))) +
  labs(x="", y= "% of ORFs") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "italic"), 
        axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") 

 ggsave(filename = paste(results_dir, "Supplementary Figure 5.png", sep = ""), 
        device = "png", dpi = 2000, width = 8.5, height = 10)



## ----supp_data_3---------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 3.xlsx", sep = "")

suppressMessages(if (suppressMessages(file.exists(fn))) {
  file.remove(fn)})

orf_output_list <- orftable %>%
  dplyr::select(`ORF-RATER name`:`Genomic stop position`,
                Rharr_min_Rnd, countRFP, fpkmRFP:floss, floss_classification) %>%
  mutate(output_group= case_when(
    `ORF type` == "Annotated" ~ "Annotated",
    `ORF type` == "Upstream" |`ORF type` == "Start overlap" ~ "Upstream",
    `ORF type` == "New" ~ "New ORFs",
    `ORF type` == "Extension" ~ "N-terminal extension",
    `ORF type` == "Isoform" ~ "Isoform",
    `ORF type` == "Stop overlap" | `ORF type` == "Downstream" ~ "Other"
  )) %>%
  dplyr::rename(`CHX RFP count` = countRFP,`CHX RFP FPKM` = fpkmRFP,`RNA FPKM` = fpkmRNA) %>%
  group_split(output_group)


write_xlsx(list(
  a = orf_output_list[[1]] %>% dplyr::select(-output_group),
  b = orf_output_list[[6]] %>% dplyr::select(-output_group),
  c = orf_output_list[[4]] %>% dplyr::select(-output_group),
  d = orf_output_list[[3]] %>% dplyr::select(-output_group),
  e = orf_output_list[[2]] %>% dplyr::select(-output_group),
  f = orf_output_list[[5]] %>% dplyr::select(-output_group)
  ), path = fn, format_headers = TRUE)


## ----fig_2a--------------------------------------------------------------------------------------------------
orftype_codon_data <- data.frame(table(orftable$`ORF type`, orftable$`Start codon`)) 


orftype_codon_data$Var2 <- factor(orftype_codon_data$Var2, 
                                                    levels = c("ATG", "CTG", "GTG", "TTG"))

orftype_codon_data <- orftype_codon_data %>%
  mutate(start_codon = gsub("T", "U", Var2)) %>%
  group_by(start_codon) %>%
  mutate(count_name_occur = Freq) %>%
  group_by(Var1) %>%
  mutate(total = sum(Freq)) 
  

summary_start <- table(orftable$`ORF type`, orftable$`Start codon`)
colnames(summary_start) <- gsub("T", "U", colnames(summary_start))
summary_start <- summary_start[order(rowSums(-summary_start)), ]

sum_table <- tableGrob(summary_start, 
                       theme = ttheme_default(base_size = 12, 
                                              padding = unit(c(1, 1), 
                                                             "mm")))

orftype_codon_data$`start_codon` <- factor(orftype_codon_data$`start_codon`, 
                                                    levels = c("AUG", "CUG", "GUG", "UUG"))
orftype_codon_data %>%
  arrange(-start_codon) %>%
  mutate(start_codon = factor(start_codon, levels = c("UUG", "GUG", "CUG", "AUG"))) %>%
  ggplot(aes(x = reorder(Var1, count_name_occur), y = Freq, fill = start_codon)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = T, direction = -1, option="D") +
  theme_minimal_vgrid() +
  scale_y_continuous(limits = c(0,5800),expand = c(0, 0), breaks=seq(0, 5800, by = 1000)) +
  coord_flip() +
  labs(y = "# Identified ORFs", x = "", fill = "Start codon") +
  theme(legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotation_custom(sum_table, xmin = 5, xmax = 1, ymin = 3600, ymax = 4400) +
  guides(fill = guide_legend(reverse=TRUE))

ggsave(filename = paste(results_dir, "Figure 2a.png", sep = ""), 
       dpi = 2000, device = "png", height =3.5, width = 6)


## ----fig2_b--------------------------------------------------------------------------------------------------
# set path to coverage
number_nt <- 1400

coverage_path="sequencing/coverage_tracks/merged/fig2b/"

padding_start <- data.frame(nt_postion=c(1:46),
             count=rep(0,46))

padding_end <- data.frame(nt_postion=c(1673:1941),
             count=rep(0,269))


# panel 1: overview pane to illustrate annotation

# placeholder for plot
annotation_placeholder <- data.frame(x1=c(0, 2000),
           y1=c(0,2))

fig2b_annotation <- annotation_placeholder %>% 
  ggplot(aes(x=x1, y=y1)) +
  lims(x=c(0,number_nt))+
  # labs(title = "Aurora kinase A", subtitle="XM_027423276.2") +
  annotate("rect", 
           xmin=c(93,144), 
           xmax = c(1323, 1323), 
           ymin=c(1.01, 0.01), 
            ymax=c(2,0.99), alpha = .7, 
            fill = c("#88CCEE", "#CC79A7")) +
  annotate("text", x=c(480, 430), 
           y=c(1.5, 0.5), 
           label = c("409 aa N-terminally extended Aurka ORF", 
                     "392 aa canonical Aurka ORF"), size=3.5, color="black") +
  theme(validate = FALSE, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), 
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
    plot.background = element_rect(fill='transparent', color=NA))

# panel 2: RNA-seq
rnaseq_data_full <- read_delim(paste0(coverage_path,"XM_027423276.2_rnaseq_f_transcriptome.wig"),
                         delim = "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, show_col_types = FALSE)  %>%
  dplyr::select(X1, X2,X4) %>%
  dplyr::rename(transcript=X1, position=X2, coverage=X4) 

# need to make a new row, because deeptools bins adjoining nt positions with an identifical value
# makes the region look aritifcally uncovered when plotting with ggplot
missing_row = data.frame(transcript = c("XM_027423276.2", "XM_027423276.2", "XM_027423276.2"),
                         position = c(219,482, 483),
                         coverage = c(0.722197, 1.36022,1.36022))

# scale coverage
rnaseq_data_full <- rnaseq_data_full %>%
  bind_rows(missing_row) %>%
  mutate(norm=coverage/sum(coverage),
         coverage_type = "rna_full")  %>%
  mutate(scaled_coverage = scale_between_0_1(coverage)) 

fig2b_rnaseq <- rnaseq_data_full %>%
  ggplot(aes(x = position, y = scaled_coverage, alpha = as.character(coverage_type))) +
  geom_bar(stat = "identity", position = "identity") +
   geom_vline(xintercept = c(83, 153), linewidth=0.1, linetype="dashed") +
  scale_fill_viridis() +
  scale_alpha_discrete(range = c(0.3, 1)) +
  labs(
    x = "Transcript coordinates (nts)",
    y = "coverage (scaled BPM)",
    title = "",
    subtitle = ""
  ) +
  lims(x = c(0,number_nt)) +
  theme_bw() +
  theme(validate = FALSE, legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill='transparent'), 
    plot.background = element_rect(fill='transparent', color=NA)) + 
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "RNA", size = 2)

# panel 3: CHX Ribo-seq
aurka_extension <- orftable %>%
    filter(`Transcript ID` == "XM_027423276.2" & `ORF type` == "Extension")

chx_data_psite <- read_delim(paste0(coverage_path, "XM_027423276.2_chx_p_transcriptome.wig"), 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)  %>%
  dplyr::rename(nt_position=X2, count=X4) %>%
  mutate(norm=count/sum(count),
         coverage_type = "chx_psite") %>%
  mutate(frame = "NA") %>%
  mutate(frame = ifelse((nt_position >= min(aurka_extension$`Transcript start position`) & nt_position <= aurka_extension$`Transcript stop position` & coverage_type == "chx_psite"), 
                        (nt_position - aurka_extension$`Transcript start position`) %% 3, NA)) %>%
  mutate(scale_bpm = scale_between_0_1(count))

chx_data_full <- read_delim(paste0(coverage_path, "XM_027423276.2_chx_f_transcriptome.wig"), 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE, show_col_types = FALSE)  %>%
  dplyr::rename(nt_position=X2, count=X4) %>%
  mutate(norm=count/sum(count),
         coverage_type = "chx_full")  %>%
  mutate(scale_bpm = scale_between_0_1(count)) %>%
  mutate(frame = NA)
 
fig2b_chx_riboseq <- bind_rows(chx_data_full,chx_data_psite) %>%
  ggplot(aes(x = nt_position, y = scale_bpm, fill = frame, alpha = as.character(coverage_type))) +
   geom_vline(xintercept = c(83, 153), linewidth=0.1, linetype="dashed") +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_viridis(option = "C") +
  scale_alpha_discrete(range = c(0.3, 1)) +
  labs(
    x = "Transcript coordinates (nts)",
    y = "coverage (scaled BPM)",
    title = "",
    subtitle = ""
  ) +
  lims(x = c(0,number_nt)) +
  theme_bw() +
  theme(validate = FALSE, legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill='transparent', color=NA)) + 
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "CHX", size = 2)


# panel 4: HARR-ND Ribo-seq
harr_coverage  <- read_delim(paste0(coverage_path, "XM_027423276.2_harr_p_transcriptome.wig"), 
                       delim = "\t", escape_double = FALSE, 
                       col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)  %>%
  dplyr::rename(nt_position=X2, count=X4) %>%
bind_rows(padding_start) %>%
bind_rows(padding_end) %>%
mutate(norm=count/sum(count))

nd_coverage  <- read_delim(paste0(coverage_path,"XM_027423276.2_nd_p_transcriptome.wig"), 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
dplyr::rename(nt_position=X2, count=X4) %>%
mutate(norm=count/sum(count))

merged_data <- merge(harr_coverage, nd_coverage, 
                     by = c("X1", "nt_position", "X3"), all = TRUE)

merged_data[is.na(merged_data)] <- 0  # Normalize Harringtonine signal by Cycloheximide signal 
merged_data$harr_nd <- (merged_data$norm.x - merged_data$norm.y)
merged_data$harr_nd[merged_data$harr_nd < 0] <- 0 
  

merged_data$harr_nd_scaled <- scale_between_0_1(merged_data$harr_nd)

fig2b_harr_full <- merged_data %>%
ggplot(aes(x=nt_position, y=harr_nd_scaled)) +
  geom_bar(stat="identity", position = "identity", alpha = 1, fill = "#404788FF") +
  geom_vline(xintercept = c(83, 153), linewidth=0.1, linetype="dashed") +
  labs(
    x = "Transcript coordinates (nts)",
    y = "coverage (scaled BPM)",
    title = "",
    subtitle = ""
  ) +
  lims(x = c(0,number_nt)) +
  theme_bw() +
  theme(validate = FALSE, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
       panel.background = element_rect(fill="transparent"), 
    plot.background = element_rect(fill='transparent', color=NA)) + 
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "HARR-ND", size = 2) 
  
  
fig2b_harr_tis <- merged_data %>%
  filter(nt_position >= 83 & nt_position <=153) %>%
ggplot(aes(x=nt_position, y=harr_nd_scaled)) +
  geom_bar(stat="identity", position = "identity", alpha = 1, fill = "#404788FF") +
   geom_vline(xintercept = c(83, 153), linewidth=0.1, linetype="dashed") +
  lims(x = c(83,153)) +
  labs(
    x = "Transcript coordinates (nts)",
    y = "coverage (scaled BPM)",
    title = "",
    subtitle = ""
  ) +
  theme_bw() +
  theme(validate = FALSE, 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="transparent"), 
    plot.background = element_rect(fill='transparent', color=NA)) + 
  geom_text_npc(aes(npcx = 0.95, npcy = 0.95), label = "HARR-ND", size = 2) 

(fig2b_annotation + theme(plot.margin = unit(c(0,0,0,0), "cm")))/
  plot_spacer() /
  (fig2b_rnaseq + theme(plot.margin = unit(c(0,0,0,0), "pt")))/
    plot_spacer() /
  (fig2b_chx_riboseq + theme(plot.margin = unit(c(0,0,0,0), "pt")))/
    plot_spacer() /
  (fig2b_harr_full + theme(plot.margin = unit(c(0,0,0,0), "pt")))/
    plot_spacer() /
  (fig2b_harr_tis + theme(plot.margin = unit(c(0,0,0,0), "pt"))) +
  plot_layout(heights = c(0.25,-0.30, 0.60,-0.30,0.60,-0.30,0.60,-0.1,1), axes = "collect_y")

ggsave(filename = paste(results_dir, "Figure 2b.png", sep = ""), 
       dpi = 2000, device = "png", height =8, width = 6)


## ------------------------------------------------------------------------------------------------------------
print("section 2.2 results complete")

