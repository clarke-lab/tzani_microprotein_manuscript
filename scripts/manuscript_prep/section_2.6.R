## ----include=FALSE-------------------------------------------------------------------------------------------------------------
package_list <- c("tidyverse", "proDA", "scales", "pheatmap", "patchwork",
  "readxl","writexl", "viridis", "cowplot", "ggvenn", "ggpubr")

lapply(package_list, require, character.only = TRUE)




## ----include=FALSE-------------------------------------------------------------------------------------------------------------
results_dir <- "manuscript/section_2.6/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ----supp_figure_12------------------------------------------------------------------------------------------------------------
ts_cell_density <- read_tsv("data/cell_density/ts_cell_density_proteomics.txt", 
                            show_col_types = FALSE) 

ggboxplot(ts_cell_density , x = "condition", y = "cell_density", 
                               color = "black", add = "jitter", fill = "condition") +
  theme_bw() +
  labs(x = "", y = "Cell Density (cells/ml)", fill = "") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#56B4E9")) +
  scale_x_discrete(labels = c("NTS", "TS 24 h", "TS 48 h"))

ggsave(filename = paste(results_dir, "Supplementary Figure 12a.png", sep = ""), 
       width = 5, height = 4, device = "png", dpi = 2000)

growth_rates <- read_tsv("data/cell_density/exponential_v_stationary_proteomics.txt", show_col_types = FALSE) %>%
  pivot_longer(names_to = "Day", cols = c("Day 1", "Day 2", "Day 3", 
                                          "Day 4", "Day 5", "Day 6", "Day 7", ))

growth_rates %>%
  ggplot(aes(x = Day, y = value, fill = Day)) +
  geom_vline(xintercept = c("Day 4", "Day 7"), linetype = 11, linewidth = 0.2) +
  geom_boxplot() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_fill_manual(values = c("#999999", "#999999", "#999999", 
                               "#8785BA", "#999999", "#999999", 
                               "#5AAA46")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Density (cells/ml)")

ggsave(filename = paste0(results_dir,"Supplementary Figure 12b.png"), 
       device = "png", dpi = 2000, width = 4.5, height = 3)




## ----supp_table_1bc------------------------------------------------------------------------------------------------------------
ts_density <- ts_cell_density %>%
  pivot_wider(names_from= c(condition),values_from=cell_density) %>%
  rename(`NTS cell density (cells/ml)`=NTS, 
         `TS 24 h cell density (cells/ml)`=`TS 24hrs`,
         `TS 48 h cell density (cells/ml)`=`TS 48hrs`)

growth_rates <- read_tsv("data/cell_density/exponential_v_stationary_proteomics.txt", show_col_types = FALSE) %>%
  rename(`Day 1 cell density (cells/ml)`=`Day 1`,
         `Day 2 cell density (cells/ml)`=`Day 2`,
         `Day 3 cell density (cells/ml)`=`Day 3`,
         `Day 4 cell density (cells/ml)`=`Day 4`,
         `Day 5 cell density (cells/ml)`=`Day 5`,
         `Day 6 cell density (cells/ml)`=`Day 6`,
         `Day 7 cell density (cells/ml)`=`Day 7`)

fn <- paste(results_dir, "Supplementary Data 1bc.xlsx", sep = "")

suppressMessages(if (file.exists(fn)) {
  file.remove(fn)})


write_xlsx(list(`B`=ts_density,
                `C`=growth_rates),
path = fn, format_headers = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
canonical_metamorpheus <- readRDS(file = "proteomics/protein_identifications/canonical_proteins_lysate.rds") %>%
  mutate(class="canonical", type="Canonical") 

paste0("_________________________________________________________________________________ ")
paste0("                         Canonical proteins identified                            ")
paste0("__________________________________________________________________________________")
canonical_metamorpheus %>%
  distinct(experiment, `Protein Accession`) %>%
  group_by(experiment) %>%
  count()
  
orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")

orftable <- orftable %>%
  rename(protein = `ORF-RATER name`)

microproteins_pepquery <- readRDS(file = "proteomics/protein_identifications/microproteins_lysate.rds")

dp_microproteins <- readRDS("proteomics/protein_identifications/microproteins_drug_product.rds")


## ------------------------------------------------------------------------------------------------------------------------------
# shared peptides for overlap microproteins, rename here and include later
microproteins_pepquery <- microproteins_pepquery %>%
  mutate(original_protein = protein) %>%
  mutate(protein = case_when(
    str_detect(protein, "XR_003483809.2_106374381") ~ "XR_003483809.2_106374381_31aa",
    str_detect(protein, "XM_027427082.2_82745680") ~ "XM_027427082.2_82745680_37aa",
    str_detect(protein, "XM_027396479.2_395237785") ~ "XM_027396479.2_395237785_62aa",
    str_detect(protein, "XM_027414474.2_215763922") ~ "XM_027414474.2_215763922_39aa",
    str_detect(protein, "XM_027417997.2_89627633") ~ "XM_027417997.2_89627633_82aa", 
    TRUE ~ protein)) 
  
paste0("_________________________________________________________________________________ ")
paste0("                         Microproteins identified.                                ")
paste0("__________________________________________________________________________________")
microproteins_pepquery %>%
  distinct(experiment, protein) %>%
  group_by(experiment) %>%
  count()

paste0("_________________________________________________________________________________ ")
paste0("                         Microprotein ORF types.                                ")
paste0("__________________________________________________________________________________")
microproteins_pepquery %>%
  left_join(orftable, by="protein") %>%
  distinct(experiment, protein, `ORF type`) %>%
  group_by(experiment, `ORF type`) %>%
  count()


## ----fig6de--------------------------------------------------------------------------------------------------------------------
canonicals_for_barplot <- canonical_metamorpheus %>%
  distinct(experiment, `Protein Accession`, .keep_all = T) %>%
  dplyr::select(experiment,`Protein Accession`,class,type) %>%
  dplyr::rename(protein = `Protein Accession`)

microproteins_for_barplot <- microproteins_pepquery %>%
  distinct(experiment, protein) %>%
  left_join(orftable, by="protein")

microproteins_for_barplot <- microproteins_for_barplot %>%
  mutate(class="microprotein") %>%
  dplyr::select(experiment, protein, class, `ORF type`) %>%
  dplyr::rename(type = `ORF type`)

identified_proteins <- bind_rows(canonicals_for_barplot, microproteins_for_barplot)

protein_counts <- identified_proteins %>%
  group_by(experiment, class,type) %>%
  summarise(count=n())

bar_labels <- data.frame(experiment =protein_counts$experiment, 
                      label=protein_counts$count,  
                      yposition=c(4500,50,38,15,
                                  4500,102,70,25))

identified_proteins$experiment <- factor(identified_proteins$experiment, 
                                                     levels = c("tempshift", "d4d7"))

facet_order.labs <- c("Canonical protein", "Microprotein")
names(facet_order.labs) <- c("canonical", "microprotein")

protein_counts %>%
  ggplot(aes(x= experiment, y=count, fill = type)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  facet_wrap(~class, nrow = 1, 
             scales = "free_y", labeller = labeller(class= facet_order.labs)) +
  labs(y = "# Identified", x = NULL, fill = "") +
  scale_fill_manual(values = c("#CC6677","#35968CFF", "#332288","#999933")) +
    scale_x_discrete(limits = c("tempshift", "d4d7"), labels= c("Temperature\nShift", 
                                                                "Growth\nPhase")) +
  geom_text(aes(x = bar_labels$experiment,
                y = bar_labels$yposition, 
                label = bar_labels$label),color = "white", size = 3) +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 10), legend.position = "top")

ggsave(filename = paste0(results_dir,"Figure 6de.png"), 
       device = "png", dpi = 5000, width = 6, height = 4)


## ----fig6fg--------------------------------------------------------------------------------------------------------------------
microproteins_with_orf_info <- microproteins_pepquery %>%
  distinct(experiment, protein) %>%
  left_join(orftable, by = "protein") %>%
  group_by(experiment,`ORF type`, `Start codon`) %>%
  summarise(orf_count=n())

experiment.labs <- c("D4D7", "Temperature shift")
names(experiment.labs) <- c("d4d7", "tempshift")

microproteins_with_orf_info$`Start codon` <- factor(microproteins_with_orf_info$`Start codon`, 
                                                    levels = c("ATG", "CTG", "GTG", "TTG"))

microproteins_with_orf_info$facet_order = factor(microproteins_with_orf_info$experiment, levels=c('tempshift',"d4d7"))                      

facet_order.labs <- c("Growth phase", "Temperature shift")
names(facet_order.labs) <- c("d4d7", "tempshift")

microproteins_with_orf_info %>%
  mutate(`Start codon` = gsub("T", "U",`Start codon`)) %>%
  arrange(`Start codon`) %>%
  mutate(codons = factor(`Start codon`, levels = c("UUG", "GUG", "CUG", "AUG"))) %>%
  ggplot(aes(x=orf_count, y = `ORF type`, fill=codons)) +
  geom_bar(stat="identity") +
  facet_wrap(~facet_order, nrow = 2, 
             scales = "free_x", 
             labeller = labeller(facet_order= facet_order.labs)) +
  labs(x = "# Microproteins", y = "", fill = "") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_fill_viridis(discrete = T, direction = -1) +
  theme_bw() +
  theme(legend.position="top", 
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 10))  +
  guides(fill = guide_legend(reverse=TRUE))


ggsave(filename = paste0(results_dir,"Figure 6fg.png"), 
       device = "png", dpi = 5000, width = 4.5, height = 3.5)


## ----fig6h---------------------------------------------------------------------------------------------------------------------
d4d7_microproteins <- microproteins_pepquery %>%
  distinct(experiment, protein) %>%
  filter(experiment=="d4d7")

ts_microproteins <- microproteins_pepquery %>%
  distinct(experiment, protein) %>%
  filter(experiment=="tempshift")

ms_identifications <- list(
  `Temperature\nShift` = ts_microproteins$protein,
  `Growth\nPhase` = d4d7_microproteins$protein,
  `Drug\nProduct` = dp_microproteins$protein
)

ggvenn(
  ms_identifications,
  stroke_size = 0.5, 
  set_name_size = 3,
  show_percentage = F) +
  scale_fill_manual(values=c("#DDCC77","#117733", "#56B4E9"))

ggsave(filename = paste0(results_dir,"Figure 6h.png"), 
       device = "png", dpi = 2000, width = 3, height = 3)


## ----supp_data_8---------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 8.xlsx", sep = "")

suppressMessages(if (suppressMessages(file.exists(fn))) {
  file.remove(fn)})

canonical_lysate_identifications <- canonical_metamorpheus %>%
  dplyr::rename(Experiment=experiment) %>%
  dplyr::select(-study) %>%
  group_by(Experiment) %>%
  group_split()

orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")
orftable <- orftable %>%
  dplyr::select(`ORF-RATER name`, `ORF type`, `Associated Gene symbol`, `Associated Gene Name`, `Start codon`) %>%
  dplyr::rename(protein = `ORF-RATER name`) 


microproteins_lysate_identifications <- microproteins_pepquery %>%
  left_join(orftable, by="protein") %>%
  dplyr::rename(Experiment=experiment, Microprotein = protein) %>%
  dplyr::select(Experiment, Microprotein, `ORF type`, `Associated Gene symbol`, `Associated Gene Name`, `Start codon`,  everything()) %>%
  group_by(Experiment) %>%
  group_split()

write_xlsx(list(
  a = canonical_lysate_identifications[[2]] %>% dplyr::select(-c("class", "type")),
  b = microproteins_lysate_identifications[[2]],
  c = canonical_lysate_identifications[[1]] %>% dplyr::select(-c("class", "type")),
  d = microproteins_lysate_identifications[[1]]), 
  path = fn, format_headers = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
proDA_analysis <- readRDS("proteomics/differential_abundance/proda_results.rds")


## ------------------------------------------------------------------------------------------------------------------------------
paste0("_________________________________________________________________________________ ")
paste0("                         Quantified microproteins tempshift                       ")
paste0("__________________________________________________________________________________")
proDA_analysis$ntsts24hrs %>%
  filter(str_detect(name, "XR_|XM_|NR_|NM_")) %>%
  count()

paste0("_________________________________________________________________________________ ")
paste0("                         Quantified microproteins d4d7                      ")
paste0("__________________________________________________________________________________")
proDA_analysis$d4d7 %>%
  filter(str_detect(name, "XR_|XM_|NR_|NM_")) %>%
  count()



## ------------------------------------------------------------------------------------------------------------------------------
proDA_results <- bind_rows(
proDA_analysis$ntsts24hrs %>% mutate(comparison="ts24hrs"),
proDA_analysis$ntsts48hrs %>% mutate(comparison="ts48hrs"),
proDA_analysis$d4d7 %>% mutate(comparison="d4d7")) %>%
  dplyr::rename(protein=name) %>%
  left_join(orftable, by = "protein") %>%
  mutate(orf_type= case_when(
    is.na(`ORF type`) ~ "Canonical",
    !is.na(`ORF type`) ~ `ORF type`
  )) %>%
  mutate(class=case_when(
    abs(diff) >= log2(1.2) & adj_pval < 0.05 ~ orf_type,
    TRUE ~ "No change"))  %>%  
   arrange(class) %>%
  mutate(point_size = case_when(
    class != "No change" & orf_type == "Canonical" ~ 0.2,
    class != "No change" & orf_type != "Canonical" ~ 1,
    TRUE ~ 0.1
  )) %>%
  mutate(point_alpha = case_when(
     class != "No change" & orf_type == "Canonical" ~ 0.2,
    class != "No change" & orf_type != "Canonical"  ~ 0.8,
    TRUE ~ 0.1
  )) 


paste0("_________________________________________________________________________________ ")
paste0("                         differential RNA                                         ")
paste0("__________________________________________________________________________________")
table(proDA_results$comparison,proDA_results$class)


proDA_results$class <- factor(proDA_results$class, 
                                                    levels = c("No change", "Canonical", "Upstream", "Start overlap", "New"))

  proDA_results$facet_order = factor(proDA_results$comparison, levels=c('ts24hrs',"ts48hrs","d4d7"))   
  facet_order.labs <- c("24 h post-TS v NTS", "48 h post-TS v NTS", "Exponential v Stationary Phase")
  names(facet_order.labs) <- c("ts24hrs","ts48hrs","d4d7")
  
  proDA_results %>%
    ggplot(aes(x=diff, y=-log10(adj_pval), color=class)) + 
    geom_vline(xintercept=c(-log2(1.5),0, log2(1.5)), color = "light grey", linetype="dashed", linewidth=0.2) +
    geom_hline(yintercept=-log10(0.05), color = "light grey", linetype="dashed", linewidth=0.2) +
    geom_point(size=proDA_results$point_size, alpha=proDA_results$point_alpha) +
    scale_color_manual(values = c( "light grey", "#CC6677", "#035968", "#332288","#999933")) +
    theme_bw() +
    xlim(-8, 8) +
    facet_wrap(~facet_order, labeller = labeller(facet_order= facet_order.labs)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "Log"[2] ~ "Protein Abunandance FC", color = "") + 
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(strip.background = element_rect(fill = "transparent"),
          strip.text = element_text(size = 10),
          panel.spacing = unit(1, "cm", data = NULL)) 
  
  ggsave(filename = paste0(results_dir,"Figure 6hij.png"), 
         device = "png", dpi = 2000, width = 11, height = 3)


## ------------------------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 9.xlsx", sep = "")

proda_outputs <- proDA_results %>%
  dplyr::select(comparison, protein, `ORF type`, `Associated Gene symbol`, `Associated Gene Name`, pval, adj_pval, diff, t_statistic, se, df, avg_abundance, n_approx, n_obs,comparison, orf_type, class) %>%
  filter(class != "No change") %>%
  mutate(orf_type_split = case_when(
    orf_type == "Canonical" ~ "Canonical",
    orf_type != "Canonical" ~ "Microprotein"
  )) %>%
  group_by(comparison, orf_type_split) %>%
  group_split()

write_xlsx(list(
  a = proda_outputs[[3]] %>% dplyr::select(-c(`ORF type`, `Associated Gene symbol`, `Associated Gene Name`, orf_type, class, orf_type_split)),
  b = proda_outputs[[4]] %>% dplyr::select(-c( orf_type, class, orf_type_split)) %>% dplyr::rename(Microprotein=protein),
  c = proda_outputs[[5]] %>% dplyr::select(-c(`ORF type`, `Associated Gene symbol`, `Associated Gene Name`, orf_type, class, orf_type_split)),
  d = proda_outputs[[6]] %>% dplyr::select(-c( orf_type, class, orf_type_split)) %>% dplyr::rename(Microprotein=protein),
  e = proda_outputs[[1]] %>% dplyr::select(-c(`ORF type`, `Associated Gene symbol`, `Associated Gene Name`, orf_type, class, orf_type_split)), 
  f = proda_outputs[[2]] %>% dplyr::select(-c( orf_type, class, orf_type_split)) %>% dplyr::rename(Microprotein=protein)),
  path = fn, format_headers = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
print("section 2.6 results complete")

