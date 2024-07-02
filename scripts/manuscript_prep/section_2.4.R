## ----setup---------------------------------------------------------------------------------------------------------------------
package_list <- c("tidyverse","viridis","gridExtra","ggpmisc", "writexl", "ggbreak",
                  "cowplot","UpSetR","ggupset", "patchwork")

lapply(package_list, require, character.only = TRUE)

# create result directory
results_dir <- "manuscript/section_2.4/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ----load_info-----------------------------------------------------------------------------------------------------------------
# orfrater
orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")

# metamorpheus
dp_canonical_proteins <- readRDS("proteomics/protein_identifications/canonical_proteins_drug_product.rds")

# pepquery
dp_microproteins <- readRDS("proteomics/protein_identifications/microproteins_drug_product.rds") %>%
  mutate(protein=str_remove(protein, "; XM_027410009.2_265668415_65aa"))


## ------------------------------------------------------------------------------------------------------------------------------
paste0("___________________________________________________________________________________")
paste0("            Number of microproteins identified across all drug products            ")
paste0("___________________________________________________________________________________")
dp_microproteins %>%
  distinct(protein, .keep_all = T) %>%
  count()

paste0("___________________________________________________________________________________")
paste0("            Number of microproteins with >=2 peptides identified                   ")
paste0("___________________________________________________________________________________")
pepides_detected <- dp_microproteins %>%
  distinct(protein, peptide, .keep_all = T) %>%
  group_by(protein) %>%
  count() %>%
  arrange(-n) 

print(paste0(sum(pepides_detected$n >= 2), " microproteins with 2 or more peptides identified"))


## ----count_hcps----------------------------------------------------------------------------------------------------------------
dp_canonical_counts <- dp_canonical_proteins %>%
  filter(!Organism %in% c("mAb", "fusion protein")) %>% 
  filter(!`Protein Accession` %in% c("clpB","Biognosys|iRT-Kit_WR_fusion", "P00924", "P00330", "P00489", "P02769")) %>%
  group_by(study, product) %>%
  summarise(count=n()) %>%
  mutate(class="canonical") 

dp_microprotein_counts <- dp_microproteins %>%
  distinct(study,product,protein, .keep_all = T) %>%
  group_by(study, product) %>%
  summarise(count=n()) %>%
  mutate(class="microprotein")

paste0("___________________________________________________________________________________")
paste0("                     Microprotein detections.                                      ")
paste0("___________________________________________________________________________________")
dp_microprotein_counts %>%
  dplyr::select(-class)


## ----fig4d---------------------------------------------------------------------------------------------------------------------
dp_hcp_counts <- bind_rows(dp_canonical_counts, dp_microprotein_counts)

bar_labels <- data.frame(product=str_to_title(dp_hcp_counts$product), 
                      label=dp_hcp_counts$count,  
                      yposition=dp_hcp_counts$count-1)

dp_hcp_counts$facet_order = factor(dp_hcp_counts$study, levels=c('tzani',"pythoud"))                      
dp_hcp_counts$facet_order2 = factor(dp_hcp_counts$class, levels=c('canonical',"microprotein"))    

facet.labs <- c("Tzani et. al", "Pythoud et al.")
names(facet.labs) <- c("tzani", "pythoud")

# barplot color palette
safe_colorblind_palette <- c("#CC6677", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#661100", "#6699CC", "#E69F00")

dp_hcp_counts %>%
  mutate(product=str_to_title(product)) %>%
  ggplot(aes(x=product, y=count, fill=product)) +
    geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_manual(values = safe_colorblind_palette) +
  
  labs(x = "", y="# HCPs detected", fill="") +
  geom_text(aes(x = bar_labels$product,
                y = bar_labels$yposition, 
                label = bar_labels$label),color = "white", size = 3) +
    facet_grid(facet_order2 ~ facet_order, scales = "free", space = "free", switch="y",
               labeller = labeller(facet_order = facet.labs)) +
  theme(legend.position="right", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 10))  +
  scale_y_continuous(expand = c(0, 0.6)) 

ggsave(filename = paste0(results_dir,"Figure 4d.png"), device = "png", dpi = 2000, width = 5, height = 4.5)


## ----supp_figure_8-------------------------------------------------------------------------------------------------------------
safe_colorblind_palette <- c("#CC6677", "#117733", "#44AA99",  "#6699CC")

dp_microprotein_counts_pythoud <- dp_microproteins %>%
  filter(study == "pythoud") %>%
  distinct(product,protein, .keep_all = T) %>%
  group_by(product, sample_prep) %>%
  summarise(count=n()) %>%
  mutate(class="microprotein") %>%
  mutate(sample_prep =case_when(
    sample_prep == "nd" ~ "Native digestion",
    sample_prep == "ond" ~ "Optimised native digestion",
    sample_prep == "old" ~ "Optimised liquid digestion",
    sample_prep == "sfp1" ~ "Semi-fractionation #1",
    sample_prep == "sfp2" ~ "Semi-fractionation #2",
    sample_prep == "sfp3" ~ "Semi-fractionation #3"
  ))

bar_labels <- data.frame(product=str_to_title(dp_microprotein_counts_pythoud$product), 
                      label=dp_microprotein_counts_pythoud$count,  
                      yposition=dp_microprotein_counts_pythoud$count-0.5)

dp_microprotein_counts_pythoud %>%
  mutate(product=str_to_title(product)) %>%
  ggplot(aes(x=product, y=count, fill=product)) +
    geom_bar(stat="identity") +
  facet_wrap(~ sample_prep) +
  scale_fill_manual(values = safe_colorblind_palette) +
  theme_bw() +
  theme(legend.position="right", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 10))  +
  labs(x = "", y="# microprotein HCPs detected", fill="Drug product") +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(x = bar_labels$product,
                y = bar_labels$yposition, 
                label = bar_labels$label),color = "white", size = 3) 

ggsave(filename = paste0(results_dir,"Supplementary Figure 8.png"), 
       device = "png", dpi = 2000, width = 7, height = 3.5)


## ------------------------------------------------------------------------------------------------------------------------------
paste0("___________________________________________________________________________________")
paste0("            Microproteins in 2 or more drug products                               ")
paste0("___________________________________________________________________________________")
num_observations <-dp_microproteins %>%
  distinct(protein, product) %>%
  group_by(protein) %>%
  summarise(count = n()) %>%
  arrange(-count)

print(paste0(sum(num_observations$count >= 2), " microproteins in 2 or more drug products"))


## ------------------------------------------------------------------------------------------------------------------------------
paste0("___________________________________________________________________________________")
paste0("            Shared across studies                                                  ")
paste0("___________________________________________________________________________________")
tzani_adalimumab <- (dp_microproteins %>% filter(str_detect(spectrum_title, "tzani_adalimumab")) %>%
                         dplyr::select(protein))$`protein`
tzani_nivolumab <- (dp_microproteins %>% filter(str_detect(spectrum_title, "tzani_nivolumab")) %>%
                         dplyr::select(protein))$`protein`
pythoud_adalimumab <- (dp_microproteins %>% filter(str_detect(spectrum_title, "pythoud_adalimumab")) %>%
                         dplyr::select(protein))$`protein`
pythoud_nivolumab <- (dp_microproteins %>% filter(str_detect(spectrum_title, "pythoud_nivolumab")) %>%
                         dplyr::select(protein))$`protein`

print("found in both adalimumab samples")
unique(intersect(tzani_adalimumab, pythoud_adalimumab))
print("found in both nivolumab samples")
unique(intersect(tzani_nivolumab, pythoud_nivolumab))


## ----fig4f---------------------------------------------------------------------------------------------------------------------
main_plot <- dp_microproteins %>%
  distinct(study, product, protein, .keep_all = T) %>%
  mutate(origin = paste(str_to_title(study), str_to_title(product), sep=" ")) %>%
  dplyr::select(origin, protein) %>%
  group_by(protein) %>%
  summarize(origins = list(origin)) %>%
  ggplot(aes(x = origins)) +
    geom_bar(size = 0.1, fill="black") +
    scale_x_upset() +
  theme_bw()+
  labs(x="", y="# microproteins") +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=2, color="white") 

upset_palette <- c("#E69F00","#CC6677","#332288", "#661100", "#44AA99", 
                            "#44AA99","#AA4499","#CC6677","#6699CC", "#117733")

side_plot <- dp_microproteins %>%
  distinct(study, product, protein, .keep_all = T) %>%
  mutate(origin = paste(str_to_title(study), str_to_title(product), sep=" ")) %>%
  dplyr::select(origin, protein) %>%
  unnest(cols = origin) %>%
  count(origin) %>%
  mutate(origin = fct_reorder(as.factor(origin), n)) %>%
  ggplot(aes(y = n, x = origin, fill=origin)) +
    geom_col() +
    coord_flip() +
  scale_fill_manual(values=upset_palette) +
    scale_y_reverse() +
    xlab("") + ylab("test") +
  theme_bw()+
  theme(plot.margin = unit(c(1, -5, -5, 1), "pt")) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position="none")

cowplot::plot_grid(
  cowplot::plot_grid(NULL, side_plot + theme(plot.margin = unit(c(1,-5, -5, 1), "pt")), ncol = 1, rel_heights = c(2.16, 2.16)),
  cowplot::plot_grid(NULL),
  main_plot + theme(plot.background = element_blank()), nrow = 1, rel_widths = c(1,0.6,3)
)

ggsave(filename = paste0(results_dir,"Figure 4f.png"), 
       device = "png", dpi = 2000, width = 8, height = 4)


## ------------------------------------------------------------------------------------------------------------------------------
proteins_identified <- dp_microproteins %>%
  distinct(protein, .keep_all = T) 

tid <- orftable %>%
  filter(`ORF-RATER name` %in% proteins_identified$protein) %>%
  select(`Transcript family`)

overlap_check <- orftable %>%
  filter(`Transcript family` %in% c(tid$`Transcript family`, "XM_027410009.2")) %>%
  filter(!`ORF type` == "Annotated") %>%
  rename(protein = `ORF-RATER name`) %>%
  left_join(proteins_identified, by="protein") %>%
  select(protein, `Associated Gene symbol`, `Transcript start position`, `Transcript stop position`, `Start codon`, peptide) %>%
  arrange(`Associated Gene symbol`,`Transcript start position`)

fn <- paste(results_dir, "overlap_check.xlsx", sep = "")
write_xlsx(overlap_check, path = fn, format_headers = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
quantdp_canonical <- readRDS("proteomics/flashlfq/result/drug_product/tzani_canonical_quant.rds")

quantplot_palette <- c("#CC6677", "#332288", "#AA4499","#E69F00")

canonical_quant <- enframe(quantdp_canonical) %>%
   unnest %>%
  arrange(-ppm) 

canonical_quant_selected <- canonical_quant %>%
  filter(!is.na(ppm)) %>%
  filter(!str_detect(`Protein Accession`, "_lc|_hc|etan|P63286")) %>%
  filter(`Number of Peptides` >= 3) %>%
  mutate(name = str_to_title(name)) %>%
  arrange(-ppm)

print("")
paste0("___________________________________________________________________________________")
paste0("                             Canonical HCP quantitation                            ")
paste0("___________________________________________________________________________________")

print("")
print("Canonical HCP ppm range")
round(range(canonical_quant_selected$ppm),2)
print("")
print(paste0("The median HCP concentration is: ", round(median(canonical_quant_selected$ppm),2)))



## ------------------------------------------------------------------------------------------------------------------------------
canonical_quant_selected %>%
   ggplot(aes(x=name, y=ppm, color=name, label=`Protein Accession`)) +
   geom_jitter(size=3,width = 0.15, alpha=0.75 ) + 
   theme_bw() +
   labs(x = "", y= "HCP abundance (ppm)", color = "Drug product", title="Canonical proteins") +
   # scale_color_manual(values = quantplot_palette) +
   guides(colour = guide_legend(override.aes = list(size=3)))


## ------------------------------------------------------------------------------------------------------------------------------
peptide_counts <- dp_microproteins %>%
  distinct(peptide, protein, product, study) %>%
  group_by(protein, product, study) %>%
  summarise(number_of_peptides=n()) %>%
  dplyr::rename(name=product, `Protein Accession`=protein) %>%
  arrange(-number_of_peptides)

quantdp <- readRDS("proteomics/flashlfq/result/drug_product/tzani_mp_quant.rds")


microprotein_quant <- enframe(quantdp) %>%
   unnest %>%
  left_join(peptide_counts, "name"="name", "Protein Accession"="Protein Accession") 

paste0("___________________________________________________________________________________")
paste0("             Abundance of microproteins with 3 peptides identified                  ")
paste0("___________________________________________________________________________________")
microprotein_quant %>% arrange(-ppm) %>% filter(number_of_peptides > 2) %>% select(name, ppm)

paste0("___________________________________________________________________________________")
paste0("             Abundance of microproteins with 2 peptides identified                  ")
paste0("___________________________________________________________________________________")
microprotein_quant %>% arrange(-ppm) %>% filter(number_of_peptides > 1 & number_of_peptides < 3) %>% select(name, ppm)

paste0("___________________________________________________________________________________")
paste0("             Abundance of microproteins with 1 peptide identified                  ")
paste0("___________________________________________________________________________________")
one_peptide_mp_quant <- microprotein_quant %>% arrange(-ppm) %>% filter(number_of_peptides == 1)  %>% select(name, ppm)

print(paste0(dim(one_peptide_mp_quant)[1], " microproteins quantified from 1 peptide"))

sum(one_peptide_mp_quant$ppm < median(canonical_quant_selected$ppm))

paste0("___________________________________________________________________________________")
paste0("        Microproteins with 1 peptide identified  above median canonical ppm        ")
paste0("___________________________________________________________________________________")
one_peptide_mp_quant[one_peptide_mp_quant$ppm > median(canonical_quant_selected$ppm),]


## ----supp_figure_9-------------------------------------------------------------------------------------------------------------
quantplot_palette <- c("#CC6677", "#332288", "#AA4499", 
                             "#44AA99", "#661100",  "#E69F00")
                             
microprotein_quant %>%
  mutate(name = str_to_title(name)) %>%
  ggplot(aes(x=name, y=ppm, color=name, label=`Protein Accession`, shape=as.factor(number_of_peptides))) +
  geom_jitter(size=3,width = 0.15, alpha=0.75, ) + 
     lims(y= c(0, 805)) +
  scale_y_break(breaks = c(8, 795)) +
  theme_bw() +
  labs(x = "", y= "HCP abundance (ppm)", color = "Drug product", 
       shape = "Microprotein\npeptides identified", title="Microproteins") +
  scale_color_manual(values = quantplot_palette) +
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(filename = paste0(results_dir,"Supplementary Figure 9.png"), device = "png", 
       dpi = 2000, width = 8, height = 5)


## ----supp_table_4--------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 4.xlsx", sep = "")

suppressMessages(if (suppressMessages(file.exists(fn))) {
  file.remove(fn)})

dp_canonical_list <- dp_canonical_proteins %>% 
  # filter(study == "tzani") %>%
  # full_join(canonical_quant, by=c("product"="name","Gene"="Gene")) %>%
  # dplyr::select(-c(contains("Intensity_tzani"),contains(".y"), "fmol","mol","g_hcp","g_hcp_per_ug_protein","g_hcp_per_mg_protein")) %>%
  dplyr::rename(Study=study, Product=product) %>%
  group_by(Study, Product) %>%
  group_split()

write_xlsx(list(
  a = dp_canonical_list[[5]],
  b = dp_canonical_list[[6]],
  c = dp_canonical_list[[7]],
  d = dp_canonical_list[[8]],
  e = dp_canonical_list[[9]],
  f = dp_canonical_list[[10]],
  g = dp_canonical_list[[1]],
  h = dp_canonical_list[[2]],
  i = dp_canonical_list[[3]],
  j = dp_canonical_list[[4]]
  ), path = fn, format_headers = TRUE)


## ----supp_table_5--------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 5.xlsx", sep = "")

suppressMessages(if (suppressMessages(file.exists(fn))) {
  file.remove(fn)})

orftable <- orftable %>%
  select(`ORF-RATER name`, `ORF type`, `Associated Gene symbol`, `Associated Gene Name`, `Start codon`) %>%
  rename(Microprotein=`ORF-RATER name`)

dp_microprotein_list <- dp_microproteins %>%
  dplyr::rename(Study=study, Product=product, Microprotein=protein) %>% 
  left_join(orftable, by = "Microprotein")   %>%
  group_by(Study, Product) %>%
  dplyr::select(Study, Product, Microprotein, `ORF type`, `Associated Gene symbol`, `Associated Gene Name`, `Start codon`,  everything()) %>%
  group_split()

write_xlsx(list(
  a = dp_microprotein_list[[5]],
  b = dp_microprotein_list[[6]],
  c = dp_microprotein_list[[7]],
  d = dp_microprotein_list[[8]],
  e = dp_microprotein_list[[9]],
  f = dp_microprotein_list[[10]],
  g = dp_microprotein_list[[1]],
  h = dp_microprotein_list[[2]],
  i = dp_microprotein_list[[3]],
  j = dp_microprotein_list[[4]]
  ), path = fn, format_headers = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
print("section 2.4 results complete")

