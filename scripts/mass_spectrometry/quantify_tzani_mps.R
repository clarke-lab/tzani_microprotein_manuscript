#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("tidyverse", "Peptides","Biostrings")

suppressMessages(invisible(lapply(package_list, require, character.only = TRUE)))

# path variables
dp_flashlfq_path = "proteomics/flashlfq/result/drug_product/"
metamorpheus_path = "proteomics/metamorpheus/drug_product/tzani/"

# import the psms for microproteins found in drug product samples 
drug_product_psms <- readRDS("proteomics/protein_identifications/microproteins_drug_product.rds")

# determine the molecular weight of each microprotein
fastaFile <- readAAStringSet("proteomics/pepquery/peptide_digestion/microproteins.fasta")

protein = names(fastaFile)
protein_sequence = paste(fastaFile)

microprotein_info <- data.frame(protein, protein_sequence, mw=mw(protein_sequence)) %>%
  dplyr::rename(`Protein Accession` = protein,
                `Protein Unmodified Mass` = mw)



canonical_quantitation <- list()
microprotein_quantitation <-list()


products <- c("adalimumab","denosumab", "etanercept",
"nivolumab","pertuzumab", "vedolizumab")

for (product in products) {
  
  # import the flashlfq file for the product
  flashlfq_file=paste0(dp_flashlfq_path, product, "/QuantifiedProteins.tsv")
  
  flashlfq_data <- read_delim(flashlfq_file, delim = "\t", escape_double = FALSE, 
     trim_ws = TRUE, show_col = F) %>% 
  dplyr::rename(`Protein Accession`=`Protein Groups`) %>%
  select(-c(`Gene Name`, `Organism`,`...10`))

  # determine which proteins are quanified and calculate the average intensity
  flashlfq_data <- flashlfq_data %>%
    rowwise() %>%
    mutate(non_zero_count = sum(across(contains("Intensity_tzani_")) != 0)) %>%
    mutate(quantified = case_when(
      non_zero_count >= 4 ~ "Yes", 
      TRUE ~ "No"
    )) %>%
    filter(quantified != "No") %>%
    mutate(`Average Abundance` = rowMeans(across(contains("Intensity_tzani_")), na.rm =T))

  # determine which proteins were confidently detected by metamorpheus
  metamorpheus_file <- paste0(metamorpheus_path,product, "/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
  
  product_metamorpheus <- read_delim(metamorpheus_file, 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE, show_col = F) %>%
  filter(`Protein QValue` < 0.01 & as.numeric(`Number of Unique Peptides`) > 1) %>%
  filter(!str_detect(replace_na(`Protein Full Name`, ''), "Keratin"))  %>%
  select(`Protein Accession`, `Gene`, `Protein Full Name`, `Sequence Coverage Fraction`, `Number of Peptides`, 
  `Number of Unique Peptides`,  `Protein Unmodified Mass`, `Best Peptide Score`, `Protein QValue`) 
  
  # merge canonical with flashlfq data
  product_metamorpheus <- product_metamorpheus %>%
  left_join(flashlfq_data, by = "Protein Accession") 
  
  # store the abundance of the hi3 standard
  hi3_standard <- product_metamorpheus %>%
  filter(`Protein Accession` == "P63286") %>%
  mutate(fmol=150) %>%
  mutate(mol=fmol/1e15) %>%
  mutate(g_hcp=(mol*as.numeric(`Protein Unmodified Mass`))) %>%
  mutate(g_hcp_per_ug_protein=g_hcp/3) %>%
  mutate(g_hcp_per_mg_protein=g_hcp_per_ug_protein/0.001) %>%
  mutate(ppm=g_hcp_per_mg_protein*1000000000)

  
  # calculate the quantity of the canonical proteins using the standard
  product_metamorpheus <- product_metamorpheus %>%
  rowwise() %>%
  mutate(fmol=((hi3_standard$fmol/hi3_standard$`Average Abundance`)*`Average Abundance`)) %>%
  mutate(mol=fmol/1e15) %>%
  mutate(g_hcp=(mol*as.numeric(`Protein Unmodified Mass`))) %>%
  mutate(g_hcp_per_ug_protein=g_hcp/3) %>%
  mutate(g_hcp_per_mg_protein=g_hcp_per_ug_protein/0.001) %>%
  mutate(ppm=g_hcp_per_mg_protein*1000000000)

  # store the canonical proteins
  canonical_quantitation[[product]] <- product_metamorpheus

  # select the microproteins
  microprotein_psms <- drug_product_psms %>%
  filter(product == get("product", envir = .env) & study == "tzani")

  # add the molecular weigth
  microprotein_flash_lfq <- flashlfq_data %>% 
    filter(`Protein Accession` %in% unique(microprotein_psms$protein)) %>%
    left_join(microprotein_info, by= "Protein Accession") %>% 
    select(-c(protein_sequence,non_zero_count)) %>%
    select(`Protein Accession`,`Protein Unmodified Mass`, everything())
  
  microprotein_flash_lfq <- microprotein_flash_lfq %>%
    mutate(fmol=((hi3_standard$fmol/hi3_standard$`Average Abundance`)*`Average Abundance`)) %>%
    mutate(mol=fmol/1e15) %>%
    mutate(g_hcp=(mol*as.numeric(`Protein Unmodified Mass`))) %>%
    mutate(g_hcp_per_ug_protein=g_hcp/3) %>%
    mutate(g_hcp_per_mg_protein=g_hcp_per_ug_protein/0.001) %>%
    mutate(ppm=g_hcp_per_mg_protein*1e9)
  
  # calculate the quantity of the microproteins using the standard
  microprotein_quantitation[[product]] <- microprotein_flash_lfq

  print(cbind(microprotein_quantitation[[product]]$ppm,microprotein_quantitation[[product]]$`Protein Accession`))
}