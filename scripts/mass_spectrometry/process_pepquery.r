#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("tidyverse", "Spectra","MsBackendMgf","foreach", "stringi")

suppressMessages(invisible(lapply(package_list, require, character.only = TRUE)))

pepquery_path="proteomics/pepquery/"

source("scripts/mass_spectrometry/process_pepquery_functions.r")

# import peptide protein mapping
pep_to_prot_map <- read_delim(paste0(pepquery_path, 
                                 "peptide_digestion/microprotein_peptide_mapping.csv"), 
                          delim = ",", escape_double = FALSE,  trim_ws = TRUE, show_col = F) %>% 
  dplyr::rename(peptide=modified_sequence) %>% 
  filter(precursor_charge == 2) %>%
  select(peptide, protein)

# import drug product pepquery psms
drug_product_psms <- bind_rows(import_pepquery(pepquery_path, "drug_product", "all", pep_to_prot_map), 
import_pepquery(pepquery_path, "drug_product", "nterm", pep_to_prot_map))

# check if duplicate spectra are indentified as acetylated & non-acetylated
dup_spectra_check <- drug_product_psms %>%
 group_by(spectrum_title) %>%
 summarise(duplicate_count = sum(duplicated(spectrum_title))) %>%
  filter(duplicate_count > 0)
print(paste0(dim(dup_spectra_check)[1], " spectra duplicated in the drug product search"))


# identify high confident psms 
confident_tzani_psms <- drug_product_psms %>%
  filter(study=="tzani") %>% 
  distinct(product, rep, protein,.keep_all = T)  %>% 
  group_by(product, protein) %>%
  mutate(num_observations= n()) %>% 
  filter(num_observations > 2) # found in at least 3 replicates

paste0(dim(confident_tzani_psms)[1], " psms confidently identified in the tzani drug product data in at least 3 replicates per product")

paste0(length(unique(confident_tzani_psms$protein)), " microproteins found in the tzani drug product data")

pythoud_psms <- drug_product_psms  %>% 
  filter(study=="pythoud") %>% 
 distinct(product, sample_prep, rep, protein, .keep_all = T)  %>% 
  group_by(product, protein, sample_prep) %>%
 mutate(num_observations= n()) %>% arrange(-num_observations) %>%
  mutate(keep = case_when(
    str_detect(sample_prep, "sfp") & num_observations > 2 ~ "yes",  # semi fractionation (n =3)
    !str_detect(sample_prep, "sfp") & num_observations > 1 ~ "yes", # account for lower reps (n =3)
    TRUE ~ "no",
  ))

confident_pythoud_psms <- pythoud_psms %>%
  filter(keep=="yes")

paste0(dim(confident_pythoud_psms)[1], " psms confidently identified in the pythoud drug product data in at least 3 replicates per product")

paste0(length(unique(confident_pythoud_psms)), " microproteins found in the pythoud drug product data")

# find psms for microproteins across both studies
drug_product_microproteins <- unique(c(confident_tzani_psms$protein, confident_pythoud_psms$protein))

paste0(length(drug_product_microproteins), " microprotein found in the drug product data")

# select all psms found for the drug products. The microproteins were found separately,  
# but if we find evidence of the existence of the microprotein in the other study via
# a PSM in one replicate with consider this significant
drug_product_psms <- drug_product_psms %>%
  filter(protein %in% drug_product_microproteins)

# preparate for flashlfq

# import the mgf files to enable extraction of retention time
mgf_files <- list.files("proteomics/pepquery/mgf_files/drug_product/", pattern = "\\.mgf$", 
recursive = TRUE, full.names = TRUE)
drug_product_sps <- Spectra(mgf_files, source = MsBackendMgf())

# add the retention time and reformat the peptide modification string
drug_product_psms  %>%
  mutate(spectrum=gsub(":\\s*", "", str_extract(spectrum_title, ":\\s*([0-9]+)"))) %>%
  mutate(filename=paste0(str_extract(spectrum_title,"^[^:]+"),".mgf")) %>%
  rowwise() %>%
  mutate(modified_peptide=add_peptide_mods(modification, peptide)) %>%
 filter(str_detect(modification,  "Acet"))
  mutate(ret_time=add_retention_time(sps, spectrum, "drug_product", filename, mz)) 