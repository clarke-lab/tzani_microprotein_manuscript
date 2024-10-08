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

source("scripts/mass_spectrometry/process_pepquery_functions.R")

# make a folder for id results
protein_id_dir <- "proteomics/protein_identifications/"
if (!dir.exists(protein_id_dir)) {
  dir.create(protein_id_dir, recursive = TRUE)
}

# import peptide protein mapping
pep_to_prot_map <- read_delim("proteomics/pepquery/protein_digestion/microprotein_peptide_mapping.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE) %>%
    dplyr::rename(peptide=sequence)

# first pass 
drug_product_psms <- bind_rows( 
import_pepquery(pepquery_path, "drug_product","all","reducing",pep_to_prot_map),
import_pepquery(pepquery_path, "drug_product","nterm","reducing",pep_to_prot_map),
import_pepquery(pepquery_path, "drug_product","all","native",pep_to_prot_map),
import_pepquery(pepquery_path, "drug_product","nterm","native",pep_to_prot_map)
)

# check if duplicate spectra are indentified as acetylated & non-acetylated
dup_spectra_check <- drug_product_psms %>%
 group_by(spectrum_title) %>%
 summarise(duplicate_count = sum(duplicated(spectrum_title))) %>%
  filter(duplicate_count > 0)

#duplicated_spectra <- drug_product_psms_tryptic %>%
#filter(spectrum_title %in% dup_spectra_check$spectrum_title)

#write.csv(duplicated_spectra, file="duplicated_spectra.csv")

print(paste0(dim(dup_spectra_check)[1], " spectra duplicated in the drug product search"))


# identify high confident psms 
confident_tzani_psms <- drug_product_psms %>%
  filter(study=="tzani") %>% 
  distinct(product, rep, protein,.keep_all = T)  %>% 
  group_by(product,protein ) %>%
  mutate(num_observations= n()) %>% 
  filter(num_observations > 2) # found in at least 3 replicates

paste0(dim(confident_tzani_psms)[1], " psms confidently identified in the tzani drug product data in at least 3 replicates per product")

paste0(length(unique(confident_tzani_psms$protein)), " microproteins found in the tzani drug product data")


pythoud_psms <- drug_product_psms  %>% 
  filter(study=="pythoud") %>% 
 distinct(product, sample_prep, rep,protein, .keep_all = T)  %>% 
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

paste0(length(unique(confident_pythoud_psms$protein)), " microproteins found in the pythoud drug product data")

# find psms for microproteins across both studies
drug_product_microproteins <- unique(c(confident_tzani_psms$protein, confident_pythoud_psms$protein))

paste0(length(drug_product_microproteins), " microprotein found in the drug product data")

# select all psms found for the drug products. The microproteins were found separately,  
# but if we find evidence of the existence of the microprotein in the other study via
# a PSM in one replicate with consider this significant
drug_product_psms<- drug_product_psms %>%
  filter(protein %in% drug_product_microproteins)


saveRDS(drug_product_psms, file = paste0(protein_id_dir,"microproteins_drug_product.rds"))

# import the canonical proteins and psms
studies <- c("tzani", "pythoud")

canonical_proteins_drug_product <- data.frame()
canonical_psms_drug_product <- data.frame()
type="drug_product"
for (study in studies) {
  if (study == "tzani") {
    samples= c("adalimumab", "denosumab", "pertuzumab","nivolumab", 
               "vedolizumab", "etanercept")
  } else {
    samples = c("adalimumab", "bevacizumab" ,"nivolumab" ,"trastuzumab")
  }
  
  for (sample in samples) {

    if (study=="pythoud") {
      imported_canonical_native <- import_canonical(type,study, sample, "native") 

      imported_canonical_reducing <- import_canonical(type,study, sample, "reducing") 

      canonical_proteins_drug_product <- bind_rows(canonical_proteins_drug_product, 
      imported_canonical_native$proteins %>% mutate(prep="native"),
      imported_canonical_reducing$proteins %>% mutate(prep="reducing"))

      
    
    } else {
      imported_canonical <- import_canonical(type,study, sample, "reducing") 

      canonical_proteins_drug_product <- bind_rows(canonical_proteins_drug_product, 
      imported_canonical$proteins %>% mutate(prep="reducing"))

      canonical_psms_drug_product <- bind_rows(canonical_psms_drug_product, imported_canonical$psms)
    }
  }
}

# output file containing protein identifications
saveRDS(canonical_proteins_drug_product, file = paste0(protein_id_dir,"canonical_proteins_drug_product.rds"))
  
# preparate for flashlfq

flashlfq_dir <- "proteomics/flashlfq/psm/drug_product"
if (!dir.exists(flashlfq_dir)) {
  dir.create(flashlfq_dir, recursive = TRUE)
}

# import the mgf files to enable extraction of retention time
mgf_files <- list.files("proteomics/pepquery/mgf_files/drug_product", pattern = "\\.mgf$", 
recursive = TRUE, full.names = TRUE)
drug_product_sps <- Spectra(mgf_files, source = MsBackendMgf())

# add the retention time and reformat the peptide modification string
dp_psms_mp_flfq <- drug_product_psms  %>%
filter(study == "tzani") %>%
  mutate(spectrum=gsub(":\\s*", "", str_extract(spectrum_title, ":\\s*([0-9]+)"))) %>%
  mutate(filename=str_extract(spectrum_title,"^[^:]+")) %>%
  rowwise() %>%
  mutate(modified_peptide=add_peptide_mods(modification, peptide)) %>%
  mutate(ret_time=add_retention_time(drug_product_sps, spectrum, "drug_product", paste0(filename, ".mgf"), mz)) %>%
  dplyr::select(study, product, filename,peptide, modified_peptide, pep_mass, ret_time, charge, protein) %>%
  dplyr::rename(`File Name` = filename, `Base Sequence`= peptide, `Full Sequence`=modified_peptide,
  `Peptide Monoisotopic Mass`= pep_mass, `Scan Retention Time` = ret_time, `Precursor Charge`=charge,
  `Protein Accession` = protein)

dp_psms_mp_flfq <- dp_psms_mp_flfq %>% 
group_by(study, product) 
  
dp_psms_mp_flfq_names <- group_keys(dp_psms_mp_flfq) %>%
    mutate(group_name = str_c(as.character(study),"_",product))
  
dp_psms_mp_flfq_names <- dp_psms_mp_flfq_names$group_name

dp_psms_mp_flfq_list <- group_split(dp_psms_mp_flfq) %>%
    setNames(dp_psms_mp_flfq_names)

# make a list for the canonical hcps 

dp_psms_canonical_flfq <- canonical_psms_drug_product %>% 
group_by(study, product) 

dp_psms_canonical_flfq_names <- group_keys(dp_psms_canonical_flfq) %>%
mutate(group_name = str_c(as.character(study),"_",product))

dp_psms_canonical_flfq_names <- dp_psms_canonical_flfq_names$group_name

dp_psms_canonical_flfq_list <- group_split(dp_psms_canonical_flfq) %>%
    setNames(dp_psms_canonical_flfq_names)

for (set_name in dp_psms_mp_flfq_names) {
    
    current_psm_set <- bind_rows(dp_psms_mp_flfq_list[[set_name]], dp_psms_canonical_flfq_list[[set_name]]) %>%
    mutate(`File Name`=paste0("proteomics/metamorpheus/drug_product/","/", study, "/", product ,"/Task2CalibrationTask/",`File Name`)) %>%
    dplyr::select(-c(study, product))

    write_tsv(current_psm_set, file = paste0(flashlfq_dir, "/", set_name, ".psms.tsv" ))
}

# import the lysate data 
# import peptide protein mapping

# import drug product pepquery psms
lysate_psms <- bind_rows(import_pepquery(pepquery_path, "lysate", "all", "reducing",pep_to_prot_map), 
import_pepquery(pepquery_path, "lysate", "nterm","reducing", pep_to_prot_map))

# check if duplicate spectra are indentified as acetylated & non-acetylated
dup_spectra_check <- lysate_psms %>%
 group_by(spectrum_title) %>%
 summarise(duplicate_count = sum(duplicated(spectrum_title))) %>%
  filter(duplicate_count > 0)
print(paste0(dim(dup_spectra_check)[1], " spectra duplicated in the lysate search"))

# identify high confident psms for the tempshift proteomics
confident_tempshift_psms <- lysate_psms %>%
  filter(experiment=="tempshift") %>% 
  distinct(condition, rep, protein,.keep_all = T)  %>% 
  group_by(condition, protein) %>%
  mutate(num_observations= n()) %>% 
  filter(num_observations > 1) # found in at least 2 replicates of a condition

paste0(dim(confident_tempshift_psms)[1], " psms confidently identified in the lysate tempshift data in at least 2 replicates per condition")

paste0(length(unique(confident_tempshift_psms$protein)), " microproteins found in the lysate tempshift data")

# identify high confident psms for the d4d7 proteomics
confident_d4d7_psms <- lysate_psms %>%
  filter(experiment=="d4d7") %>% 
  distinct(condition, rep, protein,.keep_all = T)  %>% 
  group_by(condition, protein) %>%
  mutate(num_observations= n()) %>% 
  filter(num_observations > 1) # found in at least 2 replicates of a condition

paste0(dim(confident_d4d7_psms)[1], " psms confidently identified in the lysate d4d7 data in at least 2 replicates per condition")

paste0(length(unique(confident_d4d7_psms$protein)), " microproteins found in the lysate d4d7 data")

# find psms for microproteins across both experiments
lysate_microproteins <- unique(c(confident_tempshift_psms$protein, confident_d4d7_psms$protein))

paste0(length(lysate_microproteins), " microprotein found in the lysate data")

# select all psms found for the lysate experiments. The microproteins were found separately,  
# but if we find evidence of the existence of the microprotein in the other experiment via
# a PSM in one replicate with consider this significant
lysate_psms <- lysate_psms %>%
  filter(protein %in% lysate_microproteins)

saveRDS(lysate_psms, file = paste0(protein_id_dir,"microproteins_lysate.rds"))

# import the lysate protein and psms
tempshift_canonical <- import_canonical("lysate","tzani", "tempshift", "reducing")
d4d7_canonical <- import_canonical("lysate","tzani", "d4d7", "reducing")

canonical_proteins_lysate <- bind_rows(tempshift_canonical$protein, d4d7_canonical$protein)
canonical_psms_lysate <- bind_rows(tempshift_canonical$psm, d4d7_canonical$psm)

# output file containing protein identifications
saveRDS(canonical_proteins_lysate, file = paste0(protein_id_dir,"canonical_proteins_lysate.rds"))
imported_canonical <- import_canonical("lysate","tzani", "tempshift", "reducing")

# preparate for flashlfq

flashlfq_dir <- "proteomics/flashlfq/psm/lysate"
if (!dir.exists(flashlfq_dir)) {
  dir.create(flashlfq_dir, recursive = TRUE)
}

# import the mgf files to enable extraction of retention time
mgf_files <- list.files("proteomics/pepquery/mgf_files/lysate/", pattern = "\\.mgf$", 
recursive = TRUE, full.names = TRUE)
lysate_sps <- Spectra(mgf_files, source = MsBackendMgf())

# add the retention time and reformat the peptide modification string
lysate_psms_mp_flfq <- lysate_psms  %>%
  mutate(spectrum=gsub(":\\s*", "", str_extract(spectrum_title, ":\\s*([0-9]+)"))) %>%
  mutate(filename=str_extract(spectrum_title,"^[^:]+")) %>%
  rowwise() %>%
  mutate(modified_peptide=add_peptide_mods(modification, peptide)) %>%
  mutate(ret_time=add_retention_time(lysate_sps, spectrum, "lysate", paste0(filename, ".mgf"), mz)) %>%
  dplyr::select(experiment, filename, peptide, modified_peptide, pep_mass, ret_time, charge, protein) %>%
  dplyr::rename(`File Name` = filename, `Base Sequence`= peptide, `Full Sequence`=modified_peptide,
  `Peptide Monoisotopic Mass`= pep_mass, `Scan Retention Time` = ret_time, `Precursor Charge`=charge,
  `Protein Accession` = protein)

lys_psms_mp_flfq <- lysate_psms_mp_flfq %>% 
group_by(experiment) 
  
lys_psms_mp_flfq_names <- group_keys(lys_psms_mp_flfq) %>%
    mutate(group_name = experiment)
  
lys_psms_mp_flfq_names <- lys_psms_mp_flfq_names$group_name

lys_psms_mp_flfq_list <- group_split(lys_psms_mp_flfq) %>%
    setNames(lys_psms_mp_flfq_names)

# make a list for the canonical hcps 

lys_psms_canonical_flfq <- canonical_psms_lysate %>%
select(-study) %>%
group_by(experiment) 

lys_psms_canonical_flfq_names <- group_keys(lys_psms_canonical_flfq) %>%
mutate(group_name = experiment)

lys_psms_canonical_flfq_names <- lys_psms_canonical_flfq_names$group_name

lys_psms_canonical_flfq_list <- group_split(lys_psms_canonical_flfq) %>%
    setNames(lys_psms_canonical_flfq_names)

for (set_name in lys_psms_mp_flfq_names) {
    
    current_psm_set <- bind_rows(lys_psms_mp_flfq_list[[set_name]], lys_psms_canonical_flfq_list[[set_name]]) %>%
    mutate(`File Name`=paste0("proteomics/metamorpheus/lysate/tzani/", experiment,"/Task2CalibrationTask/",`File Name`)) %>%
    dplyr::select(-c(experiment))

    write_tsv(current_psm_set, file = paste0(flashlfq_dir, "/tzani_", set_name, ".psms.tsv" ))
}
