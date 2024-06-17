#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie


package_list <- c("tidyverse", "Biostrings")

suppressMessages(invisible(lapply(package_list, require, character.only = TRUE)))

pepquery_path="proteomics/pepquery/"


# make a folder for id results
digest_path <- "proteomics/pepquery/protein_digestion/"

# functions
read_fasta <- function(fasta_file) {
  fasta_sequences <- readAAStringSet(fasta_file)
  return(fasta_sequences)
}

find_nterminal_peptides <- function(fasta_sequences, peptides) {
  results <- list()
  
  for (name in names(fasta_sequences)) {
    protein_seq <- as.character(fasta_sequences[[name]])
    unique_peptides <- character(0)  # Initialize an empty character vector for unique peptides
    
    for (peptide in peptides) {
      if (startsWith(protein_seq, peptide) && !(peptide %in% unique_peptides)) {
        unique_peptides <- c(unique_peptides, peptide)
      }
    }
    
    if (length(unique_peptides) > 0) {
      results[[name]] <- unique_peptides
    }
  }
  
  return(results)
}

# write peptides from known proteins to a FASTA file
write_fasta <- function(df, output_file) {
  
  fasta_sequences <- AAStringSet(df$sequence)
  
  # header names
  names(fasta_sequences) <- paste(df$protein, df$count, sep = "_")
  
  writeXStringSet(fasta_sequences, filepath = output_file)
}

# load microprotein fasta
microprotein_sequences <- read_fasta(paste0(digest_path,"/microproteins.fasta"))

# load chainsaw digest
microproteins_digested_peptides <- read_table(paste0(digest_path,"/microproteins.fasta_digestedPeptides.tsv")) %>%
filter(missedCleavages <= 2)

print(paste0(dim(microproteins_digested_peptides)[1], " imported from chainsaw"))

# roll up peptides to proteins
microproteins_digested_peptides <- microproteins_digested_peptides %>%
  group_by(sequence) %>%
  arrange(protein, .by_group = T) %>%
  summarise(
    across(everything(),~paste0(unique(na.omit(.x)), collapse = "; ")))

print(paste0(dim(microproteins_digested_peptides)[1], " unique peptides from microproteins"))

# write list for pepquery
writeLines(microproteins_digested_peptides$sequence, 
paste0(digest_path,"/microprotein_all_peptides.txt"))

# identify which peptides are from the n-terminal
nterminal_peptides <- find_nterminal_peptides(microprotein_sequences, 
microproteins_digested_peptides$sequence) 

print(paste0(length(unique(unlist(nterminal_peptides))), " nterminal peptides"))

# write list for pepquery
writeLines(unique(unlist(nterminal_peptides)), 
paste0(digest_path,"/microprotein_nterm_peptides.txt"))

# write mapping table for unique peptides
write.table(microproteins_digested_peptides, sep="\t", row.names = F,
file=paste0(digest_path,"microprotein_peptide_mapping.tsv"))

known_fasta_digestedPeptides <- read_table(paste0(digest_path,"/known_proteins.fasta_digestedPeptides.tsv")) %>% 
mutate(origin="FT") %>%
filter(missedCleavages <= 2)

# semi-tryptic digestion for pepquery pass 1
# select the unique peptides from the known proteome for drug product data
st_peptides_dp <- read.table(paste0(digest_path,"/metamorpheus_st_peptides_drug_product.txt"), sep = "\t", header = FALSE) %>%
dplyr::rename(sequence=V1, protein=V2) %>% 
mutate(origin="ST") 

# in some cases metamorpheus places the same peptide 
# from multiple proteins in the same line, we keep only unique
extract_unique <- function(x) {
  unique_parts <- unique(str_split(x, "\\|")[[1]])
  unique_parts[1]
}

st_peptides_dp <- st_peptides_dp %>%
  mutate(sequence = sapply(sequence, extract_unique)) 

# merge full and semi
known_peptides_dp <- bind_rows(known_fasta_digestedPeptides, st_peptides_dp) %>%
  distinct(sequence, .keep_all = T) %>%
  group_by(protein) %>%
  mutate(count = row_number())

print(paste0(dim(known_peptides_dp)[1], " unique known peptides for dp search"))

write_fasta(known_peptides_dp,paste0(digest_path,"/known_peptides_drug_product.fasta"))

# select the unique peptides from the known proteome for the lysate data
st_peptides_lys <- read.table(paste0(digest_path,"/metamorpheus_st_peptides_lysate.txt"), sep = "\t", header = FALSE) %>%
dplyr::rename(sequence=V1, protein=V2) %>% 
mutate(origin="ST")


# in some cases metamorpheus places the same peptide 
# from multiple proteins in the same line, we keep only unique
st_peptides_lys <- st_peptides_lys %>%
  mutate(sequence = sapply(sequence, extract_unique)) 

known_peptides_lys <- bind_rows(known_fasta_digestedPeptides, st_peptides_lys) %>%
  distinct(sequence, .keep_all = T) %>%
  group_by(protein) %>%
  mutate(count = row_number()) 

print(paste0(dim(known_peptides_lys)[1], " unique known peptides for lysate search"))

write_fasta(known_peptides_lys,paste0(digest_path,"/known_peptides_lysate.fasta"))