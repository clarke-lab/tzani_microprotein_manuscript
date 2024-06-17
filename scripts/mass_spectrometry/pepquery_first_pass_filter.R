#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("tidyverse")
suppressMessages(invisible(lapply(package_list, require, character.only = TRUE)))


args = commandArgs(trailingOnly=TRUE)

pepquery_path="proteomics/pepquery"

# make a folder for id results
psm_rank <- read_delim(paste0(pepquery_path,"/results/", args[1], "/", args[2],"/first_pass/", args[3], "/psm_rank.txt"), 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)

print(dim(psm_rank))

# select peptides based on PepQuery statistical evaluation and no unrestricted modifications
 first_pass_peptides <- psm_rank %>%
 mutate(keep = case_when(
     str_count(peptide) <= 8 & pvalue < 0.05 & rank == 1~ "Yes",
     str_count(peptide) > 8 & pvalue < 0.01 & rank == 1 ~ "Yes", 
     TRUE ~ "No"
 )) %>%
 filter(keep == "Yes") %>%
 distinct(peptide)

print(length(first_pass_peptides$peptide))

fileConn<-file(paste0(pepquery_path,"/results/", args[1], "/",args[2],"/first_pass/", args[3],"/selected_peptides.txt"))
 writeLines(first_pass_peptides$peptide, fileConn)
 close(fileConn)
