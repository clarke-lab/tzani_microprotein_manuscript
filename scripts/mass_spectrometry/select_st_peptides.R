#!/usr/bin/env Rscript --vanilla

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

pep_q_val_thresh <- 0.10

metamorpheus_st_path <- "proteomics/metamorpheus_semi/"

data_types=c("drug_product","lysate" )

for (data_type in data_types) {

    if (data_type == "drug_product") {
        studies=c("tzani","pythoud")
    } else {
       studies=c("tzani") 
    }

    peptides_total <- data.frame()

    for (study in studies){

        if (study == "tzani" & data_type == "drug_product") {
            
            conditions <- "reducing"

            samples <- c("adalimumab", "denosumab", "etanercept", 
            "nivolumab", "pertuzumab", "vedolizumab")

        } else if (study == "tzani" & data_type == "lysate") {
            
            conditions <- "reducing"

            samples <- c("d4d7","tempshift")

        }  else if (study == "pythoud" & data_type == "drug_product") {
            
            conditions <- c("reducing", "native")

            samples <- c("adalimumab", "bevacizumab", "nivolumab", "trastuzumab")
        }

        for (condition in conditions) {

                for (sample in samples) {
                    
                    psm_file <- paste0(metamorpheus_st_path, "/", data_type, "/", study,"/", sample,"/",
                    condition, "/Task3SearchTask/AllPeptides.psmtsv")
                                        
                    peptides <- read_delim(psm_file, 
                                delim = "\t", escape_double = FALSE, show_col=F,
                                trim_ws = TRUE) %>%
                                filter(`PEP_QValue` < pep_q_val_thresh & `Peptide Description` == "Semi") %>%
                                select(`Base Sequence`, `Protein Accession`, `PEP_QValue`)

                    print(paste0(dim(peptides)[1], " peptides found in the ", study, " ", sample, " data" ))

                    peptides_total <- bind_rows(peptides, peptides_total)
                    
                

            }
        }
    }

    selected_peptides <- peptides_total %>% 
    distinct(`Base Sequence`, .keep_all=T)

    print(paste0(dim(selected_peptides)[1], " unique semi tryptic peptides found for the ", data_type, " data"))

    write.table(selected_peptides, 
                file = paste0("proteomics/pepquery/protein_digestion/metamorpheus_st_peptides_",data_type,".txt"), 
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}