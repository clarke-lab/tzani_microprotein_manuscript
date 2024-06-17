#!/usr/bin/env Rscript --vanilla
#### Description: perform differential abundance analysis
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

package_list <- c("tidyverse", "proDA")

suppressMessages(invisible(lapply(package_list, require, character.only = TRUE)))

# path variables
lysate_flashlfq_path = "proteomics/flashlfq/result/lysate/"
metamorpheus_path = "proteomics/metamorpheus/lysate/tzani/"

experiments <- c("d4d7", "tempshift")

for (experiment in experiments) {

    # determine which proteins were confidently detected by metamorpheus
    metamorpheus_file <- paste0(metamorpheus_path, experiment, "/reducing/Task3SearchTask/AllQuantifiedProteinGroups.tsv")
  
    lysate_metamorpheus <- read_delim(metamorpheus_file, 
        delim = "\t", escape_double = FALSE, 
        trim_ws = TRUE, show_col = F) %>%
    filter(`Protein QValue` < 0.01 & as.numeric(`Number of Unique Peptides`) > 1) %>%
    filter(!str_detect(replace_na(`Protein Full Name`, ''), "Keratin"))  %>%
    select(`Protein Accession`, `Gene`, `Protein Full Name`, `Sequence Coverage Fraction`, `Number of Peptides`, 
    `Number of Unique Peptides`,  `Protein Unmodified Mass`, `Best Peptide Score`, `Protein QValue`) 

    # import the flashlfq results
    lysate_quantified_proteins <- read_delim(paste0(lysate_flashlfq_path,experiment,"/QuantifiedProteins.tsv"), 
        delim = "\t", escape_double = FALSE, 
        col_types = cols(`Gene Name` = col_skip(), 
            Organism = col_skip()), 
        trim_ws = TRUE  , show_col = F) %>%
        dplyr::select(`Protein Groups`, contains("Intensity_"))

    # cross reference flashlfq with significant canonical proteins and 
    # remove quantified proteins which were not Q < 0.01 at the protein levle
    lysate_quantified_proteins <- lysate_quantified_proteins %>%
    mutate(decision = case_when(
        str_detect(`Protein Groups`, "XR_|XM_|NM_|NR_") ~ "keep",
        (`Protein Groups` %in% lysate_metamorpheus$`Protein Accession`) ~ "keep",
        TRUE ~ "remove")) %>%
        filter(decision == "keep") %>%
        dplyr::select(-decision)
    
    if (experiment == "tempshift") {
      
      lysate_quantified_proteins <- lysate_quantified_proteins %>%
        rowwise() %>%
        mutate(
          non_zero_nts = sum(c_across(starts_with("Intensity_nts_24hr")) != 0),
          non_zero_ts_24 = sum(c_across(starts_with("Intensity_ts_24hr")) != 0),
          non_zero_ts_48 = sum(c_across(starts_with("Intensity_ts_48hr")) != 0),
        ) %>%
        # filter(str_detect(`Protein Groups`, 'XR_|XM_|NR_|XR_')) %>%
        filter(!non_zero_nts < 3 & !non_zero_ts_24 < 3 & !non_zero_ts_48 < 3) %>%
        dplyr::select(-c(non_zero_nts, non_zero_ts_24, non_zero_ts_48))
        
        print("quantified proteins in tempshift samples")
        quant_count <- lysate_quantified_proteins %>%
        mutate(class = case_when(
        str_detect(`Protein Groups`, 'XR_|XM_|NR_|XR_') ~ "microprotein",
        TRUE ~ "canonical"
        )) %>%
        group_by(class) %>%
        count()
        
        print(quant_count)
      
    } else {
      
      lysate_quantified_proteins <- lysate_quantified_proteins %>%
        rowwise() %>%
        mutate(
          non_zero_d4 = sum(c_across(starts_with("Intensity_d4")) != 0),
          non_zero_d7 = sum(c_across(starts_with("Intensity_d7")) != 0)
        ) %>%
        # filter(str_detect(`Protein Groups`, 'XR_|XM_|NR_|XR_')) %>%
        filter(!non_zero_d4 < 3 & !non_zero_d7 < 3) %>%
        dplyr::select(-c(non_zero_d4, non_zero_d7))

        # determine the number of canonical and microproteins quantitated
        print("quantified proteins in d4d7 samples")
        quant_count <- lysate_quantified_proteins %>%
        mutate(class = case_when(
        str_detect(`Protein Groups`, 'XR_|XM_|NR_|XR_') ~ "microprotein",
        TRUE ~ "canonical"
        )) %>%
        group_by(class) %>%
        count()

        print(quant_count)
    }
    
    # make the abundance matix for proDA
    abundance_matrix <- as.matrix(lysate_quantified_proteins[,-1])
    rownames(abundance_matrix) <- lysate_quantified_proteins$`Protein Groups`

    # set missing values to NA as required 
    abundance_matrix[abundance_matrix == 0] <- NA
    abundance_matrix <- log2(abundance_matrix)
    normalized_abundance_matrix <- median_normalization(abundance_matrix)

    # order the column names 
    normalized_abundance_matrix  <- normalized_abundance_matrix[ , order(colnames(normalized_abundance_matrix))]

    # create sample info with condition replicate
    sample_info_df <- data.frame(
    name = colnames(normalized_abundance_matrix),
    stringsAsFactors = FALSE)
    
    # create sample info based on experiment
    if (experiment == "tempshift") {
        
        sample_info_df$condition <- c("nts_24hr", "nts_24hr", "nts_24hr", "nts_24hr", "nts_24hr", "nts_24hr", 
        "ts_24hr", "ts_24hr", "ts_24hr", "ts_24hr", "ts_24hr", "ts_24hr", 
        "ts_48hr", "ts_48hr", "ts_48hr", "ts_48hr", "ts_48hr", "ts_48hr") 
        sample_info_df$replicate <- c(1, 1, 2, 2, 3, 3,1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3)
        comparator_level = "nts_24hr"

        print(sample_info_df)

    } else {
        
        sample_info_df$condition <- c("d4", "d4", "d4", "d4", "d7", "d7", "d7", "d7")
        sample_info_df$replicate <- c(1, 2, 3, 4, 1, 2, 3, 4)
        comparator_level = "d4"

        print(sample_info_df)

    }
    
    # run proDA
    fit <- proDA(normalized_abundance_matrix,
    design = ~condition,
    col_data = sample_info_df, 
    reference_level = comparator_level)

    # select results based on current experiment
    if (experiment == "tempshift") {
        
        ntsts_24hrs <- test_diff(fit, "conditionts_24hr") 

        ntsts_48hrs <- test_diff(fit, "conditionts_48hr")

    } else {
        result_d4d7 <- test_diff(fit, "conditiond7") 
    }
}

# save results
proda_results <- list()
proda_results[["ntsts24hrs"]] <- ntsts_24hrs
proda_results[["ntsts48hrs"]] <- ntsts_48hrs
proda_results[["d4d7"]] <- result_d4d7

# make a folder for da results
diff_abund_dir <- "proteomics/differential_abundance/"
if (!dir.exists(diff_abund_dir)) {
  dir.create(diff_abund_dir, recursive = TRUE)
}

saveRDS(proda_results, file = paste0(diff_abund_dir,"proda_results.rds"))