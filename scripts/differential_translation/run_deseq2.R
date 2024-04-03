#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentiall translated genes  
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

## 1. Prepare for analysis
package_list <- c(
  "tidyverse", "DESeq2", "fuzzyjoin", "readr", "data.table", 
  "ggrepel"
)

invisible(lapply(package_list, require, character.only = TRUE))

source("scripts/orf_filtering_functions.R")

results_dir <- "differential_translation/deseq_output"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# load the counts
count_data <- read.delim("differential_translation/gene_cds_counts.txt",
  sep = "\t", header = TRUE, row.names = "region")
rownames(count_data)[which(rownames(count_data) == "LOC100760815_1")] <- "Maoa"

sample_table <- read.csv("data/experimental_design/sequencing_design.txt",
  header = TRUE,  row.names = "sample_name")

# import the annotation
reference_annotation_plastid_merged <- read_delim("differential_translation/plastid_reference/annotation_merged.txt",
  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 27
) %>% mutate(X3 = str_remove(X1, "gene_")) %>%
  mutate(`product_accession` = str_remove(X3, "_1"))

# import the ncbi annotation
CriGri_PICRH_1_0_annotation <- read_delim( "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt",
  "\t", escape_double = FALSE, trim_ws = TRUE)

# load the mouse annotation - used to match gene names to complete LOC to symbol conversion
mouse_feature_table <- read_delim("reference_genome/GCF_000001635.27_GRCm39_feature_table.txt",
  "\t",   escape_double = FALSE, trim_ws = TRUE)

# Differential expression and translation
# set the differential expression/translation criteria
fold_change_threshold <- 1.5
pval_threshold <- 0.05
base_mean_threshold <- 0
average_counts <- 20

# 2. RNA-seq only differential expression

differential_expression <-  run_deseq(count_data,
                                       sample_table,
                                       fold_change_threshold, 
                                       pval_threshold, 
                                       base_mean_threshold, 
                                       average_counts,
                                       "rnaseq",
                                       CriGri_PICRH_1_0_annotation,
                                       mouse_feature_table)

differential_rpf <-  run_deseq(count_data,
                                      sample_table,
                                      fold_change_threshold, 
                                      pval_threshold, 
                                      base_mean_threshold, 
                                      average_counts,
                                      "riboseq",
                                      CriGri_PICRH_1_0_annotation,
                                      mouse_feature_table)

differential_translation <-  run_deseq(count_data,sample_table, fold_change_threshold, 
          pval_threshold, 
          base_mean_threshold, 
          average_counts,
          "translation",
          CriGri_PICRH_1_0_annotation,
          mouse_feature_table)

# save(differential_expression,differential_rpf,differential_translation, file = paste(results_dir, "/deseq_output.RData", sep = ""))