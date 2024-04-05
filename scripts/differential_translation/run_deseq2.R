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

source("./scripts/differential_translation/deseq_utility_functions.R")

results_dir <- "differential_translation/deseq_output"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

############################################################
# 1. import annotation files
############################################################

# plastid
reference_annotation_plastid_merged <- read_delim(
  "differential_translation/plastid_reference/annotation_merged.txt",
  "\t", 
  escape_double = FALSE, 
  col_names = FALSE, 
  trim_ws = TRUE, 
  skip = 27, 
  show_col =F) %>% 
  mutate(X3 = str_remove(X1, "gene_")) 

# Chinese hamster  
CriGri_PICRH_1_0_annotation <- read_delim(
  "reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_feature_table.txt",
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE, 
  show_col =F)

# load the mouse annotation - used to match gene names to LOC ids
mouse_feature_table <- read_delim(
  "reference_genome/GCF_000001635.27_GRCm39_feature_table.txt",
  "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE, 
  show_col =F)

############################################################
# 2. Determine differential expression and translation
############################################################

# load the counts
count_data <- read.delim("differential_translation/combined_cds_counts.txt",
  sep = "\t", header = TRUE, row.names = "region")

# gene name causes a problem with the regex in fuzzy match, so we change it here manually
rownames(count_data)[which(rownames(count_data) == "LOC100760815")] <- "Maoa"

# experimental design
sample_table <- read.csv("data/differential_translation/experimental_design.txt",
  header = TRUE,  row.names = "sample_name")

sample_table$condition <- as.factor(sample_table$condition)
sample_table$assay <- as.factor(sample_table$assay)

# set the differential expression/translation criteria
fold_change_threshold <- 1.5
pval_threshold <- 0.05
base_mean_threshold <- 0
average_counts <- 10

# 1. differential gene expression
differential_expression <-  run_deseq(count_data, sample_table,
fold_change_threshold, pval_threshold, base_mean_threshold, average_counts,
"rnaseq", CriGri_PICRH_1_0_annotation,mouse_feature_table)

# 2. differential RPF 
differential_rpf <-  run_deseq(count_data, sample_table, fold_change_threshold, 
pval_threshold, base_mean_threshold, average_counts,"riboseq", CriGri_PICRH_1_0_annotation, mouse_feature_table)

# 3. differential translation efficency
differential_translation <-  run_deseq(count_data,sample_table, fold_change_threshold, pval_threshold, base_mean_threshold, 
average_counts,"translation", CriGri_PICRH_1_0_annotation,mouse_feature_table)




# save(differential_expression,differential_rpf,differential_translation, file = paste(results_dir, "/deseq_output.RData", sep = ""))