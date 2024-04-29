#!/usr/bin/env Rscript --vanilla
#### Description: Use DESeq2 to finds
####              1. DE genes from the RNA-seq data
####              2. DE genes from the Ribo-seq data
####              3. Differentially translated genes 
####              4. Classify events
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
rna_deseq <-  run_deseq(count_data, sample_table,
fold_change_threshold, pval_threshold, base_mean_threshold, average_counts,
"rnaseq", CriGri_PICRH_1_0_annotation,mouse_feature_table)

# 2. differential RPF 
rpf_deseq <-  run_deseq(count_data, sample_table, fold_change_threshold, 
pval_threshold, base_mean_threshold, average_counts,"riboseq", CriGri_PICRH_1_0_annotation, mouse_feature_table)

# 3. differential translation efficency
te_deseq <-  run_deseq(count_data,sample_table, fold_change_threshold, pval_threshold, base_mean_threshold, 
average_counts,"translation", CriGri_PICRH_1_0_annotation,mouse_feature_table)

# write the scale factors to file to create coverage tracks

scale_factors <- sizeFactors(te_deseq$deseq_object) %>% 
as_tibble(rownames="sample") %>%
 mutate(value = (value^-1))

scale_factors_rnaseq <-scale_factors %>% filter(str_detect(sample, "rnaseq"))
write_tsv(scale_factors_rnaseq, file = paste(results_dir, "joint_rna_scale_factors.txt", sep = "/"), col_names =F)

scale_factors_riboseq <-scale_factors %>% filter(str_detect(sample, "riboseq"))
write_tsv(scale_factors_riboseq, file = paste(results_dir, "joint_rpf_scale_factors.txt", sep = "/"), col_names =F)

# save the deseq objects, res, and significant
saveRDS(rna_deseq, file = paste(results_dir, "/rna_deseq.rds", sep = ""))
saveRDS(rpf_deseq, file = paste(results_dir, "/rpf_deseq.rds", sep = ""))
saveRDS(te_deseq, file = paste(results_dir, "/te_deseq.rds", sep = ""))

############################################################
# 2. Delta TE
############################################################

# classify regulation
res_te <- te_deseq$deseq_res
res_rna <- rna_deseq$deseq_res[rownames(res_te),]
res_ribo <- rpf_deseq$deseq_res[rownames(res_te),]

fold_change_threshold <- 0
transcriptional_select=which(res_te$padj >= 0.05 & 
res_ribo$padj < 0.05 & abs(res_ribo$log2FoldChange) >= log2(fold_change_threshold) & 
res_rna$padj < 0.05 & abs(res_rna$log2FoldChange) >= log2(fold_change_threshold))
transcriptional = rownames(res_te)[transcriptional_select]

translation_exclusive_select=which(res_te$padj < 0.05 & abs(res_te$log2FoldChange) >= log2(fold_change_threshold) & 
res_ribo$padj < 0.05 & abs(res_ribo$log2FoldChange) >= log2(fold_change_threshold) &
res_rna$padj >= 0.05)
translation_exclusive = rownames(res_te)[translation_exclusive_select]

both = which(res_te$padj < 0.05 & abs(res_te$log2FoldChange) >= log2(fold_change_threshold) & 
res_ribo$padj < 0.05 & abs(res_ribo$log2FoldChange) >= log2(fold_change_threshold) &
res_rna$padj < 0.05 & abs(res_rna$log2FoldChange) >= log2(fold_change_threshold))

intensified = rownames(res_te)[both[which(res_te[both,2]*res_rna[both,2] > 0)]]
buffered = rownames(res_te)[both[which(res_te[both,2]*res_rna[both,2] < 0)]]

special_buffered_select=which(res_te$padj < 0.05 & abs(res_te$log2FoldChange) >= log2(fold_change_threshold) &
res_ribo$padj >= 0.05 & 
res_rna$padj < 0.05 & abs(res_rna$log2FoldChange) >= log2(fold_change_threshold ))

buffered = c(rownames(res_te)[special_buffered_select],buffered)

# join the 3 results objects
classified_deseq <- bind_cols(gene_id=rownames(res_rna), rna_fc =res_rna$log2FoldChange, rna_padj=res_rna$padj,
rpf_fc = res_ribo$log2FoldChange, rpf_padj=res_ribo$padj,te_fc=res_te$log2FoldChange, te_padj=res_te$padj)

classified_events <- classified_deseq %>%
mutate(class=case_when(
 gene_id %in% transcriptional  ~ "forwarded",
 gene_id %in% translation_exclusive  ~ "translation exclusive",
 gene_id %in% intensified  ~ "intensified",
 gene_id %in% buffered  ~ "buffered",
 TRUE ~ "undetermined"
))

print(table(classified_events$class))

# maoa is causing a problem    
moao_annoation <- classified_events %>% 
filter(gene_id=="Maoa") %>%
mutate(symbol = "Maoa",
GeneID = 100760815,
name = "amine oxidase [flavin-containing] A")


pcg_classified <- classified_events %>%
  filter(!str_detect(gene_id, "NR|XR"))%>% 
  filter(class != "undetermined") %>%
  mutate(symbol = gsub("_.*", "", gene_id)) %>%
  left_join(CriGri_PICRH_1_0_annotation, by = "symbol", relationship = "many-to-many") %>%
  filter(`# feature` == "mRNA") %>%
  mutate(name = gsub(",.*", "", name)) %>%
  distinct(gene_id, .keep_all = T) %>%
  dplyr::rename(class=class.y) %>%
  dplyr::select("symbol", "name", "GeneID", c(colnames(classified_events)))

# determine genes successfully assigned to gene symbols
pcg_without_loc_ids <- pcg_classified %>%
  filter(!str_detect(symbol, "LOC"))

# determine genes with no symbol annotation
pcg_with_loc_ids <- pcg_classified %>%
  filter(str_detect(symbol, "LOC"))

# match the long gene name against mouse using fuzzy match
loc_id_gene_symbol <- mouse_feature_table %>%
  mutate(ncbi_symbol = symbol) %>%
  fuzzy_inner_join(pcg_with_loc_ids, by = "name", match_fun = str_detect) %>%
  distinct(symbol.y, .keep_all = T)

pcg_with_loc_ids <- pcg_with_loc_ids %>%
mutate("symbol.y" = symbol) %>%
left_join(loc_id_gene_symbol, by = "symbol.y") %>%
mutate(ncbi_symbol = coalesce(ncbi_symbol, symbol.y)) %>%
mutate(
  "gene_id" = gene_id.x,
  "rna_fc" = rna_fc.x,
  "rna_padj" = rna_padj.x,
  "rpf_fc" = rpf_fc.x, 
  "rpf_padj" = rpf_padj.x,
  "te_fc" = te_fc.x,
  "te_padj" = te_padj.x,
  "symbol" = ncbi_symbol
) %>%
dplyr::select(gene_id, symbol, name, GeneID, rna_fc, rna_padj, rpf_fc, rpf_padj, te_fc, te_padj, class)

# annotate new ORFs
nc_sig_res <- classified_events %>%
filter(class != "undetermined") %>%
filter(str_detect(gene_id, "NR|XR")) %>%
mutate(symbol = "New", name = "New",
        GeneID = 0) %>%
dplyr::select(gene_id, symbol, name, GeneID, rna_fc, rna_padj, rpf_fc, rpf_padj, te_fc, te_padj, class )

classified_event_annotation <-
bind_rows(moao_annoation, pcg_without_loc_ids, pcg_with_loc_ids, nc_sig_res) %>%
dplyr::select(c("gene_id", "GeneID", "symbol", "name")
) 


saveRDS(classified_events, file = paste(results_dir, "/classified_events.rds", sep = ""))
saveRDS(classified_event_annotation, file = paste(results_dir, "/classified_event_annotation.rds", sep = ""))