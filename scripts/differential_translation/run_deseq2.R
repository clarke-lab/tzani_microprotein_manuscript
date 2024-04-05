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

############################################################
# 2. Assess deltaTE classes
############################################################

te_res <- differential_translation$deseq_res
rna_res <- differential_expression$deseq_res[rownames(te_res),]
rpf_res <- differential_rpf$deseq_res[rownames(te_res),]

forwarded = rownames(te_res)[which(te_res$padj >= 0.05 & rpf_res$padj < 0.05 & rna_res$padj < 0.05)]  # transcription only
exclusive = rownames(te_res)[which(te_res$padj < 0.05 & rpf_res$padj < 0.05 & rna_res$padj >= 0.05)]  # translation only 
both = rownames(te_res)[which(te_res$padj < 0.05 & rpf_res$padj < 0.05 & rna_res$padj < 0.05)] # change in both transcription and translation

intensified = rownames(te_res[both[which(te_res[both,2]*rna_res[both,2] > 0)],]) # transcription only

buffered = rownames(te_res[both[which(te_res[both,2]*rna_res[both,2] < 0)],])
buffered = c(rownames(te_res)[which(te_res$padj < 0.05 & rpf_res$padj >= 0.05 & rna_res$padj < 0.05)], buffered)

rna_res <- rna_res[rownames(rpf_res),]
joined <- bind_cols(rpf_fc = rpf_res$log2FoldChange, rna_fc =rna_res$log2FoldChange)
joined <- joined %>% dplyr::mutate(gene_id=rownames(rpf_res))

joined <- joined %>%
mutate(class=case_when(
 gene_id %in% forwarded  ~ "forwarded",
 gene_id %in% exclusive  ~ "exclusive",
 gene_id %in% intensified  ~ "intensified",
 gene_id %in% buffered  ~ "buffered",
 TRUE ~ "undetermined"
))

joined <- joined %>%
summarise(pmax(rna_fc, rpf_fc))

colorblind_palette <- c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#999999")
max_fc <- max(c(abs(joined$rna_fc), abs(joined$rpf_fc)), na.rm = T)
max_fc =6
joined %>%
arrange(desc(class)) %>%
ggplot(aes(x=rna_fc, y=rpf_fc, color=class)) +
geom_abline(intercept = 0, slope = 1, color = "light grey") +  # Add diagonal line through origin
geom_vline(xintercept=0, color = "light grey") +
geom_hline(yintercept=0, color = "light grey") +
geom_point(size=0.8, alpha=1) +
scale_color_manual(values=colorblind_palette) +
geom_text_repel(data = subset(joined, str_detect(gene_id, "XR") & class != "undetermined" & class != "forwarded"), aes(label = gene_id), box.padding = 0.5) +
labs(x = expression("Log"[2] ~ "mRNA FC (TS/NTS)"), 
y = expression("Log"[2] ~ "RPF FC (TS/NTS)"), color = "", alpha = "") +
xlim(-max_fc, max_fc) +
ylim(-max_fc, max_fc) +
theme_bw() +
guides(color = guide_legend(override.aes = list(size = 3)))




b <- joined %>%
filter(class != "undetermined") %>%
filter(str_detect(gene_id, "XR")

ggscatterhist(
  joined, x = "rna_fc", y = "rpf_fc",
  color = "class", size = 0.1, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07", "green", "blue"),
  margin.params = list(fill = "class", color = "black", size = 0.2)
  )

# save(differential_expression,differential_rpf,differential_translation, file = paste(results_dir, "/deseq_output.RData", sep = ""))