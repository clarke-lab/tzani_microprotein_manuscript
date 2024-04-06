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
scale_factors_rnaseq <- data.frame(colnames(rna_deseq$deseq_object),sizeFactors(rna_deseq$deseq_object)^-1)
write_tsv(scale_factors_rnaseq, file = paste(results_dir, "rna_scale_factors.txt", sep = "/"), col_names =F)

scale_factors_riboseq <- data.frame(colnames(rpf_deseq$deseq_object),sizeFactors(rpf_deseq$deseq_object)^-1)
write_tsv(scale_factors_riboseq, file = paste(results_dir, "rpf_scale_factors.txt", sep = "/"), col_names =F)

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

special_buffered_select=which(res_te$padj < 0.05 & abs(res_te$log2FoldChange) >= log2(fold_change_threshold ) &
res_ribo$padj >= 0.05 & 
res_rna$padj < 0.05 & abs(res_rna$log2FoldChange) >= log2(fold_change_threshold ))

buffered = c(rownames(res_te)[special_buffered_select],buffered)

# join the 3 results objects
classified_deseq <- bind_cols(gene_id=rownames(res_rna), rna_fc =res_rna$log2FoldChange, rna_padj=res_rna$padj,
rpf_fc = res_ribo$log2FoldChange, rpf_padj=res_ribo$padj,te_fc=res_te$log2FoldChange, te_pdj=res_te$padj)

classified_events <- classified_deseq %>%
mutate(class=case_when(
 gene_id %in% transcriptional  ~ "forwarded",
 gene_id %in% translation_exclusive  ~ "translation exclusive",
 gene_id %in% intensified  ~ "intensified",
 gene_id %in% buffered  ~ "buffered",
 TRUE ~ "undetermined"
))

saveRDS(classified_events, file = paste(results_dir, "/classified_events.rds", sep = ""))



# colorblind_palette <- c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#999999")
# max_fc =max(abs(c(joined$rna_fc, joined$rpf_fc)))

# joined %>%
# arrange(desc(class)) %>%
# ggplot(aes(x=rna_fc, y=rpf_fc, color=class)) +
# geom_abline(intercept = 0, slope = 1, color = "grey") +  # Add diagonal line through origin
# geom_vline(xintercept=c(-log2(1.5),0, log2(1.5)), color = "grey") +
# geom_hline(yintercept=c(-log2(1.5),0, log2(1.5)), color = "grey") +
# geom_point(size=0.6, alpha=0.8) +
# scale_color_manual(values=colorblind_palette) +
# geom_text_repel(data = subset(joined, str_detect(gene_id, "XR") & class != "undetermined" & class != "forwarded"), aes(label = gene_id), box.padding = 0.5) +
# labs(x = expression("Log"[2] ~ "mRNA FC (TS/NTS)"), 
# y = expression("Log"[2] ~ "RPF FC (TS/NTS)"), color = "", alpha = "") +
# xlim(-max_fc, max_fc) +
# ylim(-max_fc, max_fc) +
# theme_bw() +
# guides(color = guide_legend(override.aes = list(size = 3)))

# b <- joined %>%
# filter(class != "undetermined") %>%
# filter(str_detect(gene_id, "XR|NR"))

# ggscatterhist(
#   joined, x = "rna_fc", y = "rpf_fc",
#   color = "class", size = 0.1, alpha = 0.6,
#   palette = c("#00AFBB", "#E7B800", "#FC4E07", "green", "blue"),
#   margin.params = list(fill = "class", color = "black", size = 0.2)
#   )

# # save(differential_expression,rpf_deseq,te_deseq, file = paste(results_dir, "/deseq_output.RData", sep = ""))


# ggplot(data.frame(res_rna), aes(x = log2FoldChange, y = -log10(pvalue))) +
#   geom_point(size=0.1) +
#   geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
#   scale_color_manual(values = c("black" = "black", "red" = "red")) +
#   labs(x = "log2(Fold Change)", y = "-log10(p-value)", title = "Volcano Plot") +
#   theme_minimal()