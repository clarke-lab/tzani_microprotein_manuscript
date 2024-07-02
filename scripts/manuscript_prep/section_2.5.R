## ----setup---------------------------------------------------------------------------------------------------------------------
package_list <- c(
  "tidyverse", "DESeq2", "writexl", "viridis", "ggpp", "scales",
  "heatmaply", "ggvenn", "WebGestaltR", "ggrepel", "cowplot",
  "patchwork", "ggpubr", "GenomicFeatures", "wiggleplotr"
)

suppressMessages(lapply(package_list, require, character.only = TRUE))


## ------------------------------------------------------------------------------------------------------------------------------
source(file = "scripts/manuscript_prep/utility_functions.R")


## ----make_results_dir----------------------------------------------------------------------------------------------------------
results_dir <- "manuscript/section_2.5/"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


## ----load_data-----------------------------------------------------------------------------------------------------------------
orftable <- readRDS("orf_identification/orf_filtered/final_orfs.rds")

rna_analysis <- readRDS(file = "differential_translation/deseq_output/rna_deseq.rds")

rpf_analysis <- readRDS(file = "differential_translation/deseq_output/rpf_deseq.rds")

te_analysis <- readRDS(file = "differential_translation/deseq_output/te_deseq.rds")

classified_events <- readRDS(file = "differential_translation/deseq_output/classified_events.rds")

classified_event_annotation <- readRDS(file = "differential_translation/deseq_output/classified_event_annotation.rds")


## ----supp_figure_10------------------------------------------------------------------------------------------------------------
vsd <- vst(te_analysis$deseq_object, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

ngs_dendrogram <- heatmaply::ggheatmap(sampleDistMatrix)

ggsave(ngs_dendrogram, filename = paste(results_dir, "Supplementary Figure 10.png", sep = ""), 
       width = 8, height = 8, device = "png", dpi = 2000)


## ----new_orf_summary-----------------------------------------------------------------------------------------------------------
new_orf_classification <- orftable %>%
  filter(`ORF type` == "New") %>%
  mutate(rna_type = case_when(
    str_detect(`ORF-RATER name`, "XR_|NR_") ~ "ncRNA",
    !str_detect(`ORF-RATER name`, "XR_|NR_") ~ "mRNA"
  )) %>%
  mutate(orf_length = case_when(
    `Length (AAs)` >= 5 & `Length (AAs)` & `Length (AAs)` < 100 &
      rna_type == "ncRNA" ~ "lncRNA ORF < 100aa",
  rna_type == "mRNA" ~ "mRNA",
    `Length (AAs)` < 5 | `Length (AAs)` &
      `Length (AAs)` >= 100 &
      rna_type == "ncRNA" ~ "lncRNA ORF >= 100aa"
  )) %>%
  mutate(
    type = "New ORFs",
    orf_length = factor(orf_length, levels = c("mRNA", "lncRNA ORF >= 100aa", "lncRNA ORF < 100aa"))
  )

paste0("_________________________________________________________________________________ ")
paste0("                         New ORF length/transcript type                           ")
paste0("__________________________________________________________________________________")
table(new_orf_classification$rna_type, new_orf_classification$orf_length)


## ----figa5a--------------------------------------------------------------------------------------------------------------------
new_orf_classification %>%
  ggplot(aes(x = type, fill = orf_length)) +
  geom_bar(aes(y = after_stat(count) / sum(..count..)), position = "stack", width = 1) +
  theme_minimal() +
  scale_fill_manual(values = c("#999999", "#CC79A7", "#35968CFF")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(fill = "", x = "", y = "Proportion of ORFs")

ggsave(filename = paste(results_dir, "Figure 5a.png", sep = ""), 
       dpi = 2000, device = "png", width = 3, height = 3)


## ----fig5b---------------------------------------------------------------------------------------------------------------------
lncRNA_sorfs <- orftable %>%
  filter(`ORF type` == "New" &
    str_detect(`ORF-RATER name`, "XR|NR") &
    `Length (AAs)` >= 5 & `Length (AAs)` < 100) %>%
  group_by(`Transcript ID`) %>%
  mutate(num_ORFs = n()) %>%
  ungroup()

length(lncRNA_sorfs$`Transcript ID`)


lncRNA_sorfs %>%
  ggplot(aes(x = sort(`Length (AAs)`))) +
  geom_bar(size=0.1, fill = "#35968CFF", color = "black") +
  labs(y = "# lncRNA sORFs", x = "ORF Length (aa)") +
  geom_text_npc(aes(npcx = 0.9, npcy = 0.9, label = paste0("Mean length:\n", round(mean(`Length (AAs)`)), " aa")), size = 3) +
  theme_bw()

ggsave(filename = paste(results_dir, "Figure 5b.png", sep = ""), 
       dpi = 2000, device = "png", width = 3, height = 3)


## ----fig5c---------------------------------------------------------------------------------------------------------------------
new_ORFs_per_transcript_plot <- lncRNA_sorfs %>%
  dplyr::select(`Transcript ID`, num_ORFs) %>%
  ggplot(aes(x = num_ORFs)) +
  geom_bar(stat = "count", fill = "#35968CFF", color = "black", size=0.1) +
  scale_x_discrete(limits = as.character(c(1:max(lncRNA_sorfs$num_ORFs)))) +
  theme_bw() +
  labs(x = "ORFs per transcript", y = "Number of transcripts") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.2, color = "black", size = 3)
new_ORFs_per_transcript_plot

ggsave(plot = new_ORFs_per_transcript_plot, filename = paste(results_dir, "Figure 5c.png", sep = ""), 
       dpi = 2000, device = "png", width = 3, height = 3)


## ------------------------------------------------------------------------------------------------------------------------------
paste0("_________________________________________________________________________________ ")
paste0("                         TE total features                                         ")
paste0("__________________________________________________________________________________")
dim(te_analysis$deseq_res)

paste0("_________________________________________________________________________________ ")
paste0("                         differential TE                                         ")
paste0("__________________________________________________________________________________")
dim(te_analysis$significant)


## ------------------------------------------------------------------------------------------------------------------------------
paste0("_________________________________________________________________________________ ")
paste0("                         classes identified by deltaTE                            ")
paste0("__________________________________________________________________________________")

dim(classified_events %>% filter(class != "undetermined") )[1]

classified_events %>%
  dplyr::count(class) %>%
  filter(class != "undetermined") %>%
  mutate(proportion = n / sum(n))



## ----fig5d---------------------------------------------------------------------------------------------------------------------
translation_regulation_palette <- c("#377eb8",  "#4daf4a","#e41a1c", "#984ea3", "#999999")

max_fc =max(abs(c(classified_events$rna_fc, classified_events$rpf_fc)))

classified_events %>%
   arrange(desc(class)) %>%
   ggplot(aes(x=rna_fc, y=rpf_fc, color=class)) +
   geom_abline(intercept = 0, slope = 1, color = "grey") +  # Add diagonal line through origin
   geom_vline(xintercept=0, color = "light grey", linetype="dashed", size=0.2) +
   geom_hline(yintercept=0, color = "light grey", linetype="dashed", size=0.2) +
   geom_point(size=0.1, alpha=0.5) +
   scale_color_manual(values= translation_regulation_palette, labels=label_wrap(10)) +
   labs(x = expression("Log"[2] ~ "mRNA FC (TS/NTS)"), 
        y = expression("Log"[2] ~ "RPF FC (TS/NTS)"), color = "", alpha = "") +
   theme_bw() +
   scale_x_continuous(expand = c(0, 0)) +
   scale_y_continuous(expand = c(0, 0)) +
   xlim(-7.6, 7.6) +
   ylim(-7.6,7.6) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   guides(color = guide_legend(override.aes = list(size = 2, alpha=c(1,1,1,1,1))))
 
ggsave(filename = paste0(results_dir, "Figure 5d.png"), 
        dpi = 2000, device = "png", width = 5, height = 3.5)


## ----de_summary----------------------------------------------------------------------------------------------------------------
paste0("_________________________________________________________________________________ ")
paste0("                         DEG with |FC| > 1.5 fold.                                ")
paste0("__________________________________________________________________________________")
DEG <- classified_events %>%
  mutate(`ORF type` = case_when(
    str_detect(gene_id, "XR_|NR_") ~ "sORF",
    TRUE ~ "canonical"
  )) %>%
filter( class %in% c("forwarded")) %>%
filter(abs(rna_fc) > log2(1.5) & rna_padj < 0.05 &
abs(rpf_fc) > log2(1.5) & rpf_padj < 0.05) %>%
left_join(classified_event_annotation, by="gene_id") %>%
filter(!is.na(symbol)) %>%
dplyr::select(gene_id,GeneID, symbol,name,class,everything())
print(paste0(dim(DEG)[1], " DEGs are significant"))

paste0("_________________________________________________________________________________ ")
paste0("                         DTEG with |TE| > 1.5 fold.                                ")
paste0("__________________________________________________________________________________")
DTEG <- classified_events %>%
  mutate(`ORF type` = case_when(
    str_detect(gene_id, "XR_|NR_") ~ "sORF",
    TRUE ~ "canonical"
  )) %>%
filter( class %in% c("buffered", "intensified", "translation exclusive")) %>%
filter(abs(te_fc) > log2(1.5) & te_padj < 0.05) %>%
left_join(classified_event_annotation, by="gene_id") %>%
filter(!is.na(symbol)) %>%
dplyr::select(gene_id,GeneID, symbol,name,class,everything())
print(paste0(dim(DTEG)[1], " DEGs are significant"))

merged_significant <- bind_rows(DEG, DTEG)


## ----supp_data_6---------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 6.xlsx", sep = "")
suppressMessages(if (suppressMessages(file.exists(fn))) {
file.remove(fn)})

differential_trannslation_out_list <- merged_significant %>%
  dplyr::select(gene_id, GeneID, symbol, name, `ORF type`, everything()) %>%
  dplyr::rename("Gene ID"=gene_id, "NCBI Gene ID" = GeneID, "Gene symbol" = symbol, "Gene Name" = name) %>%
  arrange(desc(`ORF type`)) %>%
group_split(class)

write_xlsx(list(
  a = differential_trannslation_out_list[[2]], # forwarded
  b = differential_trannslation_out_list[[1]], # buffered
  c = differential_trannslation_out_list[[3]], # intensified
  d = differential_trannslation_out_list[[4]]  # translation exclusive
), path = fn, format_headers = TRUE)


## ----go_analysis---------------------------------------------------------------------------------------------------------------
enrichdir <- paste0(results_dir, "enrichment_analysis/")
suppressMessages(if (file.exists(enrichdir)) {
  unlink(enrichdir)
})
if (!dir.exists(enrichdir)) {
  dir.create(enrichdir)
}

go_gene_symbols <- merged_significant$gene_id

write(go_gene_symbols, file = paste0(enrichdir, "gene_symbols.txt"))

enrich_result <- WebGestaltRBatch(
  enrichMethod = "ORA",
  organism = "mmusculus",
  hostName = "https://www.webgestalt.org",
  enrichDatabase = c("geneontology_Biological_Process_noRedundant"),
  enrichDatabaseType = "genesymbol",
  interestGeneFolder = enrichdir,
  interestGeneType = "genesymbol",
  referenceSet = "genome",
  minNum = 10,
  maxNum = 500,
  sigMethod = "fdr",
  fdrMethod = "BH",
  fdrThr = 0.05,
  topThr = 10,
  reportNum = 20,
  perNum = 1000,
  projectName = "diff_trans",
  isOutput = TRUE,
  outputDirectory = enrichdir,
  dagColor = "continuous",
  setCoverNum = 10,
  networkConstructionMethod = NULL,
  neighborNum = 10,
  highlightType = "Seeds",
  highlightSeedNum = 10,
  nThreads = 32
)

te_gene_enrichment <- enrich_result[[1]]$enrichResult[, c(1, 2, 4, 5, 7, 8, 9, 11)] %>%
  filter(FDR < 0.05)


## ------------------------------------------------------------------------------------------------------------------------------
class_proportions <- data.frame()
for (i in 1:dim(te_gene_enrichment)[1]) {
  gene_vector <- unlist(strsplit(te_gene_enrichment$userId[i], ";"))
  
 current_data <- merged_significant %>%
    filter(gene_id %in% gene_vector) %>%
    group_by(class) %>%
    summarise(count = n()) %>%
   mutate(proportion = count / sum(count)) %>%
   mutate(process=te_gene_enrichment$description[i]) %>%
   mutate(FDR=te_gene_enrichment$FDR[i],
          geneSet=te_gene_enrichment$geneSet[i]) %>% 
    dplyr::select(process,class, count, proportion, FDR, geneSet)
 
 current_data <- current_data  %>%
   mutate(forwarded_prop = current_data[current_data$class == "forwarded",]$proportion)
 
 class_proportions <- bind_rows(class_proportions, current_data)
}


top_processes <- class_proportions %>% 
  arrange(forwarded_prop) %>%
  distinct(process,.keep_all = T)
paste0("_____________________________________________________________________________________ ")
paste0("GO biological processes with > 25%  of genes with evidence of translational regulation")
paste0("_____________________________________________________________________________________ ")
length(top_processes$process[top_processes$forwarded_prop < 0.75])


## ----fig5e---------------------------------------------------------------------------------------------------------------------
translation_regulation_palette <- c("#4daf4a",  "#e41a1c","#377eb8", "#984ea3")

class_proportions %>%
  arrange(forwarded_prop) %>%
  filter(process %in% top_processes$process[1:10]) %>% 
  ggplot(aes(x=reorder(process, -forwarded_prop), y=proportion, fill=reorder(class, -proportion))) +
  geom_col(alpha=0.95) +
  theme_minimal_vgrid() +
  coord_flip() +
  labs(x="", y = "Proportion of genes in GO category", fill="") +
  scale_fill_manual(values=translation_regulation_palette, labels=label_wrap(10)) + 
  theme(legend.position = "right", axis.text.y = element_text(size = 14)) +
  scale_x_discrete(labels = label_wrap(30))  

ggsave(filename = paste0(results_dir, "Figure 5e.png"), 
        dpi = 2000, device = "png", width = 9, height = 6)


## ----supp_data_7---------------------------------------------------------------------------------------------------------------
fn <- paste(results_dir, "Supplementary Data 7.xlsx", sep = "")

suppressMessages(if (suppressMessages(file.exists(fn))) {
file.remove(fn)})

translation_portions <- class_proportions %>%
  mutate(`% translationally regulated` = round((1-forwarded_prop)*100,2)) %>%
  distinct(process, .keep_all = T) %>%
  dplyr::select(geneSet, `% translationally regulated`)

te_gene_enrichment <- te_gene_enrichment %>%
left_join(translation_portions, by="geneSet") %>%
dplyr::select(geneSet, description, size, overlap, `% translationally regulated`, everything()) %>%
rename(GO_ID = geneSet, `Biological Procress`=description, `Genes in category` = size,`DTG or DTEG` = overlap) %>%
  arrange(-`% translationally regulated`)

write_xlsx(list(
a = te_gene_enrichment
), path = fn, format_headers = TRUE)



## ----fig5f---------------------------------------------------------------------------------------------------------------------
# for DTG select only those with a RNA and RPF FC > 1.5 
rpf_analysis_classified <- as_tibble(rpf_analysis$deseq_res, rownames="gene_id") %>% 
  filter(gene_id %in% classified_events$gene_id) %>%
  mutate(class = case_when(
    (gene_id %in% DEG$gene_id) & !str_detect(gene_id, "XR_|NR_") ~ "canonical ORF",
    gene_id %in% DEG$gene_id & str_detect(gene_id, "XR_|NR_") ~ "sORF",
    TRUE ~ "Not significant"
  )) %>%
   arrange(class) %>%
  mutate(point_size = case_when(
    class == "canonical ORF" ~ 0.2,
    class == "sORF" ~ 1,
    TRUE ~ 0.1
  )) %>%
  mutate(point_alpha = case_when(
    class == "canonical ORF" ~ 0.8,
    class == "sORF" ~ 0.8,
    TRUE ~ 0.1
  ))

paste0("_________________________________________________________________________________ ")
paste0("                        Significant with sORFs                                           ")
paste0("__________________________________________________________________________________")
table(rpf_analysis_classified$class)

rpf_analysis_classified %>%
ggplot(aes(x=log2FoldChange, y=-log10(pvalue), color=class)) + 
    geom_point(size=rpf_analysis_classified$point_size, 
               alpha=rpf_analysis_classified$point_alpha) +
  geom_vline(xintercept=c(-log2(1.5),0, log2(1.5)), 
             color = "light grey", 
             linetype="dashed", 
             size=0.2) +
  geom_hline(yintercept=-log10(0.05), 
             color = "light grey", 
             linetype="dashed", 
             size=0.2) +
  scale_color_manual(values = c("#CC6677","light grey",  "#035968"), labels=label_wrap(10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  labs(x = "Log"[2] ~ "RPF FC (TS/NTS)", color = "")

ggsave(filename = paste0(results_dir, "Figure 5f.png"), 
        dpi = 2000, device = "png", width = 3.1, height = 3.5)


## ----supp_figure_11------------------------------------------------------------------------------------------------------------
forwarded_events <- classified_events %>%
  filter(str_detect(gene_id,"XR_|NR_")) %>%
  filter(class == "forwarded") %>%
  filter(abs(rna_fc) > log2(1.5) & rna_padj < 0.05 & abs(rpf_fc) > log2(1.5) & rpf_padj < 0.05) %>%
  dplyr::select(gene_id) 

goi <- forwarded_events$gene_id

tcounts <- t(log2((counts(te_analysis$deseq_object[goi, ], normalized = TRUE, replaced = FALSE) + .5))) %>%
  merge(colData(te_analysis$deseq_object), ., by = "row.names")

tcounts <- tcounts %>% 
  pivot_longer(cols = -c(condition, assay, Row.names, sizeFactor)) %>%
  mutate(deltaTE_class = "forwarded") %>%
  mutate(label_facet=paste0(name,"\n",deltaTE_class)) %>%
  mutate(label_facet=gsub(",",",\n", label_facet)) %>%
  mutate(x_value= paste0(condition,"\n",assay))

 tcounts %>% 
  ggplot(aes(x = x_value, y = value, fill = condition)) +
  geom_boxplot() +
  geom_jitter(width=0.05) +
  facet_wrap(~label_facet, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(x = "", y = expression("Log"[2] ~ "DESeq2 normalised counts"), fill = "") +
  scale_fill_manual(values = c("#D55E00", "#56B4E9")) +
  theme(strip.text = element_text(face = "italic", size = 8), 
        strip.background = element_blank(),
        axis.title = element_text(size = 9), 
        legend.position = "bottom", 
        axis.title.y = element_text(size = 12)
        ) + 
   scale_x_discrete(limits=c("nts\nrnaseq", "ts\nrnaseq", "nts\nriboseq", "ts\nriboseq"),
                    labels=c("NTS\nRNA-seq","TS\nRNA-seq", "NTS\nRibo-seq", "TS\nRibo-seq" ))

ggsave(filename = paste0(results_dir, "Supplementary Figure 11.png"), 
        dpi = 1200, device = "png", width = 8.5, height = 12)


## ------------------------------------------------------------------------------------------------------------------------------
translation_regulation_events <- classified_events %>%
  filter(str_detect(gene_id,"XR_|NR_")) %>%
  filter(!class == "forwarded" & !class == "undetermined") %>%
  dplyr::select(gene_id)

goi <- translation_regulation_events$gene_id

tcounts <- t(log2((counts(te_analysis$deseq_object[goi, ], normalized = TRUE, replaced = FALSE) + .5))) %>%
  merge(colData(te_analysis$deseq_object), ., by = "row.names")

sorf_te_data <- bind_rows(
  tcounts[tcounts$condition == "nts" & tcounts$assay == "riboseq", goi] /
    tcounts[tcounts$condition == "nts" & tcounts$assay == "rnaseq", goi],
  tcounts[tcounts$condition == "ts" & tcounts$assay == "riboseq", goi] /
    tcounts[tcounts$condition == "ts" & tcounts$assay == "rnaseq", goi]
) %>%
  mutate(condition = toupper(tcounts$condition[1:8])) %>%
  pivot_longer(cols = -condition) %>%
  mutate(deltaTE_class = case_when(
    name=="XR_003481490.2_407599334_59aa" ~ "buffered",
    name=="XR_003483976.2_53501071_63aa" ~ "buffered",
    name=="XR_003484005.2_59719147_46aa" ~ "intensified",
    name=="XR_004771714.1_255645447_91aa" ~ "translation exclusive",
    name=="XR_004770567.1_82548222_17aa" ~ "translation exclusive",
  )) %>%
  mutate(label_facet=paste0(name,"\n",deltaTE_class))
  
  
 sorf_te_data$name <- factor(sorf_te_data$name, levels = c('XR_003481490.2_407599334_59aa', 'XR_003483976.2_53501071_63aa', 
                                                           'XR_003484005.2_59719147_46aa',  'XR_004771714.1_255645447_91aa',
                                                           'XR_004770567.1_82548222_17aa'), 
                   labels = unique(sorf_te_data$label_facet))
 
sorf_te_data %>% 
ggplot(aes(x = condition, y = value, fill = condition)) +
geom_boxplot() +
geom_point(size = 0.5) +
facet_wrap(~name, scales = "free_y", ncol = 3) +
theme_bw() +
labs(x = "", y = expression("Log"[2] ~ "TE (TS/NTS)"), fill = "") +
scale_fill_manual(values = c("#D55E00", "#56B4E9")) +
theme(strip.text = element_text(face = "italic", size = 8), 
      strip.background = element_blank(),
      axis.title = element_text(size = 6), 
      legend.position = "bottom", 
      axis.title.y = element_text(size = 12),
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.ticks.x = element_blank())

ggsave(filename = paste(results_dir, "Figure 5g.png", sep = ""), 
       width = 7, height = 5, device = "png", dpi = 2000)


## ------------------------------------------------------------------------------------------------------------------------------
print("section 2.5 results complete")

