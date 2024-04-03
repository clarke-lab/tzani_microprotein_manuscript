run_deseq <-
  function(count_data,
           sample_table,
           fold_change_thresh,
           pval_thresh,
           base_mean_thresh,
           average_counts,
           type,
           CriGri_PICRH_1_0_annotation,
           mouse_feature_table) {
    
    
    if (type == "translation") {
      rna_counts <- count_data %>%
        dplyr::select(contains("rna"))
      ribo_counts <- count_data %>%
        dplyr::select(contains("ribo"))
      
      detected.rna <-
        rownames(rna_counts)[rowMeans(rna_counts) >= average_counts]
      
      detected.ribo <-
        rownames(rna_counts)[rowMeans(ribo_counts) >= average_counts]
      detected.both <- intersect(detected.rna, detected.ribo)
      
      print(length("detected.both"))
      
      counts <- count_data[detected.both, ]
      
    } else {
      counts <- count_data %>%
        dplyr::select(contains(type))
      
      detected.genes <-
        rownames(counts)[rowMeans(counts) >= average_counts]
      
      counts <- counts[detected.genes,]
    }
    
    if (type == "translation") {
      dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = sample_table,
        design = ~ condition + assay + condition:assay
      )
      
      dds$condition <- relevel(dds$condition, ref = "nts")
      dds$assay <- relevel(dds$assay, ref = "rnaseq")
    } else {
      design <- sample_table %>%
        filter(assay == type)
      
      dds <-
        DESeqDataSetFromMatrix(
          countData = counts,
          colData = design,
          design = ~ condition
        )
      
      dds$condition <- relevel(dds$condition, ref = "nts")
    }
    
    print(head(counts))
    print(design)
    
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    dds <- DESeq(dds)
    
    # track plotting scaling
    # bam_coverage_scaling_rnaseq <- 1 / sizeFactors(dds_rnaseq)
    # write_lines(bam_coverage_scaling_rnaseq, paste(results_dir, "rnaseq_size_factors.txt", sep = ""))
    
    if (type == "translation") {
      res = results(dds, name = "conditionts.assayriboseq")
    } else {
      res <- results(dds)
    }
    
    sig_res <- res %>%
      as_tibble(rownames = "geneid") %>%
      dplyr::filter(
        abs(log2FoldChange) >= log2(fold_change_thresh) &
          padj < pval_thresh & baseMean >= base_mean_thresh
      )
    
    # annotate protein coding genes
    pcg_sig_res <- sig_res %>%
      filter(!str_detect(geneid, "NR|XR")) %>%
      mutate(symbol = gsub("_.*", "", geneid)) %>%
      left_join(CriGri_PICRH_1_0_annotation, by = "symbol") %>%
      filter(`# feature` == "mRNA") %>%
      mutate(name = gsub(",.*", "", name)) %>%
      distinct(geneid, .keep_all = T) %>%
      dplyr::select("symbol", "name", "GeneID", c(colnames(sig_res)))
    
    # gene names causes a problem with the regex in fuzzy match
    # pcg_sig_re <- pcg_sig_res %>%
    #   mutate(symbol = replace(symbol, name == "amine oxidase [flavin-containing] A", "Maoa"))
    #
    
    pcg_without_loc_ids <- pcg_sig_res %>%
      filter(!str_detect(symbol, "LOC"))
    
    pcg_with_loc_ids <- pcg_sig_res %>%
      filter(str_detect(symbol, "LOC"))
    
    # match the names against mouse
    loc_id_gene_symbol <- mouse_feature_table %>%
      mutate(ncbi_symbol = symbol) %>%
      fuzzy_inner_join(pcg_with_loc_ids, by = "name", match_fun = str_detect) %>%
      distinct(symbol.y, .keep_all = T)
    
    pcg_with_loc_ids <- pcg_with_loc_ids %>%
      mutate("symbol.y" = symbol) %>%
      left_join(loc_id_gene_symbol, by = "symbol.y") %>%
      mutate(ncbi_symbol = coalesce(ncbi_symbol, symbol.y)) %>%
      mutate(
        "geneid" = geneid.x,
        "baseMean" = baseMean.x,
        "log2FoldChange" = log2FoldChange.x,
        "lfcSE" = lfcSE.x,
        "stat" = stat.x,
        "pvalue" = pvalue.x,
        "padj" = padj.x,
        "symbol" = ncbi_symbol
      ) %>%
      dplyr::select(symbol, name, GeneID, c(colnames(sig_res)))
    
    nc_sig_res <- sig_res %>%
      filter(str_detect(geneid, "NR|XR")) %>%
      mutate(symbol = "New",
             name = "New",
             GeneID = 0) %>%
      dplyr::select("symbol",
                    "name",
                    "GeneID",
                    c(colnames(sig_res)))
    
    sig_res <-
      bind_rows(pcg_without_loc_ids, pcg_with_loc_ids, nc_sig_res) %>%
      arrange(-log2FoldChange) %>%
      dplyr::select(
        c(
          "geneid",
          "GeneID",
          "symbol",
          "name",
          "baseMean",
          "log2FoldChange",
          "lfcSE",
          "stat",
          "pvalue",
          "padj"
        )
      ) %>%
      dplyr::rename(Plastid_ID = "geneid")
    
    return(list(
      deseq_object = dds,
      deseq_res = res,
      significant = sig_res
    ))
    
  }