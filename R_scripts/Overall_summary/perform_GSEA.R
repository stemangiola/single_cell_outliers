# Execute `GSEA` using `MSigDB`, `clusterProfiler` and `enrichplot` 
# import COVID_19 dataset all_de.rds
all_de <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/all_de.rds")

de_all_gene_rank = all_de %>%
    # symbol_to_entrez(.transcript = transcript,
    #                  .sample = sample) %>% 
    mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys = transcript,
                                          keytype = "SYMBOL",
                                          column = "ENTREZID",
                                          multiVals = "first"
    )) %>% filter(method == "edgeR_quasi_likelihood") %>% 
    filter(edgerQLT_PValue   %>% is.na %>% `!`) %>%
    tidybulk::test_gene_rank(
        .entrez = entrez,
        .arrange_desc = edgerQLT_logFC,
        # .arrange_desc = logFC,
        .sample = sample,
        species="Homo sapiens",
        gene_sets = c("H", "C2", "C5")
    )


# Examine significantly enriched gene sets

de_all_gene_rank %>%
    filter(gs_cat  == "C2" ) %>%
    dplyr::select(-fit) %>%
    unnest(test) %>% 
    filter(p.adjust < 0.05)


# Visualise enrichment


de_all_gene_rank_plots = 
    de_all_gene_rank  %>% 
    unnest(test) %>% 
    
    # Select top 10
    slice(1:10) %>% 
    mutate(plot = pmap(
        list(fit, ID, idx_for_plotting, p.adjust), 
        ~ enrichplot::gseaplot2(
            ..1, 
            geneSetID = ..3, 
            title = sprintf("%s \nadj pvalue %s", ..2, round(..4, 2)),
            base_size = 6, rel_heights = c(1.5, 0.5), subplots = c(1, 2)
        )  
    ))
# Visualise the first plot
de_all_gene_rank_plots %>% 
    pull(plot) %>% 
    wrap_plots()

