library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(plotly)
library(ggrepel)
library(GGally)
library(tidyverse)
library(tidybulk)
library(tidySummarizedExperiment)
library(patchwork)
library(org.Hs.eg.db)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]
method = args[3] # input method for DE analysis
pval_column_name = args[4]
logFC_column_name = args[5]
ppcseq_file <- args[6]
cell_type = args[7]

counts = readRDS(in_file) 
counts_exclude_outliers = counts # replicate the counts
# Filtering out lowly expressed counts
counts_gene_rank = counts %>% 
    mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys = transcript,
                                          keytype = "SYMBOL",
                                          column = "ENTREZID",
                                          multiVals = "first"
    )) %>% keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = response) %>% # factor of interest = response
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) %>% 
    
    reduce_dimensions(method = "PCA") %>% 
    # Test differential gene transcript abundance
    test_differential_abundance(~ response,
                                method = method) %>% 
    
    # filter(PValue   %>% is.na %>% `!`) %>%
    filter(!!pval_column_name %>% is.na %>% `!`) %>%
    test_gene_rank(
        .sample = sample,
        .entrez = entrez,
        # .arrange_desc = logFC ,
        .arrange_desc = !!sym(logFC_column_name),
        species="Homo sapiens",
        gene_sets = c("H", "C2", "C5")
    )

counts_gene_rank_list <- counts_gene_rank %>%
    filter(gs_cat  == "C2" ) %>%
    dplyr::select(-fit) %>%
    unnest(test) %>%
    filter(p.adjust < 0.05)
#
# # Filter out the outlier transcripts and perform GSEA
# # Read ppcseq rds first and filter out the outlier transcripts
#
ppcseq_file <- readRDS(ppcseq_file)
# if (nrow(ppcseq_file) == 0) {
#     summary_tibble <- tibble(
#         overlapped_prop_in_top_10_pathways = 1,
#         cell_type = cell_type,
#         method = method) %>% saveRDS(file = out_file)
# } else {
#     
    ppcseq_file <- ppcseq_file %>%
        filter(tot_deleterious_outliers > 0) %>% pull(transcript)
    #
    counts_exclude_outliers <- counts_exclude_outliers %>%
        filter(!transcript %in% (ppcseq_file)) %>%
        impute_missing_abundance(
            ~ response,
            .sample = sample,
            .transcript = transcript,
            .abundance = abundance_RNA
        )
    
    # # Run GSEA
    counts_exclude_outliers_gene_rank <- counts_exclude_outliers %>%
        # Annotate
        mutate(
            entrez = AnnotationDbi::mapIds(
                org.Hs.eg.db::org.Hs.eg.db,
                keys = transcript,
                keytype = "SYMBOL",
                column = "ENTREZID",
                multiVals = "first"
            )
        ) %>%
        keep_abundant(
            .sample = sample,
            .transcript = transcript,
            .abundance = abundance_RNA,
            factor_of_interest = response
        ) %>%  # factor of interest = response
        
        scale_abundance(.sample = sample,
                        .transcript = transcript,
                        .abundance = abundance_RNA) %>%
        
        reduce_dimensions(method = "PCA") %>%
        # Test differential gene transcript abundance
        test_differential_abundance( ~ response,
                                     method = method) %>%
        
        filter(!!pval_column_name %>% is.na %>% `!`) %>%
        test_gene_rank(
            .sample = sample,
            .entrez = entrez,
            # .arrange_desc = logFC ,
            .arrange_desc = !!sym(logFC_column_name),
            species = "Homo sapiens",
            gene_sets = c("H", "C2", "C5")
        )
    
    counts_exclude_outliers_gene_rank_list <-
        counts_exclude_outliers_gene_rank %>%
        filter(gs_cat  == "C2") %>%
        dplyr::select(-fit) %>%
        unnest(test) %>%
        filter(p.adjust < 0.05)
    
    #
    # Compute the proportions of the overlapped pathways
    list_including_outliers <-
        counts_gene_rank_list %>% arrange(p.adjust) %>% pull(ID) %>% 
    list_excluding_outliers <-
        counts_exclude_outliers_gene_rank_list %>% arrange(p.adjust) %>% pull(ID) %>% head(10)
    
    # Generate fraction:
    summary_tibble <- tibble(
        overlapped_prop_in_top_10_pathways = c(sum(
            as.numeric(list_including_outliers %in% list_excluding_outliers)
        ) / 10),
        cell_type = cell_type,
        method = method)
    summary_tibble %>% saveRDS(file = out_file) 
    # }