# Execute `GSEA` using `MSigDB`, `clusterProfiler` and `enrichplot` 
# import COVID_19 dataset all_de.rds
# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
# tidyverse-friendly packages
library(plotly)
library(ggrepel)
library(GGally)
library(tidyverse)
library(tidybulk)
library(tidySummarizedExperiment) # we'll load this below to show what it can do
library(patchwork)
library(org.Hs.eg.db)


# all_de <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/all_de.rds")

# edgeR_quasi_likelihood_Secretory_DE <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/de_data/edgeR_quasi_likelihood_Secretory_DE.rds")

# edgeR_quasi_likelihood_Secretory_DE_consider_outliers = edgeR_quasi_likelihood_Secretory_DE
Secretory <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/all_raw_data/Secretory_raw.rds")
Secretory_exclude_outliers <- Secretory

###-----------------------------------------------------------------------------------------------------------------------------------------
Secretory_ppcseq <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/ppcseq_data/edgeR_quasi_likelihood_Secretory_ppcseq.rds")

Secretory_ppcseq  <- Secretory_ppcseq %>% 
    filter(tot_deleterious_outliers > 0) %>% pull(transcript)

Secretory_exclude_outliers <- Secretory_exclude_outliers %>% 
    filter(!transcript %in% (Secretory_ppcseq)) %>% 
    impute_missing_abundance(~ response,
                            .sample = sample,
                            .transcript = transcript,
                            .abundance = abundance_RNA)

Secretory_exclude_outliers_gene_rank <- Secretory_exclude_outliers %>% 
    # Annotate
    mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys = transcript,
                                          keytype = "SYMBOL",
                                          column = "ENTREZID",
                                          multiVals = "first"
    )) %>% 
    keep_abundant(
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA,
        factor_of_interest = response) %>%  # factor of interest = response
    
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) %>% 
    
    reduce_dimensions(method = "PCA") %>% 
    # Test differential gene transcript abundance
    test_differential_abundance(~ response) %>% 
    
    filter(PValue   %>% is.na %>% `!`) %>%
    # filter(pvalue%>% is.na %>% `!`) %>%
    test_gene_rank(
        .sample = sample,
        .entrez = entrez,
        .arrange_desc = logFC ,
        species="Homo sapiens",
        gene_sets = c("H", "C2", "C5")
    )


# Examine significantly enriched gene sets
# de_all_gene_rank_list1 <- de_all_gene_rank %>%
Secretory_exclude_outliers_gene_rank_list <- Secretory_exclude_outliers_gene_rank %>%
    filter(gs_cat  == "C2" ) %>%
    dplyr::select(-fit) %>%
    unnest(test) %>% 
    filter(p.adjust < 0.05)
    # filter(padj < 0.05)

a <- Secretory_exclude_outliers_gene_rank_list %>% arrange(p.adjust) %>% pull(ID) %>% head(10)
b <- Secretory_gene_rank_list %>% arrange(p.adjust) %>% pull(ID) %>% head(10)

# Generate fraction: 
c <- tibble(overlapped_prop_in_top_10_pathways = c(sum(as.numeric(a %in% b))/10),
            cell_type = c("Secretory"),
            method = "edgeQLT")
c %>% 
    ggplot(aes(x = method, y = overlapped_prop_in_top_10_pathways)) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 0.1)
# To compare with the dataset including the outliers

Secretory_gene_rank <- Secretory %>% 
    # Annotate
    
    mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys = transcript,
                                          keytype = "SYMBOL",
                                          column = "ENTREZID",
                                          multiVals = "first"
    )) %>% 
    keep_abundant(
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA,
        factor_of_interest = response) %>%  # factor of interest = response
    
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) %>% 
    
    reduce_dimensions(method = "PCA") %>% 
    # Test differential gene transcript abundance
    test_differential_abundance(~ response) %>% 
    
    filter(PValue   %>% is.na %>% `!`) %>%
    test_gene_rank(
        .sample = sample,
        .entrez = entrez,
        .arrange_desc = logFC ,
        species="Homo sapiens",
        gene_sets = c("H", "C2", "C5")
    )


Secretory_gene_rank_list <- Secretory_gene_rank %>%
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

