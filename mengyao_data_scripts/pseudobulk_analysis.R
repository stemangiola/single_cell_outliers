# devtools::install_github("stemangiola/tidySingleCellExperiment")
library(tidySummarizedExperiment)
library(tidysc)
library(tidyverse)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(dittoSeq)
library(plotly)
library(ggrepel)
library(GGally)
library(tidybulk)

#----------------------------------------------------------------------------
# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )


setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data")
GSE120575_melanoma_bulk =
  readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/GSE120575_melanoma.rds") %>%
  tidysc::aggregate_cells(c(sample, cell_type), slot = "counts")  # 11 different cell types
'''
# A tibble: 11 x 1
cell_type                           
<fct>                               
  1 G5-Lymphocytes                      
2 G7-Regulatory T cells               
3 G9-Exhausted/HS CD8+ T cells        
4 G4-Dendritic cells                  
5 G10-Memory T cells                  
6 G1-B cells                          
7 G11-Lymphocytes exhausted/cell cycle
8 G3-Monocytes/Macrophages            
9 G8-Cytotoxicity Lymphocytes         
10 G6-Exhausted CD8+ T cells           
11 G2-Plasma cells    
'''



# TidyBulk 
# set up data
counts_G5 <-  
  GSE120575_melanoma_bulk %>% 
  filter(cell_type == "G5-Lymphocytes")

# Pre-processing data # factor of interest = time (pre/post)
counts_scaled_G5 <- counts_G5 %>% 
  keep_abundant(.sample = sample,
                .transcript = transcript,
                .abundance = abundance_RNA,
                factor_of_interest = condition) %>% 
  scale_abundance(.sample = sample,
                  .transcript = transcript,
                  .abundance = abundance_RNA)
  
  
counts_scaled_G5 %>%
  # Reshaping
  pivot_longer(cols = c("abundance_RNA", "abundance_RNA_scaled"), names_to = "source", values_to = "abundance") %>%
  
  # Plotting
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  custom_theme




# t-SNE 
counts_G5.tSNE =
  counts_scaled_G5 %>%
  identify_abundant() %>%
  reduce_dimensions(
    method = "tSNE",
    perplexity=10,
    pca_scale =TRUE
  )
counts_G5.tSNE %>%
  pivot_sample() %>%
  select(contains("tSNE"), everything()) 
# plot
counts_G5.tSNE %>%
  pivot_sample() %>%
  ggplot(aes(x = `tSNE1`, y = `tSNE2`, color = condition)) + geom_point() + custom_theme


'''
# Dimensional Reduction_ PCA
counts_scal_G5_PCA <-
  counts_scaled_G5 %>%
  reduce_dimensions(method = "PCA")

counts_scal_G5_PCA # into 2 new columns (PC1, PC2)

# plot sample-wise information by using pivot_sample
counts_scal_G5_PCA %>% pivot_sample()

# Plot sample with dimensional reduction:
# PCA plot-----------
counts_scal_G5_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = condition )) + 
  geom_point() +
  custom_theme

# Differential Expression: ~ condition + sample
de_all_G5_tidy <- counts_scal_G5_PCA %>%
  # edgeR QLT
  test_differential_abundance(.sample = sample,
                              .transcript = transcript,
                              .abundance = abundance_RNA,
    ~ condition ,
    method = "edgeR_quasi_likelihood",
    prefix = "edgerQLT_"
  ) 
#  
  # edgeR LRT
  test_differential_abundance(
    ~ dex + sample,
    method = "edgeR_likelihood_ratio",
    prefix = "edgerLR_"
  ) %>%
  
  # limma-voom
  test_differential_abundance(
    ~ dex + sample,
    method = "limma_voom",
    prefix = "voom_"
  ) %>%
  
  # DESeq2
  test_differential_abundance(
    ~ dex + sample,
    method = "deseq2",
    prefix = "deseq2_"
  )

counts_de_G5 =
  counts_scal_G5_PCA %>%
  identify_abundant(factor_of_interest = condition) %>%
  test_differential_abundance(
    ~ 0 + condition,                  
    .contrasts = c( "Non-responder - Responder"),
    action="get"
  )

'''

# Test abundance
counts_G5.de =
  counts_G5.tSNE %>%
  test_differential_abundance( ~ time)
counts_G5.de


### volcano plot, minimal
# run default edgeR method
counts_de_G5 <- counts_G5.tSNE %>%
  test_differential_abundance(~ time)
counts_de_G5 %>%
  pivot_transcript() %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point() +
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme



##### ppcseq

fileConn<-file("/home/users/allstaff/ma.m/single_cell_outliers/parse_data_scripts/pseudobulk_analysis.R")
writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
close(fileConn)

devtools::install_github("stemangiola/ppcseq")
library(ppcseq)

counts_de_G5

# Import libraries
counts_de_G5 <- counts_de_G5  %>% 
  mutate(abundance_RNA = abundance_RNA %>% as.integer)

counts.ppc = 
  counts_de_G5 %>%
  mutate(is_significant = FDR < 0.05) %>%
  identify_outliers(
    formula = ~ condition,
    .sample = sample, 
    .transcript = transcript,
    .abundance = abundance_RNA,
    .significance = PValue,
    .do_check = is_significant,
    percent_false_positive_genes = 5
  )

counts.ppc # contain plots for each gene

# visualise the top five differentially transcribed genes
counts.ppc_plots = 
  counts.ppc %>% 
  plot_credible_intervals() 

counts.ppc_plots %>%
  pull(plot) %>% 
  .[1:2]


