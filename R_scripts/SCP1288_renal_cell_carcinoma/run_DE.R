# Analyze SCP1288_renal_cell_carcinoma
# use makeflow to run DE
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
library(glue)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]
method = args[3]



prefix = glue("{args[3]}")
prefix = case_when(prefix == "edgeR_quasi_likelihood" ~ "edgerQLT_",
                   prefix == "edger_robust_likelihood_ratio" ~ "edgerRobust_",
                   prefix == "deseq2" ~ "deseq2_")


# command to run DE analysis and save to rds
counts = readRDS(in_file) 

# Filtering out lowly expressed counts
counts = counts %>% keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance,
    factor_of_interest = ICB_Exposed
) %>% # factor of interest =  ICB_Exposed(ICB/NoICB)
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance) 

# Differential expression analysis
counts = counts %>% 
    # reduce dimension
    reduce_dimensions(method = "PCA") %>% 
    test_differential_abundance(~ ICB_Exposed,
                                method = method,
                                prefix = prefix) %>% 
    mutate(method = method) %>% 
    saveRDS(glue("{out_file}"))







