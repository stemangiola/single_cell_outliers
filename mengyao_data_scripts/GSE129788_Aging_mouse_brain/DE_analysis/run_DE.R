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

# command to run DE analysis and save to rds
counts = readRDS(in_file) 

# Filtering out lowly expressed counts
counts = counts %>% keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = condition
) %>% # factor of interest = condition (young/old)
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) 

# Differential expression analysis
counts = counts %>% 
    # reduce dimension
    reduce_dimensions(method = "PCA") %>% 
    test_differential_abundance(~ condition) %>% 
    saveRDS(glue("{out_file}"))







