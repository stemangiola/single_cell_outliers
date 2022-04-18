# read {cell_name}_DE.rds to run ppcseq function
library(tidyverse)
library(glue)
library(ppcseq)
# devtools::install_github("stemangiola/ppcseq@dev")

args = commandArgs(trailingOnly = TRUE)
filename <- args[1]
out_file <- args[2]

# command to run ppcseq 
input_file = readRDS(filename)
input_file <- input_file %>%
    mutate(abundance_RNA = abundance_RNA %>% as.integer) %>% # run ppcseq
    mutate(is_significant = FDR < 0.05) 

input_file %>%
    identify_outliers(
        formula = ~ condition,
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA,
        .significance = PValue,
        .do_check = is_significant,
        percent_false_positive_genes = 5
    ) %>%
    mutate(sample_wise_data = map(sample_wise_data, ~ {attr(.x, "fit") = NULL; .x})) %>%
    saveRDS(file = out_file)
