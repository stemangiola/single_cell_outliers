# read {cell_name}_DE.rds to run ppcseq function
library(tidyverse)
library(glue)
library(ppcseq)
# install.packages("rlang")
# devtools::install_github("stemangiola/ppcseq@dev")

args = commandArgs(trailingOnly = TRUE)
filename <- args[1]
out_file <- args[2]
method <- args[3]
FDR_column_name <- args[4]
pvalue_column_name <- args[5]

# command to run ppcseq 
input_file = readRDS(filename)
input_file = input_file %>%
    # edgeR_quasi_likelihood_Lymphocytes_exhausted_cell_cycle_DE %>% 
    mutate(abundance_RNA = abundance_RNA %>% as.integer) %>% # run ppcseq
    mutate(is_significant = !!sym(FDR_column_name) < 0.05)
    # mutate(is_significant = edgerQLT_FDR < 0.05) # edgerQLT_FDR
# filter(FDR_coloumn_name %>% is.na()) %>%
# filter(edgerQLT_FDR %>% is.na()) %>% 

input_file %>%
    identify_outliers(
        formula = ~ condition,
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA,
        # .significance = edgerQLT_PValue,
        .significance = !!pvalue_column_name,
        # .significance = edgerQLT_PValue,
        .do_check = is_significant,
        percent_false_positive_genes = 5
    ) %>%
    mutate(sample_wise_data = map(sample_wise_data, ~ {attr(.x, "fit") = NULL; .x})) %>%
    mutate(method = method) %>%
    saveRDS(file = out_file)
