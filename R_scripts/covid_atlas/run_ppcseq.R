# read {cell_name}_DE.rds to run ppcseq function
library(tidyverse)
library(glue)
library(ppcseq)
# install.packages("rlang")
# devtools::install_github("stemangiola/ppcseq@dev")

args = commandArgs(trailingOnly = TRUE)
filename <- args[1]
out_file <- args[2]

# command to run ppcseq 
input_file = readRDS(filename)
input_file <- input_file %>%
    identify_outliers(
        formula = ~ 0 + disease_severity_standard + filename + sex_standard,
        .sample = sample_name,
        .transcript = transcript,
        .abundance = abundance_RNA,
        .significance = .significance,
        .do_check = is_significant,
        percent_false_positive_genes = 5,
        pass_fit = TRUE) %>% 
    mutate(sample_wise_data = map(sample_wise_data, ~{attr(.x, "fit") = NULL; .x})) %>% 
    saveRDS(file = out_file)
