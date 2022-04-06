# create summary table
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
input_ppcseq <- args[1]
input_de <- args[2]
output_file <- args[3]

ppcseq_data <- readRDS(input_ppcseq)
de_data <- readRDS(input_de)
cell_name = de_data$cell_type[1]

if (nrow(ppcseq_data) != 0) {
    outlier_data <- ppcseq_data %>% left_join(
        de_data %>% pivot_transcript(transcript) %>%
            arrange(PValue) %>% rowid_to_column(var = "rank")
    )
    
    outlier_summary_table <-
        tibble(cell_type = cell_name,
               num_genes_with_outliers = sum(outlier_data$tot_deleterious_outliers > 0),
               total_num_genes = nrow(outlier_data),
               num_genes_with_outliers_in_top_PValue = sum(
                   outlier_data$tot_deleterious_outliers > 0 &
                       outlier_data$rank <= 10
               )
        )
    outlier_summary_table %>% saveRDS(file = output_file)
} else{
    my_tibble = tibble(cell_type = character(),
                       num_genes_with_outliers = integer(),
                       total_num_genes = integer(),
                       num_genes_with_outliers_in_top_PValue = integer()) %>% 
        saveRDS(file = output_file)
}