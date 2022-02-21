# Makeflow execute ppcseq 
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
library(ppcseq)
library(glue)
library(ppcseq)
input_address = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/"
# read DE results
input_de = readRDS(glue("{input_address}All_DE_results.rds"))
input_de %>% slice(1) %>% 
  mutate(data_cell_type = map(
  data_cell_type,
  ~ .x %>% mutate(abundance_RNA = abundance_RNA %>% as.integer)
)) %>%   
  # execute ppcseq for all cell types (new column: ppcseq_result)
  mutate(ppcseq_result = map(
    data_cell_type,
    ~ .x %>% mutate(is_significant = FDR < 0.05) %>%
      identify_outliers(
        formula = ~ response,
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA,
        .significance = PValue,
        .do_check = is_significant,
        percent_false_positive_genes = 5
      )
  )) %>% 
  # generate ppcseq plots
  # mutate(ppcseq_plot = map(
  #   ppcseq_result, 
  #   ~ .x %>% plot_credible_intervals()
  # )) %>% 
  # pull(ppcseq_plot) %>%
  # .[1:2] %>% 

  # save ppcseq results
  mutate(ppcseq_result = map(
    ppcseq_result, ~.x %>% mutate(`sample wise data` = map(`sample wise data`, ~ {
      attr(.x, "fit") = NULL
      .x
    }))
  )) %>%
  select(cell_type, ppcseq_result) %>% saveRDS(glue("{input_address}All_ppcseq.rds"))

# Generate rank
  # mutate(rank_table = map2(
  #   ppcseq_result, data_cell_type,
  #   ~ .x %>% left_join({.y} %>%
  #                        pivot_transcript(transcript) %>% arrange(PValue) %>%
  #                        rowid_to_column(var ="rank"))
  # )) %>%
  # 
  # # Create summary table
  # mutate(summary_table = map(
  #   rank_table, ~.x %>% tibble(
  #     num_genes_with_outliers = sum({.x}$`tot deleterious outliers` > 0),
  #     total_num_genes = nrow({.x}),
  #     num_genes_with_outliers_in_top_PValue = sum(
  #       {.x}$`tot deleterious outliers` > 0 &
  #         {.x}$rank <= 10
  #     )
  #   )
  # )) %>%
  # 
  # # select cell_type and summary table , unnest
  # select(cell_type,ppcseq_plot,summary_table) %>% saveRDS(glue("{input_address}tot_table.rds"))
  
  