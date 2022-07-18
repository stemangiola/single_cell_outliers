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
library(glue)

#----------------------------------------------------------------------------
# load data
in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE130973_aging_human_skin/data/raw_data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE130973_aging_human_skin/data/tidy_data/"

df = readRDS(glue("{in_dir}GSE130973_seurat.rds")) %>% 
    tidysc::aggregate_cells(c(sample, cell_type), slot = "counts")  # 25 different cell types

df = df %>% nest(data_cell_type = -cell_type) %>% 
    mutate(saved = map2(
        data_cell_type, cell_type,
        ~ .x %>%
            mutate(cell_type = .y) %>%
            saveRDS(glue("{out_dir}{.y}_tidy.rds"))
    ))
