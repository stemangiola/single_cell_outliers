# save pseudo_bulk_data for each cell_type

# load packages
library(tidySummarizedExperiment)
library(tidysc)
library(tidyverse)
library(seurat)
library(tidyseurat)
library(tidybulk)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(dittoSeq)
library(plotly)
library(GGally)
library(glue)

#----------------------------------------------------------------------------
# load data
in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/tidy_data/"

df <-  readRDS(glue("{in_dir}SCP1288_renal_cell_carcinoma.rds")) %>% 
    tidysc::aggregate_cells(c(sample, cell_type), slot = "counts")  # 34 different cell types

# Modify 
dataframe <- df %>% 
    mutate(abundance=abundance_originalexp) 
dataframe$cell_type <- gsub(" ", "_", dataframe$cell_type)
dataframe$cell_type <- gsub("/", "_", dataframe$cell_type)

# remove NA values

# A tibble: 3 × 2
# ICB_Exposed       n
# <chr>         <int>
#     1 ICB         9275931
# 2 NoICB       5638311
# 3 NA            60627

dataframe <- dataframe %>% filter(ICB_Exposed != "NA")

dataframe %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/SCP1288_renal_cell_carcinoma_tidyseruat.rds")

dataframe = dataframe %>% nest(data_cell_type = -cell_type) %>% 
    mutate(saved = map2(
        data_cell_type, cell_type,
        ~ .x %>%
            mutate(cell_type = .y) %>%
            saveRDS(glue("{out_dir}{.y}_tidy.rds"))
    ))



