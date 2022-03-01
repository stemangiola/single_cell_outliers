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
in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/raw_data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/tidy_data/"

df = readRDS(glue("{in_dir}GSE120575_melanoma.rds")) %>% 
    tidysc::aggregate_cells(c(sample, cell_type), slot = "counts")  # 11 different cell types

# # A tibble: 11 x 1
# cell_type                           
# <fct>                               
# 1 G5-Lymphocytes                      
# 2 G7-Regulatory T cells               
# 3 G9-Exhausted/HS CD8+ T cells        
# 4 G4-Dendritic cells                  
# 5 G10-Memory T cells                  
# 6 G1-B cells                          
# 7 G11-Lymphocytes exhausted/cell cycle
# 8 G3-Monocytes/Macrophages            
# 9 G8-Cytotoxicity Lymphocytes         
# 10 G6-Exhausted CD8+ T cells           
# 11 G2-Plasma cells    

df = df %>% separate(cell_type, c(NA, "cell_type"), sep = "-")
df$cell_type <- gsub(" ", "_", df$cell_type)
df$cell_type <- gsub("/", "_", df$cell_type)

# A tibble: 11 × 1
# cell_type                       
# <chr>                           
#     1 Lymphocytes                     
# 2 Regulatory_T_cells              
# 3 Exhausted_HS_CD8+_T_cells       
# 4 Dendritic_cells                 
# 5 Memory_T_cells                  
# 6 B_cells                         
# 7 Lymphocytes_exhausted_cell_cycle
# 8 Monocytes_Macrophages           
# 9 Cytotoxicity_Lymphocytes        
# 10 Exhausted_CD8+_T_cells          
# 11 Plasma_cells         

df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/raw_data/GSE120575_raw.rds")

df = df %>% nest(data_cell_type = -cell_type) %>% 
    mutate(saved = map2(
        data_cell_type, cell_type,
        ~ .x %>%
            mutate(cell_type = .y) %>%
            saveRDS(glue("{out_dir}{.y}_tidy.rds"))
    ))
