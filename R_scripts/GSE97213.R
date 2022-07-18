## 	Single Cell and Open Chromatin Analysis Reveals Molecular Origin of Epidermal Cells of the Skin
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97213


library(tidyverse)
library(dplyr)
library(tidybulk)
library(Seurat)
setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers/GSE97213_Mouse_skin_epidermal_cell")

## load data
library(readr)
GSE97213_Mouse_skin_epidermal_cell <- read_delim("data/GSE97213_RAW.tar", 
                                                       delim = "\t", escape_double = FALSE, 
                                                       trim_ws = TRUE)

GSE97213 = 
  GSE97213_Mouse_skin_epidermal_cell[!duplicated(GSE72056_melanoma_single_cell_revised_v2$Cell),] %>%
  
  tidybulk:::as_matrix(rownames = "Cell") %>%
  CreateSeuratObject(counts = ., min.cells = 3, min.genes = 200, project = "GSE97213.R")

GSE72056 %>% saveRDS("project/GSE72056/GSE72056.rds")