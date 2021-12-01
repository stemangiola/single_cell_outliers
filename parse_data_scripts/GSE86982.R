## 	REGION-SPECIFIC NEURAL STEM CELL LINEAGES REVEALED BY SINGLE-CELL RNA-SEQ FROM HUMAN EMBRYONIC STEM CELLS
library(tidyverse)
library(dplyr)
library(tidybulk)

setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers/GSE86982_hESC")
library(readr)
GSE86982_smartseq_tpm_csv <- read_csv("data/GSE86982_smartseq.tpm.csv.gz")
GSE86982_smartseq_tpm_csv
## unzip----------------

GSE86982_smartseq_tpm_csv %>%
  tidybulk::as_matrix(rownames = "X1") %>%
  CreateSeuratObject() %>%
  saveRDS("GSE86982_hESC/data/GSE86982.rds") 
