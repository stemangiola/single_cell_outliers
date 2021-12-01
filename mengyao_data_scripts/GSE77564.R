## Coupled electrophysiological recording and 
# single-cell transcriptome analyses revealed molecular mechanisms underlying neuronal maturation Organism	
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77564
## article: https://link.springer.com/article/10.1007%2Fs13238-016-0247-8 

library(tidyverse)
library(dplyr)
library(tidybulk)
library(Seurat)
setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers")

library(readr)
GSE77564_single_neuron_exp_txt <- read_csv("GSE77564_single_neuron/data/GSE77564_single_neuron_exp.txt.gz")
GSE77564_single_neuron_exp_txt

# figure out cell_1 : cell type
'''
GSE86982_smartseq_tpm_csv %>%
  tidybulk::as_matrix(rownames = "X1") %>%
  CreateSeuratObject() %>%
  saveRDS("GSE77564_single_neuron//data/GSE77564.rds") 
saveRDS(GSE77564_single_neuron_exp_txt,"GSE77564_single_neuron/data/GSE77564.rds")
'''

GSE77564_single_neuron_exp_txt %>% 
  tidybulk::as_matrix(rownames = "X1") %>% 
  CreateSeuratObject() %>% 
  saveRDS("GSE77564_single_neuron/data/GSE77564.rds")

# This study compares this condition (3 samples) with this condition (2 samples) across 20 cells 
# human embryonic stem (ES) cell cultures & human induced pluripotent stem cell (hiPSC)
# -- Just indicated 20 neuron cells


# 20000 genes 4000 cells
# 1000 cells are t cells
# 4 cells are Mac M1


# 200 are cluster 1
# Cluster 2

# cell | annotation | cluster | cell_type 
# cell_1 | sample_1 | 2.       | erythrocyte

# Many cells
# Many samples
