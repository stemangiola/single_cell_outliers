### data:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130973
### Original paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7181753/
library(tidyverse)
library(tidybulk)
library(SingleCellExperiment)
library(Seurat) #V3
library(tidyseurat)

## Load data
dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE130973_aging_human_skin/data/raw_data/"
counts <- Read10X(data.dir = dir)
GSE130973_seurat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
GSE130973_seurat <- GSE130973_seurat %>% 
    mutate(dataset_id = "GSE130973") %>% 
    mutate(sample = orig.ident) %>% 
    saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE130973_aging_human_skin/data/raw_data/GSE130973_seurat.rds")

