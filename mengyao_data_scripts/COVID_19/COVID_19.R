### DA-seq on COVID-19 data
### Original paper: https://www.nature.com/articles/s41587-020-0602-4
### This script reproduces analysis presented in Figure 4
### orginal article:https://www.nature.com/articles/s41587-020-0602-4


library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyseurat)
# source("convenience.R")

## Set Python and GPU
#python2use <- "/data/henry/henry_env/venv/bin/python"
#GPU <- 3

## Set path for FIt-SNE R wrapper
#fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


##=============================================##
## Data

## Load data
main_S <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/22927382.1")

# main_S <- read.delim("./data/22927382.1")
DefaultAssay(main_S) <- "integrated"
head(main_S@meta.data)
table(main_S@meta.data$severity)

# check commanes
names(main_S@commands)

distinct(main_S,sample)

main_S %>% 
  rename(cell_type = celltype, sample = sample) %>%
  mutate(dataset_id = "s41587-020-0602-4") %>% 
  saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/s41587-020-0602-4_COVID_19.rds")
## Prepare object according to commands

# analysis
#main_S <- ScaleData(main_S)
#main_S <- RunPCA(
#  main_S, npcs = 90, verbose = F
#)
#main_S <- runFItSNE(
 # main_S, dims.use = 1:90, seed.use = 3, fast.R.path = fitsneR,
  #ann_not_vptree = FALSE, nthreads = 12
#)
#TSNEPlot(main_S, group.by = "severity")
#TSNEPlot(main_S, group.by = "celltype", label = T)



## Separate immune cells from patient samples

# get immune cells
'''
sort(unique(main_S@meta.data$celltype))
epi_type <- c("Basal","Ciliated","Ciliated-diff","FOXN4","Ionocyte","IRC",
              "Secretory","Secretory-diff","Squamous","outliers_epithelial","unknown_epithelial")
immune_type <- setdiff(sort(unique(main_S@meta.data$celltype)), epi_type)

immune_S <- subset(x = main_S, cells = which(main_S@meta.data$celltype %in% immune_type))

# remove control cells
immune_S <- subset(immune_S, cells = which(immune_S@meta.data$severity != "control"))
immune_S@meta.data$severity <- factor(as.character(immune_S@meta.data$severity))

TSNEPlot(immune_S, group.by = "severity")
TSNEPlot(immune_S, group.by = "celltype", label = T)


# assign a number label for each cell type
immune_S@meta.data$celltype_num <- as.numeric(factor(immune_S@meta.data$celltype))
'''