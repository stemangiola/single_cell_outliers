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


