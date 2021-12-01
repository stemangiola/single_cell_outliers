# Transcriptomic Profiling of the Developing Cardiac Conduction System at Single-Cell Resolution
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132658

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE132658_Mouse_heart/")

library(Seurat) #V3
library(Matrix)
library(reshape2) 
library(ggplot2)
library(cowplot)
library(tidyseurat)
library(readr)
library(tibble)
library(SingleCellExperiment)
library(Matrix)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DropletUtils")
library(DropletUtils)
# Load data
untar("./data/GSE132658_RAW.tar") 

# Read in `matrix.mtx`
counts_SAN <- readMM("data/filtered_gene_bc_matrices/hg19/matrix.mtx")

# Read in `genes.tsv`
genes_SAN <- read_tsv("data/filtered_gene_bc_matrices/hg19/genes.tsv", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cell_ids_SAN <- read_tsv("data/filtered_gene_bc_matrices/hg19/barcodes.tsv", col_names = FALSE)$X1

