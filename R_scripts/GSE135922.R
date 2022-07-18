## 	Single-cell transcriptomics of the human retinal pigment epithelium and 
## choroid in health and macular degeneration
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135922 

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE135922_Human_Retina")

library(Seurat) #V3
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyseurat)

###--------------------------------------------------------------------------
# Load data
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135922&format=file", 
             # destfile = "data/GSE135922_RAW.tar")
untar("./data/GSE135922_RAW.tar") 

library(readr)
GSM4037981_macula_donor_1_expression_tsv <- read_table2("data/GSM4037981_macula_donor_1_expression.tsv.gz")
GSM4037981_macula_donor_1_expression_tsv
