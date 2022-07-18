### download file
### description: Single Cell and Open Chromatin Analysis Reveals Molecular Origin of Epidermal Cells of the Skin
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102086 
### orginal paper: https://www.cell.com/developmental-cell/fulltext/S1534-5807(18)30680-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1534580718306804%3Fshowall%3Dtrue

library(tidyverse)
library(tidySeurat) 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("DropletUtils")
getwd()
setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers/GSE102086_mouse_epidermis_cell")

## load data
data_gene <- read.table("./data/GSE102086_genes.tsv.gz")
data_gene
data_exp <- read.table("./data/GSE102086_RAW.tar",header=FALSE, sep="",na.strings = "NA",fill=TRUE, quote='',
                       row.names=,skipNul = TRUE)

data_S <- CreateSeuratObject(
  counts = data_exp, project = "Epidermal Cells" # project name
)

saveRDS(data_S, "/data/seurat_GSE102086.rds")
