## 	Single-cell RNA-seq reveals dynamic paracrine control of cellular variation
## mouse dendritic cells (DCs) stimulated with three pathogenic components
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48968

library(tidyverse)
library(dplyr)
library(tidybulk)

setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers/GSE48968_mouse_dendritic cells")
## load data 
## 1,700 primary mouse dendritic cells
library(readr)
GSE48968_allgenesTPM_GSM1406531_GSM1407094_txt <- read_delim("data/GSE48968_allgenesTPM_GSM1406531_GSM1407094.txt.gz", 
                                                             "\t", escape_double = FALSE, trim_ws = TRUE)
GSE48968_allgenesTPM_GSM1406531_GSM1407094_txt %>% 
CreateSeuratObject()


### Stimulated data
GSE48968_allgenesTPM_GSM1189042_GSM1190902_txt <- read_delim("data/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz", 
                                                             "\t", escape_double = FALSE, trim_ws = TRUE)
GSE48968_allgenesTPM_GSM1189042_GSM1190902_txt

saveRDS(GSE48968_allgenesTPM_GSM1189042_GSM1190902_txt,"data/GSE48968_GSM1190902.rds") 
saveRDS(GSE48968_allgenesTPM_GSM1406531_GSM1407094_txt,"data/GSE48968_GSM1406531.rds") 
