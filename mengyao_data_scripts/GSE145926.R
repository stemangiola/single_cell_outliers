## 	Single-cell landscape of bronchoalveolar immune cells in COVID-19 patients
## load data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926
Test 
setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE145926_BALF_COVID_19/")
library(tidyverse)
library(tibble)
library(dplyr)
library(Seurat)
# install.packages("hdf5r")
library("hdf5r")
library(readr)
library(tidyseurat)
# GSE145926_RAW <- read_csv("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE145926_BALF_COVID_19/data/GSE145926_RAW.tar")
# GSE145926_RAW


# ID C141
GSM4385990_C141_filtered_contig_annotations_csv <- 
  read_csv("./data/GSM4385990_C141_filtered_contig_annotations.csv.gz") %>%
  mutate(cell=barcode)
## Google _matrix .h5m
C141 = 
  Read10X_h5("./data/GSM4339769_C141_filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject() %>%
  left_join(GSM4385990_C141_filtered_contig_annotations_csv %>% nest(TCR = -cell))
  

###--------------------------------------
# ID C142
GSM4385991_C142_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4385991_C142_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C142 = 
  Read10X_h5("./data/GSM4339770_C142_filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject() %>%
  left_join(GSM4385991_C142_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C143
GSM4385992_C143_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4385992_C143_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C143 = 
  Read10X_h5("./data/GSM4339771_C143_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4385992_C143_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C144
GSM4385993_C144_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4385993_C144_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C144 = 
  Read10X_h5("./data/GSM4339772_C144_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4385993_C144_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C145
GSM4385994_C145_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4385994_C145_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell =barcode)

C145 = 
  Read10X_h5("./data/GSM4339773_C145_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4385994_C145_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C146
GSM4385995_C146_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4385995_C146_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C146 = 
  Read10X_h5("./data/GSM4339774_C146_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4385995_C146_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C148
GSM4475054_C148_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4475054_C148_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C148 = 
  Read10X_h5("./data/GSM4475051_C148_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4475054_C148_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C149
GSM4475055_C149_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4475055_C149_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C149 = 
  Read10X_h5("./data/GSM4475052_C149_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4475055_C149_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# ID C152
GSM4475056_C152_filtered_contig_annotations_csv <- 
  read_csv("data/GSM4475056_C152_filtered_contig_annotations.csv.gz") %>% 
  mutate(cell=barcode)

C152 = 
  Read10X_h5("./data/GSM4475053_C152_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() %>% 
  left_join(GSM4475056_C152_filtered_contig_annotations_csv %>% nest(TCR = -cell))

# C51
C51 = 
  Read10X_h5("./data/GSM4475048_C51_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject() 

# C52
C52 = 
  Read10X_h5("./data/GSM4475049_C52_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject()

# C100
C100 =
  Read10X_h5("./data/GSM4475050_C100_filtered_feature_bc_matrix.h5") %>% 
  CreateSeuratObject()

###--------------------------------------
# Left_join adds columns
# bind_rows adds rows
counts_nested = 
  tibble(sample = c("C141", "C142","C143","C144","C145","C146","C148","C149","C152","C51","C52")) %>%
  mutate(seurat = list(
    C141, C142, C143, C144, C145, C146, C148, C149, C152, C51, C52
  ))

annotation_nested =
  cell_annotations %>% 
  mutate(cell = gsub("_", "-", cell)) %>%
  nest(cell_type_array = cell_type)%>%
  nest(annotation = -sample) 

counts_nested %>%
  left_join(annotation_nested) %>%
  mutate(seurat = map2(seurat, annotation, ~ left_join(.x, .y, by="cell"))) %>%
  select(-annotation) %>%
  unnest_seurat(seurat)
  



# all_cell_annotation
all_cell_annotation_meta <- read_delim("data/all.cell.annotation.meta.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(sample = sample, cell_type = celltype, cell_cluster = cluster) %>% 
  mutate(dataset_id = "GSE145926")

all_cell_annotation_meta
distinct(all_cell_annotation_meta, cell_type)
# myeloid_cell_annotation
myeloid_cell_annotation_meta <- read_delim("data/myeloid.cell.annotation.meta.txt", 
                                           "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(sample = sample, cell_type = celltype, cell_cluster = cluster) %>% 
  mutate(dataset_id = "GSE145926")

myeloid_cell_annotation_meta
distinct(myeloid_cell_annotation_meta, cell_type)


# NKT cell annotation
NKT_cell_annotation_meta <- read_delim("data/NKT.cell.annotation.meta.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(sample = sample, cell_type = celltype, cell_cluster = cluster) %>% 
  mutate(dataset_id = "GSE145926")

NKT_cell_annotation_meta

#####

## 122,877 x 9
cell_annotations =
  all_cell_annotation_meta %>% 
  bind_rows(myeloid_cell_annotation_meta) %>% 
  bind_rows(NKT_cell_annotation_meta) %>% 
  rename(cell = ID)
  

cell_annotations %>% 
  left_join (counts)


saveRDS(GSE145926, "./data/GSE145926_BALF_COVID_19.rds")  

