### DA-seq on aging mouse brain data
### Original paper: https://www.nature.com/articles/s41593-019-0491-3

library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

# Please go to https://singlecell.broadinstitute.org/single_cell/study/SCP263/aging-mouse-brain#/ to download
# processed data.

## Load data
# exp mat
data_exp <- read.table(
  "./expression_Aging_mouse_brain_portal_data_updated.txt", sep = "\t", header = T, 
  row.names = 1, stringsAsFactors = F
)

data_exp <- Matrix(as.matrix(data_exp), sparse = T)


# meta data
data_meta <- read.table(
  "./meta_Aging_mouse_brain_portal_data.txt", sep = "\t", header = T, stringsAsFactors = F
)
data_meta <- data_meta[-1,]
data_meta$cell_type <- gsub("NEUR_mature","mNEUR",data_meta$cell_type)
data_meta$cell_type <- gsub("NEUR_immature","ImmN",data_meta$cell_type)

table(data_meta$cell_type)
table(data_meta$all_cells_by_age)
rownames(data_meta) <- data_meta[,1]

celltype2label <- c(
    "1-OPC","2-OLG","3-OEG","4-NSC","5-ARP","6-ASC","7-NRP","8-ImmN","9-mNEUR","10-NendC",
    "11-EPC","12-HypEPC","13-TNC","14-CPC","15-EC","16-PC","17-VSMC","18-Hb_VC","19-VLMC","20-ABC",
    "21-MG","22-MNC","23-MAC","24-DC","25-NEUT"
)
names(celltype2label) <- sapply(celltype2label, function(x) unlist(strsplit(x,"-"))[2])

## Seurat
# create object
data_S <- CreateSeuratObject(
  counts = data_exp, names.field = 6, names.delim = "_"
)
table(data_S@meta.data$orig.ident)


# add metadata
data_S@meta.data$cell_type <- data_meta[colnames(data_S),"cell_type"]
data_S@meta.data$age <- data_meta[colnames(data_S),"all_cells_by_age"]
data_S@meta.data$condition <- data_S@meta.data$age
data_S@meta.data$condition[data_S@meta.data$age == "2-3mo"] <- "young"
data_S@meta.data$condition[data_S@meta.data$age == "21-22mo"] <- "old"
table(data_S@meta.data$condition)

data_S@meta.data$cell_type_label <- factor(celltype2label[data_S@meta.data$cell_type], levels = celltype2label)
data_S@meta.data$cell_type_num <- as.numeric(data_S@meta.data$cell_type_label)

# change cell_type name
data_S = 
    data_S %>%
    mutate(cell_type =case_when(
        cell_type == "OPC" ~ "Oligodendrocyte_precursor_cells",
        cell_type == "OLG" ~ "Oligodendrocytes",
        cell_type == "OEG" ~ "Olfactory_ensheathing_glia",
        cell_type == "NSC" ~ "Neural_stem_cells",
        cell_type == "ARP" ~ "Astrocyte_restricted_precursors",
        cell_type == "ASC" ~ "Astrocytes",
        cell_type == "NRP" ~ "Neuronal_restricted_precursors",
        cell_type == "ImmN" ~ "NEUR_immature",
        cell_type == "mNEUR" ~ "NEUR_mature",
        cell_type == "NendC" ~ "Neuroendocrine_cells",
        cell_type == "EPC" ~ "Ependymocytes",
        cell_type == "HypEPC" ~ "Hypendymal_cells",
        cell_type == "TNC" ~ "Tanycytes",
        cell_type == "CPC" ~ "Choroid_plexus_epithelial_cells",
        cell_type == "EC" ~ "Endothelial_cells",
        cell_type == "PC" ~ "Pericytes",
        cell_type == "VSMC" ~ "Vascular_smooth_muscle_cells",
        cell_type == "Hb_VC" ~ "Hemoglobin_expressing_vascular_cells",
        cell_type == "VLMC" ~ "Vascular_and_leptomeningeal_cells",
        cell_type == "ABC" ~ "Arachnoid_barrier_cells",
        cell_type == "MG" ~ "Microglia",
        cell_type == "MNC" ~ "Monocytes",
        cell_type == "MAC" ~ "Macrophages",
        cell_type == "DC" ~ "Dendritic_cells",
        cell_type == "NEUT" ~ "Neutrophils"
    ))

# final formate the dataset
data_S <- data_S %>% 
    mutate(dataset_id = "GSE129788") %>% 
    mutate(sample = orig.ident) %>% 
    saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE129788_aging_mouse_brain/data/GSE129788_aging_mouse_brain.rds")

# analysis
#data_S <- ScaleData(data_S)
#data_S <- FindVariableFeatures(
#  data_S, selection.method = "mvp", mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf)
#)
#data_S <- RunPCA(
#  data_S, npcs = 20, verbose = F
#)
#data_S <- runFItSNE(
#  data_S, dims.use = 1:20, seed.use = 3, fast.R.path = fitsneR, 
#  ann_not_vptree = FALSE, nthreads = 12
#)

#TSNEPlot(data_S, group.by = "age", pt.size = 0.1)
#TSNEPlot(data_S, group.by = "cell_type_label", label = T, pt.size = 0.1)