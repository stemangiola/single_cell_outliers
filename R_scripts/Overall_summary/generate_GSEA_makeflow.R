# deseq2 :
# pval column name = pvalue  
# p.adj column name = padj

library(tidyverse)
library(glue)

script_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/Overall_summary/"
data_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/"

tidy_files = list.files(path = glue("{data_dir}all_raw_data/"))
cell_name = tools::file_path_sans_ext(tidy_files) %>% basename() %>% str_replace("_raw","")

all_command <- tidyr::expand_grid(
    GSEA_dir = glue("{data_dir}GSEA_data/"),
    input_dir = glue("{data_dir}all_raw_data/"), 
    Rscript = glue("Rscript {script_dir}perform_GSEA.R"), # including 3 methods
    files = glue("{cell_name}_raw.rds"),
    method = c("deseq2", "edgeR_quasi_likelihood","edger_robust_likelihood_ratio")) %>%
    
    mutate(cell_name = str_replace(files, "_raw.rds", "")) %>% 
    unite("file_name", c(method, cell_name),sep = "_", remove = FALSE) %>% 
    unite("GSEA_file", c(GSEA_dir, file_name), sep = "", remove = FALSE) %>% 
    mutate(GSEA_file = glue("{GSEA_file}_GSEA.rds")) %>% 
    unite("input_file", c(input_dir, files), sep = "", remove = FALSE) %>% 
    
    mutate(logFC_column_name = case_when(
        method == "deseq2" ~ "log2FoldChange",
        method == "edgeR_quasi_likelihood" ~ "logFC",
        method == "edger_robust_likelihood_ratio" ~ "logFC"
    )) %>% 
    mutate(pval_column_name = case_when(
        method == "deseq2" ~ "pvalue",
        method == "edgeR_quasi_likelihood" ~ "PValue",
        method == "edger_robust_likelihood_ratio" ~ "PValue"
    )) %>% 
    
    mutate(ppcseq_dir = glue("{data_dir}ppcseq_data/")) %>% 
    unite("ppcseq_file", c(ppcseq_dir, file_name), sep = "", remove = FALSE) %>% 
    mutate(ppcseq_file = glue("{ppcseq_file}_ppcseq.rds")) %>% 
    
    mutate(command_GSEA = sprintf("%s:%s\n\t%s %s %s %s %s %s %s %s",
                                  GSEA_file, input_file, Rscript, input_file, GSEA_file, method, 
                                  pval_column_name, logFC_column_name, ppcseq_file, cell_name
                                  ))
    
all_command <- all_command %>% pull(command_GSEA) %>% 
    prepend("CATEGORY=COVID19_GSEA\nMEMORY=30024\nCORES=4") %>% 
    write_lines(glue("{script_dir}GSEA.makeflow"))

