# 1. tidy
# 2. filtered
# 3. de.rds
# 4. ppcseq
# 5. summary 

# input tidy.rds --> filtered.rds
# generate fiter_sample_makefile.R
library(tidyverse)
library(glue)

script_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/GSE120575_melanoma/"
data_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/GSE120575_melanoma/data/"

# Create result directory. Google how to create a directory if directory exists.
# GOOGLE_THE_FUNCTION("DIR", warning=FALSE)

tidy_files = list.files(path = glue("{data_dir}tidy_data/"))
cell_name = tools::file_path_sans_ext(tidy_files) %>% basename() %>% str_replace("_tidy","")

command_filter = tibble(output_dir = glue("{data_dir}filtered_data/"), 
                        cell_name = cell_name,
                        Rscript = glue("Rscript {script_dir}filter_low_abundance_samples.R"),
                        input_dir = glue("{data_dir}tidy_data/"),
                        files = tidy_files) %>% 
    unite("output_file", c(output_dir, cell_name), sep = "") %>%
    mutate(output_file = glue("{output_file}_filtered.rds")) %>%
    unite("input_file", c(input_dir, files), sep = "") %>% 
    # unite("args", c(input_dir, files), sep = "", remove = TRUE) %>% 
    mutate(command = sprintf("%s:%s\n\t%s %s %s",
                             output_file,input_file, Rscript, input_file, output_file)) %>% 
    select(command)

# DE analysis ------------------------------------------------------------------------------
# filtered_files as input
command <- tidyr::expand_grid(
    de_dir = glue("{data_dir}de_data/"),
    input_dir = glue("{data_dir}filtered_data/"),
    Rscript_de = glue("Rscript {script_dir}run_DE.R"), # including 3 methods
    files = glue("{cell_name}_filtered.rds"),
    method = c("deseq2", "edgeR_quasi_likelihood","edger_robust_likelihood_ratio")) %>%
    
    mutate(cell_name = str_replace(files, "_filtered.rds", "")) %>% 
    unite("file_name", c(method, cell_name),sep = "_", remove = FALSE) %>% 
    unite("de_file", c(de_dir, file_name), sep = "", remove = FALSE) %>% 
    mutate(de_file = glue("{de_file}_DE.rds")) %>% 
    unite("input_file", c(input_dir, files), sep = "", remove = FALSE) %>% 
    mutate(command_de = sprintf("%s:%s\n\t%s %s %s %s",
                                de_file,input_file, Rscript_de, input_file, de_file, method)) %>% 
    
    # generate running ppcseq command
    mutate(FDR_coloumn_name = case_when(
        method == "deseq2" ~ "deseq2_padj",
        method == "edgeR_quasi_likelihood" ~ "edgerQLT_FDR",
        method == "edger_robust_likelihood_ratio" ~ "edgerRobust_FDR"
    )) %>% 
    mutate(pvalue_column_name = case_when(
        method == "deseq2" ~ "deseq2_pvalue",
        method == "edgeR_quasi_likelihood" ~ "edgerQLT_PValue",
        method == "edger_robust_likelihood_ratio" ~ "edgerRobust_PValue"
    )) %>% 

    mutate(ppcseq_dir = glue("{data_dir}ppcseq_data/")) %>% 
    unite("ppcseq_file", c(ppcseq_dir, file_name), sep = "", remove = FALSE) %>% 
    mutate(ppcseq_file = glue("{ppcseq_file}_ppcseq.rds")) %>% 
    # mutate(Rscript_ppcseq = case_when(method == "deseq2" ~ glue("Rscript {script_dir}run_ppcseq_for_deseq2.R"),
    #                                   TRUE ~ glue("Rscript {script_dir}run_ppcseq.R"))) %>% 
    mutate(Rscript_ppcseq = glue("Rscript {script_dir}run_ppcseq.R")) %>%
    mutate(command_ppcseq = sprintf("%s:%s\n\t%s %s %s %s %s %s",
                             ppcseq_file, de_file, Rscript_ppcseq, de_file, ppcseq_file, method, FDR_coloumn_name, pvalue_column_name)) %>% 
        
    # generate summary table command
        mutate(summary_dir = glue("{data_dir}summary_table/")) %>% 
        unite("summary_file", c(summary_dir, file_name), sep = "", remove = FALSE) %>% 
        mutate(summary_file = glue("{summary_file}_table.rds")) %>%
        unite("two_files", c(ppcseq_file,de_file), sep = " ", remove = FALSE) %>%
        mutate(Rscript_summary = glue("Rscript {script_dir}create_summary.R")) %>% 
        mutate(command_summary = sprintf("%s:%s\n\t%s %s %s %s %s %s",
                                 summary_file, two_files, Rscript_summary,ppcseq_file, de_file, summary_file, method, pvalue_column_name)) %>% 
    pivot_longer(c(command_de, command_ppcseq, command_summary), values_to = "command") %>% 
    select(command)
        
 
all_command <- command_filter %>% rbind(command) 
all_command %>% pull(command) %>% 
            prepend("CATEGORY=GSE120575\nMEMORY=30024\nCORES=4") %>% 
            write_lines(glue("{script_dir}GSE120575.makeflow"))
        
    


# ppcseq -----------------------------------------------------------------------------------------------
# de_files = list.files(path = glue("{data_dir}de_data/"))
# For edgeQLT & edgeLT
# command_ppcseq_edgeR <- tidyr::expand_grid(
#     output_dir = glue("{data_dir}ppcseq_data/"),
#     input_dir = glue("{data_dir}de_data/"),
#     Rscript = glue("Rscript {script_dir}run_ppcseq.R"),
#     cell_name = cell_name,
#     method = c("edgeR_likelihood_ratio","edgeR_quasi_likelihood")) %>% 
#     unite("file_name", c(method, cell_name),sep = "_", remove = FALSE) %>% 
#     unite("output_file", c(output_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(output_file = glue("{output_file}_ppcseq.rds")) %>% 
#     unite("input_file", c(input_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(input_file = glue("{input_file}_DE.rds")) %>% 
#     mutate(command = sprintf("%s:%s\n\t%s %s %s %s",
#                              output_file,input_file, Rscript, input_file, output_file, method)) %>% 
#     select(command)
# 
# command_ppcseq_deseq2 <- tidyr::expand_grid(
#     output_dir = glue("{data_dir}ppcseq_data/"),
#     input_dir = glue("{data_dir}de_data/"),
#     Rscript = glue("Rscript {script_dir}run_ppcseq_for_deseq2.R"),
#     cell_name = cell_name,
#     method = c("deseq2")) %>% 
#     unite("file_name", c(method, cell_name),sep = "_", remove = FALSE) %>% 
#     unite("output_file", c(output_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(output_file = glue("{output_file}_ppcseq.rds")) %>% 
#     unite("input_file", c(input_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(input_file = glue("{input_file}_DE.rds")) %>% 
#     mutate(command = sprintf("%s:%s\n\t%s %s %s %s",
#                              output_file,input_file, Rscript, input_file, output_file, method)) %>% 
#     select(command)

# Summary table-------------------------------------------------------------------------------------
# ppcseq_files = list.files(path = glue("{data_dir}ppcseq_data/"))
# command_summary_table <- tidyr::expand_grid(
#     output_dir = glue("{data_dir}summary_table/"),
#     ppcseq_dir = glue("{data_dir}ppcseq_data/"),
#     de_dir = glue("{data_dir}de_data/"),
#     Rscript = glue("Rscript {script_dir}create_summary.R"),
#     cell_name = cell_name,
#     method = c("deseq2", "edgeR_likelihood_ratio","edgeR_quasi_likelihood")) %>% 
#     unite("file_name", c(method, cell_name),sep = "_", remove = FALSE) %>% 
#     unite("output_file", c(output_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(output_file = glue("{output_file}_table.rds")) %>% 
#     unite("de_file", c(de_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(de_file = glue("{de_file}_DE.rds")) %>%
#     unite("ppcseq_file", c(ppcseq_dir, file_name), sep = "", remove = FALSE) %>% 
#     mutate(ppcseq_file = glue("{ppcseq_file}_ppcseq.rds")) %>% 
#     unite("input_file", c(ppcseq_file, de_file), sep = " ", remove = FALSE) %>% 
#     mutate(command = sprintf("%s:%s\n\t%s %s %s %s",
#                              output_file, input_file, Rscript, ppcseq_file, de_file, output_file)) %>%
#     select(command)


