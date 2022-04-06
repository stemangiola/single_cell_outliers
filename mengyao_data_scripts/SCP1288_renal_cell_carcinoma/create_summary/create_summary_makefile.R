# create_summary_makefile.R
library(tidyverse)
library(glue)

script_dir = "/stornext/HPCScratch/home/ma.m/mengyao_data_scripts/SCP1288_renal_cell_carcinoma/create_summary/"
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/ppcseq_data/"
de_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/de_data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/summary_table/"


ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
de_files = list.files(path = glue("{de_dir}"))
cell_name = tools::file_path_sans_ext(de_files) %>% basename() %>% str_replace_all("_DE","")

command_df = tibble(out_dir = out_dir, 
                    cell_name = cell_name,
                    Rscript = glue("Rscript {script_dir}create_summary.R"),
                    ppcseq_dir = ppcseq_dir,
                    de_dir = de_dir,
                    ppcseq_files = ppcseq_files,
                    de_files = de_files) %>% 
    unite("output_file", c(out_dir, cell_name), sep = "") %>% 
    mutate(output_file = glue("{output_file}_table.rds")) %>% 
    unite("arg1", c(ppcseq_dir, ppcseq_files), sep = "") %>% 
    unite("arg2", c(de_dir, de_files), sep = "") %>% 
    mutate(command = sprintf("%s:\n\t%s %s %s %s",
                             output_file, Rscript, arg1, arg2, output_file))
# pull command and add SLURM requirements
command_df %>% pull(command) %>% 
    prepend("CATEGORY=create_summary_table\nMEMORY=30024\nCORES=4") %>% 
    write_lines(glue("{script_dir}create_summary.makeflow"))
