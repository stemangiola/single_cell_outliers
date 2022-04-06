# generate DE makefile.R
library(tidyverse)
library(glue)

script_dir = "/stornext/HPCScratch/home/ma.m/mengyao_data_scripts/SCP1288_renal_cell_carcinoma/DE_analysis/"
input_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/filtered_data/"
output_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/de_data/"


files = list.files(path = glue("{input_dir}"))
cell_name = tools::file_path_sans_ext(files) %>% basename() %>% str_replace("_filtered","")

command_df = tibble(output_dir = output_dir, 
                    cell_name = cell_name,
                    Rscript = glue("Rscript {script_dir}run_DE.R"),
                    input_dir = input_dir,
                    files = files) %>% 
    unite("output_file", c(output_dir, cell_name), sep = "") %>% 
    mutate(output_file = glue("{output_file}_DE.rds")) %>% 
    unite("args", c(input_dir, files), sep = "") %>% 
    mutate(command = sprintf("%s:\n\t%s %s %s",
                             output_file, Rscript, args, output_file))

# pull command and add SLURM requirements
command_df %>% pull(command) %>% 
    prepend("CATEGORY=run_DE_analysis\nMEMORY=30024\nCORES=4") %>% 
    write_lines(glue("{script_dir}run_de.makeflow"))