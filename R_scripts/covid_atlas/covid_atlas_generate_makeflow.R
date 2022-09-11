library(tidyverse)
library(glue)

script_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/covid_atlas/"
data_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/"

# Create result directory. Google how to create a directory if directory exists.
# GOOGLE_THE_FUNCTION("DIR", warning=FALSE)

tidy_files = list.files(path = glue("{data_dir}tidy_data/"))
# cell_name = tools::file_path_sans_ext(tidy_files) %>% basename() %>% str_replace("_tidy","")

all_command <- tibble(
    input_dir = glue("{data_dir}tidy_data/"),
    de_dir = glue("{data_dir}de_data/"),
    Rscript_de = glue("Rscript {script_dir}run_DE.R"),
    files = tidy_files) %>%
    
    mutate(cell_name = str_replace(files, "_tidy.rds", "")) %>% 
    unite("de_file", c(de_dir, cell_name), sep = "", remove = FALSE) %>% 
    mutate(de_file = glue("{de_file}_DE.rds")) %>% 
    unite("input_file", c(input_dir, files), sep = "", remove = FALSE) %>% 
    mutate(command_de = sprintf("%s:%s\n\t%s %s %s",
                                de_file,input_file, Rscript_de, input_file, de_file)) %>% 

    mutate(ppcseq_dir = glue("{data_dir}ppcseq_data/")) %>% 
    unite("ppcseq_file", c(ppcseq_dir, cell_name), sep = "", remove = FALSE) %>% 
    mutate(ppcseq_file = glue("{ppcseq_file}_ppcseq.rds")) %>% 
    mutate(Rscript_ppcseq = glue("Rscript {script_dir}run_ppcseq.R")) %>%
    mutate(command_ppcseq = sprintf("%s:%s\n\t%s %s %s",
                                    ppcseq_file, de_file, Rscript_ppcseq, de_file, ppcseq_file)) %>% 
    
    # generate summary table command
    mutate(summary_dir = glue("{data_dir}summary_table/")) %>% 
    unite("summary_file", c(summary_dir, cell_name), sep = "", remove = FALSE) %>% 
    mutate(summary_file = glue("{summary_file}_table.rds")) %>%
    unite("two_files", c(ppcseq_file,de_file), sep = " ", remove = FALSE) %>%
    mutate(Rscript_summary = glue("Rscript {script_dir}create_summary.R")) %>% 
    mutate(command_summary = sprintf("%s:%s\n\t%s %s %s %s",
                                     summary_file, two_files, Rscript_summary,ppcseq_file, de_file, summary_file)) %>% 
    pivot_longer(c(command_de, command_ppcseq, command_summary), values_to = "command") %>% 
    select(command)

all_command %>% pull(command) %>% 
    prepend("CATEGORY=covid_atlas\nMEMORY=30024\nCORES=4") %>% 
    write_lines(glue("{script_dir}covid_atlas.makeflow"))