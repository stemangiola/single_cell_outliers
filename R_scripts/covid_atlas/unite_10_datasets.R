library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ppcseq)
library(plyr)
library(tidySummarizedExperiment)
 
data_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/pseudobulk_data/"
out_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/tidy_data/"
#  
files <- list.files(path = glue("{data_dir}"))
covid_atlas <- tibble(filename = files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(glue("{data_dir}"), .)) %>%
                                   as_tibble()) # a new data column
    )

covid_atlas$filename <- gsub(".rds", "", covid_atlas$filename)
covid_atlas$filename <- gsub("-","_", covid_atlas$filename)

# unnest the file_contents
covid_atlas <- unnest(covid_atlas)


# Check the distinct cell types
covid_atlas$predicted.celltype.l2 %>% unique()

covid_atlas <- covid_atlas %>% mutate(cell_type = predicted.celltype.l2)
covid_atlas$cell_type <- gsub(" ", "_", covid_atlas$cell_type)
covid_atlas %>% saveRDS(glue("{data_dir}united_10_datasets.rds"))

# covid_atlas <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/pseudobulk_data/united_10_datasets.rds")

# total 7308 samples
# covid_atlas %>% dplyr::count(transcript) %>% dplyr::count(n) %>% arrange(desc(n))
covid_atlas <- covid_atlas %>% add_count(transcript) %>% filter(n == 7604)

covid_atlas = covid_atlas %>% nest(data_cell_type = -cell_type) %>% 
    mutate(saved = map2(
        data_cell_type, cell_type,
        ~ .x %>%
            mutate(cell_type = .y) %>%
            saveRDS(glue("{out_dir}{.y}_tidy.rds"))
    ))
