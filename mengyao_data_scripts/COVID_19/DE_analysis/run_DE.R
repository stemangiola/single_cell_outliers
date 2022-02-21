# use makeflow to run DE
library(tidySummarizedExperiment)
library(tidysc)
library(tidyverse)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(dittoSeq)
library(plotly)
library(ggrepel)
library(GGally)
library(tidybulk)
library(glue)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

friendly_cols <- dittoSeq::dittoColors()
# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(
          t = 10,
          r = 10,
          b = 10,
          l = 10
        )),
        axis.title.y = element_text(margin = margin(
          t = 10,
          r = 10,
          b = 10,
          l = 10
        )),
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          vjust = 1
        )
      )
  )

# command to run DE analysis and save to rds
counts = readRDS(in_file) 

# Filtering out lowly expressed counts
counts = counts %>% keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = response
  ) %>% # factor of interest = response
  scale_abundance(.sample = sample,
                  .transcript = transcript,
                  .abundance = abundance_RNA) 

# Differential expression analysis
counts = counts %>% 
  # reduce dimension
  reduce_dimensions(method = "PCA") %>% 
  test_differential_abundance(~ response) %>% 
  saveRDS(glue("{out_file}"))
  
  
  
  
  
  
  
  
  

