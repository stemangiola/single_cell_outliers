# Pesudobulk analysis for COVID_19 dataset
# cited in https://www.nature.com/articles/s41587-020-0602-4

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
#----------------------------------------------------------------------------
# Use colourblind-friendly colours
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

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/")
counts_COVID_19 =
  readRDS(
    "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/s41587-020-0602-4_COVID_19.rds"
  ) %>%
  tidysc::aggregate_cells(c(sample, cell_type), slot = "counts") # 24 distinct cell types

#---------------------------------------------------------------------------------------------------------
counts_COVID_19 %>% distinct(severity) # critical, moderate, control
counts_COVID_19 %>% count(sample) # total samples 

# how many samples in each group
counts_COVID_19 %>% 
  group_by(severity) %>% 
  summarise(n_sample = n_distinct(sample),.groups = "drop")

# A tibble: 3 × 2
# severity n_sample
# <fct>       <int>
# 1 control         5
# 2 moderate       14
# 3 critical       13

# total sample 
n_distinct(counts_COVID_19$sample) # 32

# create a new column (response) to rename the factor of interests 
counts_COVID_19 <- counts_COVID_19 %>% 
  mutate(response = fct_recode(severity, 
                              "Not Critical" = "control",
                              "Not Critical" = "moderate",
                              "Critical" = "critical"))
  
##--------------------------------------------------------------------------------
# select one cell type:Squamous
counts_Squamous <-
  counts_COVID_19 %>%
  filter(cell_type == "Squamous")

# Filtering out lowly expressed counts
counts_scaled_Squamous <- counts_Squamous %>%
  keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = response
  ) %>% # factor of interest = response
  scale_abundance(.sample = sample,
                  .transcript = transcript,
                  .abundance = abundance_RNA)

# visualize abundance difference before and after scaling -- density plot
counts_scaled_Squamous %>%
  # Reshaping
  pivot_longer(
    cols = c("abundance_RNA", "abundance_RNA_scaled"),
    names_to = "source",
    values_to = "abundance"
  ) %>%
  # Plotting
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~ source) +
  scale_x_log10() +
  custom_theme 
# ggsave("Squamous_density_plot.pdf", path = "/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers/plots/",width=6, height=4,dpi=300)

# PCA dimensional reduction
counts_scaled_Squamous_PCA <-
  counts_scaled_Squamous %>%
  reduce_dimensions(method = "PCA")

counts_scaled_Squamous_PCA # into 2 new columns (PC1, PC2)

# # A tibble: 2 × 2
# `Fraction of variance`    PC
# <dbl> <int>
#   1                  0.318     1
# 2                  0.143     2
# tidybulk says: to access the raw results do `attr(..., "internals")$PCA`

# PCA plot
counts_scaled_Squamous_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = response)) +
  geom_point() +
  custom_theme

# Differential expression analysis -----------------------------------------------------------------
counts_de_Squamous <- 
  counts_scaled_Squamous_PCA %>%
  test_differential_abundance(~ response)  # ~ factor_of_interest + unsplicable_cluster_separation
counts_de_Squamous


# Volcano plot -- label out the significantly expressed genes 
counts_de_Squamous %>%
  pivot_transcript() %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point() +
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme

# ppcseq
fileConn<-file("/home/users/allstaff/ma.m/single_cell_outliers/mengyao_data_scripts/pseudobulk_analysis_COVID_19.R")
writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
close(fileConn)

# Import libraries
library(ppcseq)

counts_de_Squamous
# Import libraries
counts_de_Squamous <- counts_de_Squamous  %>%
  mutate(abundance_RNA = abundance_RNA %>% as.integer) # changed into integer value

counts_ppc_Squamous <-
  counts_de_Squamous %>%
  mutate(is_significant = FDR < 0.05) %>%
  identify_outliers(
    formula = ~ response,
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    .significance = PValue,
    .do_check = is_significant,
    percent_false_positive_genes = 5
  )

counts_ppc_Squamous # contain plots for each gene

# visualise the top five differential transcribed genes
counts_ppc_Squamous_plots <- 
  counts_ppc_Squamous %>%
  plot_credible_intervals()


counts_ppc_Squamous_plots %>%
  pull(plot) %>% 
  .[1:2]
