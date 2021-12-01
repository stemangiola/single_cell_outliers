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

# total sample --
n_distinct(counts_COVID_19$sample)

# select one cell type
counts_41BB_Hi_CD8_T_cell <-
  counts %>%
  filter(cell_type == "41BB-Hi CD8+ T cell")
# Filtering out lowly expressed counts
counts_scaled_41BB_Hi_CD8_T_cell <- counts_41BB_Hi_CD8_T_cell %>%
  keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_originalexp,
    factor_of_interest = ICB_Exposed
  ) %>% # factor of interst = severity
  scale_abundance(.sample = sample,
                  .transcript = transcript,
                  .abundance = abundance_originalexp)

# visualize abundance difference before and after scaling
counts_scaled_41BB_Hi_CD8_T_cell %>%
  # Reshaping
  pivot_longer(
    cols = c("abundance_originalexp", "abundance_originalexp_scaled"),
    names_to = "source",
    values_to = "abundance"
  ) %>%
  # Plotting
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap( ~ source) +
  scale_x_log10() +
  custom_theme

# PCA dimensional reduction
counts_scaled_CD8_Tcell_PCA <-
  counts_scaled_41BB_Hi_CD8_T_cell %>%
  reduce_dimensions(method = "PCA")

# # A tibble: 2 x 2
# `Fraction of variance`    PC
# <dbl> <int>
#   1                  0.695     1
# 2                  0.223     2
# tidybulk says: to access the raw results do `attr(..., "internals")$PCA`

counts_scaled_CD8_Tcell_PCA # into 2 new columns (PC1, PC2)

# plot sample-wise information by using pivot_sample
# counts_scaled_B_cell_PCA %>% pivot_sample()

# Plot sample with dimensional reduction:
# PCA plot-----------
counts_scaled_CD8_Tcell_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = ICB_Exposed)) +
  geom_point() +
  custom_theme

# DE
counts_B_cell_de =
  counts_scaled_B_cell_PCA %>%
  test_differential_abundance(~ severity)  # ~ factor_of_interest + unesplicable_cluster_separation
counts_B_cell_de
# Volcano plot
counts_B_cell_de %>%
  pivot_transcript() %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point() +
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme

# ppcseq
##### ppcseq

fileConn <-
  file(
    "/home/users/allstaff/ma.m/single_cell_outliers/parse_data_scripts/pseudobulk_analysis_COVID19.R"
  )
writeLines(
  c(
    "CXX14FLAGS += -O3",
    "CXX14FLAGS += -DSTAN_THREADS",
    "CXX14FLAGS += -pthread"
  ),
  fileConn
)
close(fileConn)

#devtools::install_github("stemangiola/ppcseq")
library(ppcseq)

counts_B_cell_de

# Import libraries
counts_B_cell_de <- counts_B_cell_de  %>%
  mutate(abundance_RNA = abundance_RNA %>% as.integer)

counts.ppc =
  counts_B_cell_de %>%
  mutate(is_significant = FDR < 0.05) %>%
  identify_outliers(
    formula = ~ severity,
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    .significance = PValue,
    .do_check = is_significant,
    percent_false_positive_genes = 5
  )

counts.ppc # contain plots for each gene

# visualise the top five differentially transcribed genes
counts.ppc_plots =
  counts.ppc %>%
  plot_credible_intervals()

counts.ppc_plots %>%
  pull(plot) %>%
  .[1:2]
