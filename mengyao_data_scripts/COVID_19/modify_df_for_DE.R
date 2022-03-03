# Pesudobulk analysis for COVID_19 dataset
# cited in https://www.nature.com/articles/s41587-020-0602-4
### Latest version 
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
library(ppcseq)
library(glue)

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/")
counts_COVID_19 =
  readRDS("./s41587-020-0602-4_COVID_19.rds") %>%
  tidysc::aggregate_cells(c(sample, cell_type), slot = "counts") # 24 distinct cell types

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
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

# Rename the factor of interest column 
counts_COVID_19 <- counts_COVID_19 %>% 
  mutate(response = fct_recode(severity, 
                               "Not Critical" = "control",
                               "Not Critical" = "moderate",
                               "Critical" = "critical"))

# DE analysis 
counts_COVID_19 <- counts_COVID_19 %>%
  nest(data_cell_type = -cell_type) %>%
  # slice(1) %>%
  mutate(data_cell_type = map(
    data_cell_type,
    ~ .x %>% keep_abundant(
      .sample = sample,
      .transcript = transcript,
      .abundance = abundance_RNA,
      factor_of_interest = response
    ) %>%
      scale_abundance(
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA
      )
  )) 
# Density plot
counts_COVID_19 <- counts_COVID_19 %>%
  mutate(
    density_plot = map2(
      data_cell_type, cell_type,
      ~ .x %>% pivot_longer(
        cols = c("abundance_RNA", "abundance_RNA_scaled"),
        names_to = "source",
        values_to = "abundance"
      ) %>%
        # Plotting
        # ggsave(glue(
        #   "Density plot of {.y}"
        # ), path = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/plots/"),
      ggplot(aes(x = abundance + 1, color = sample)) +
        ggtitle(glue(
          "Density plot of {.y}"
        )) +
        geom_density() +
        facet_wrap(~ source) +
        scale_x_log10() +
        custom_theme
    )
  ) %>%
  
# %>% map2(
#     density_plot,
#     cell_type,
#     ~ .x %>% pull() %>% ggsave(glue("Density plot of {.y}"),
#                                path = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/plots/")
#   )

#   # pull(plot_density) %>% .[[1]] # to see the first row 
#   #-------save for 24 types 
  # PCA plot
  mutate(data_cell_type = map(data_cell_type,
                              ~ .x %>% reduce_dimensions(method = "PCA"))) %>% # PCA dimensional reduction
  mutate(
    PCA_plot = map2(
      data_cell_type, cell_type,
      ~ .x %>% pivot_sample() %>%
        ggplot(aes(
          x = PC1, y = PC2, colour = response)) + ggtitle(glue("PCA plot of {.y}")) +
        geom_point() +
        custom_theme
    )
  ) %>% # pull(PCA_plot) %>% .[[3]]
#   #--------save PCA plots for all 24 cell types
#   
  mutate(data_cell_type = map(
    data_cell_type,
    ~ .x %>% test_differential_abundance(~ response)
  )) %>%
#
  mutate(
    volcano_plot = map2(
      data_cell_type, cell_type,
      ~ .x %>% pivot_transcript() %>%
        ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
        ggtitle(glue("Volcano plot of {.y}")) +
        geom_point() +
        scale_y_continuous(trans = "log10_reverse") +
        custom_theme
    )
  ) 
 
counts_COVID_19 %>% saveRDS("./All_DE_results.rds")


for (i in 1:nrow(counts_COVID_19)) {
  ggsave(glue("Density plot of {counts_COVID_19$cell_type[i]}.pdf"),
         plot = counts_COVID_19$density_plot[[i]],
         path = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/plots/Density_plot/",
         width=297, height=210,units = "mm",dpi=300)
}

for (i in 1:nrow(counts_COVID_19)) {
  ggsave(glue("PCA plot of {counts_COVID_19$cell_type[i]}.pdf"),
         plot = counts_COVID_19$PCA_plot[[i]],
         path = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/plots/PCA_plot/",
         width=297, height=210,units = "mm",dpi=300)
}

for (i in 1:nrow(counts_COVID_19)) {
  ggsave(glue("Volcano of {counts_COVID_19$cell_type[i]}.pdf"),
         plot = counts_COVID_19$volcano_plot[[i]],
         path = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/plots/Volcano_plot/",
         width=297, height=210,units = "mm",dpi=300)
}

  # # ppcseq
  # mutate(data_cell_type = map(
  #   data_cell_type,
  #   ~ .x %>% mutate(abundance_RNA = abundance_RNA %>% as.integer)
  # )) %>%
  # mutate(data_cell_type = map(
  #   data_cell_type,
  #   ~ .x %>% mutate(is_significant = FDR < 0.05) %>%
  #     identify_outliers(
  #       formula = ~ response,
  #       .sample = sample,
  #       .transcript = transcript,
  #       .abundance = abundance_RNA,
  #       .significance = PValue,
  #       .do_check = is_significant,
  #       percent_false_positive_genes = 5
  #     )
  # )) %>%
  # 
  # # ppcseq plot
  # mutate(ppcseq_plot = map(
  #   data_cell_type,
  #   ~ .x %>% plot_credible_intervals()
  # )) %>%
  # pull(ppcseq_plot) %>%
  # .[1:2]
