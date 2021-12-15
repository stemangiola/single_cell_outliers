# Pesudobulk analysis for COVID_19 dataset
# cited in https://www.nature.com/articles/s41587-020-0602-4
# install.packages("job")
library(job)
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
  readRDS("./s41587-020-0602-4_COVID_19.rds") %>%
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
# save_rds for differential expression results
job::job({counts_de_Squamous %>% saveRDS("./Squamous_de_result_7_12.rds", compress = "xz")})


# Volcano plot -- label out the significantly expressed genes 
counts_de_Squamous %>%
  pivot_transcript() %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point() +
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme

# topgenes_symbols <-
counts_de_Squamous %>%
  filter(FDR < 0.05) %>%
  distinct(transcript)  # 1043 differentiated transcripts


# ppcseq
# dir.create(file.path("~/", ".R"), showWarnings = FALSE)
# fileConn<-file("~/.R/Makevars")
# writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
# close(fileConn)


# Import libraries
library(ppcseq)

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


# job::job({counts_ppc_Squamous %>% saveRDS("./Squamous_ppcseq_7_12.rds")})

# counts_ppc_Squamous %>%
#   saveRDS(
#     "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/Squamous_ppcseq_result_06_12_21.rds"
#   )

# counts_ppc_Squamous %>%
#   count(how_many_genes_include_outliers = `tot deleterious outliers` > 0)

# # A tibble: 2 × 2
# how_many_genes_include_outliers     n
# <lgl>                           <int>
#   1 FALSE                             761
# 2 TRUE                              282

# Summary tibble
outlier_analysis_Squamous <- counts_ppc_Squamous %>%
  left_join(
    counts_de_Squamous %>% pivot_transcript(transcript) %>% arrange(PValue) %>% rowid_to_column(var =
                                                                                                  "rank")
  )

outlier_analysis_Squamous %>% saveRDS("./outlier_analysis_Squamous_9_12.rds")


# summary tibble
# num_genes_with_outliers | total_num_genes | num_genes_with_outliers_in_top_10_PValue

outlier_summary_table <-
  tibble(
    num_genes_with_outliers = sum(outlier_analysis_Squamous$`tot deleterious outliers` > 0),
    total_num_genes = nrow(outlier_analysis_Squamous),
    num_genes_with_outliers_in_top_PValue = sum(
      outlier_analysis_Squamous$`tot deleterious outliers` > 0 &
        outlier_analysis_Squamous$rank <= 10
    )
  )

outlier_summary_table %>% saveRDS("./outlier_summary_table.rds")


# visualise the top five differential transcribed genes
counts_ppc_Squamous_plots <- 
  counts_ppc_Squamous %>%
  plot_credible_intervals()

counts_ppc_Squamous_plots %>%
  pull(plot) %>% 
  .[1:2]
