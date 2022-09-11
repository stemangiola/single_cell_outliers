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
library(tidybulk)
library(glue)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]


counts = readRDS(in_file)
counts <- counts %>% 
    mutate(disease_severity_standard = case_when(
        disease_severity_standard == "mild" ~ "mild",
        disease_severity_standard == "moderate" ~ "moderate",
        disease_severity_standard == "severe" ~ "severe",
        TRUE ~ "healthy"
    ))

# Filtering out lowly expressed counts
counts <- counts %>%
    unite("sample_name", c(sample, filename), sep = "_", remove = FALSE)

counts <- counts %>% 
    tidybulk::keep_abundant(
    .sample = sample_name,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = disease_severity_standard) %>%
    tidybulk::scale_abundance(.sample = sample_name,
                    .transcript = transcript,
                    .abundance = abundance_RNA) 


# Generate rectangular dataset by imputation
# counts <- counts %>%
#     nest(data = -c(transcript, disease_severity_standard)) %>%
#     add_count(transcript) %>%
#     # pull(n) %>%
#     # table()
#     filter(n == n_distinct(disease_severity_standard)) %>%
#     select(-n) %>%
#     unnest(data) %>%
#     impute_missing_abundance(~ disease_severity_standard,
#                              .sample = sample_name,
#                              .transcript = transcript,
#                              .abundance = abundance_RNA)



# Differential expression analysis
counts = counts %>%
    # reduce dimension
    reduce_dimensions(
        .element = sample_name,
        .feature = transcript,
        .abundance = abundance_RNA,
        method = "PCA"
    ) %>%
    tidybulk::test_differential_abundance(
        ~ 0 + disease_severity_standard + filename + sex_standard,
        .sample = sample_name,
        .transcript = transcript,
        .abundance = abundance_RNA,
        contrasts = c(
            "(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy",
            "(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)",
            "disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)"
        )
    ) %>%
    mutate(abundance_RNA = abundance_RNA %>% as.integer()) %>%
    mutate(
        is_significant = `FDR___(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy` < 0.05 |
            `FDR___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)` < 0.05 |
            `FDR___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)` < 0.05
    ) %>% 

    mutate(
        .significance = (
            `PValue___(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy` +
                `PValue___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)` +
                `PValue___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)`
        ) / 3
    ) %>%
    saveRDS(glue("{out_file}"))



