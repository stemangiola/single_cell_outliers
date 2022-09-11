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



ASDC_tidy <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/tidy_data/ASDC_tidy.rds")
# Factor the disease_severity_standard
counts <- ASDC_tidy %>% # 94 samples
# counts <- B_memory_tidy %>% # total 245 distinct samples
    mutate(disease_severity_standard = case_when(
        disease_severity_standard == "mild" ~ "mild",
        disease_severity_standard == "moderate" ~ "moderate",
        disease_severity_standard == "severe" ~ "severe",
        TRUE ~ "healthy"
        ))

# For solving duplicated sample/gene pair in multiple dataset:
counts <- counts %>% 
    unite("sample_name", c(sample, filename), sep = "_", remove = FALSE)

counts <- counts %>% keep_abundant(
    .sample = sample_name,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = disease_severity_standard) %>%
    scale_abundance(.sample = sample_name,
                    .transcript = transcript,
                    .abundance = abundance_RNA) 

# Generate rectangular dataset
counts <- counts %>%
    # nest(data = -transcript) %>%
    nest(data = -c(transcript, disease_severity_standard)) %>%
    add_count(transcript) %>%
    # pull(n) %>%
    # table()
    filter(n == n_distinct(disease_severity_standard)) %>%
    select(-n) %>%
    unnest(data) %>%
    impute_missing_abundance(~ disease_severity_standard,
                             .sample = sample_name,
                             .transcript = transcript,
                             .abundance = abundance_RNA)

# Differential expression analysis
counts = counts %>%
    # reduce dimension
    reduce_dimensions(.element = sample_name,
                                .feature = transcript,
                                .abundance = abundance_RNA,
        method = "PCA") %>%
    tidybulk::test_differential_abundance(~ 0 + disease_severity_standard + filename + sex_standard,
                                .sample = sample_name,
                                .transcript = transcript,
                                .abundance = abundance_RNA,
                                .contrasts = c(
                                        "(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy",
                                        "(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)",
                                        "disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)"
                                    )
                                ) 

counts %>% 
    saveRDS(file = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/ASDC_de.rds"))
# Run ppcseq

# # counts_test <- counts_test %>%
# counts_test_2 = counts_test_2 %>% 
#     # Group by contrast. Comparisons both ways.
#     pivot_longer(
#         cols = contains("___"),
#         names_to = c("stats", "contrast"),
#         values_to = ".value",
#         names_sep="___"
#     ) %>%
#     
#     # Markers selection within each pair of contrast
#     nest(stat_df = -contrast) %>%
#     
#     # Reshape inside each contrast
#     mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>% 
#     unnest(stat_df)
    
    
counts <- counts %>%
    mutate(abundance_RNA = abundance_RNA %>% as.integer()) %>% # run ppcseq
    # mutate(is_significant = !!sym(FDR_column_name) < 0.05)
    mutate(is_significant = `FDR___(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy`< 0.05 |
               `FDR___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)` < 0.05 |
               `FDR___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)`< 0.05)


library(ppcseq)
# devtools::install_github("stemangiola/ppcseq@dev")
counts <- counts %>%
    mutate(.significance = (`PValue___(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy` +
                                 `PValue___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)` +
                                 `PValue___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)`)/3) %>%
    # dplyr::select(-c(TMM, multiplier, .abundant, abundance_RNA_scaled)) %>% 
    identify_outliers(
        formula = ~ 0 + disease_severity_standard + filename + sex_standard,
        # formula = ~ disease_severity_standard,
        .sample = sample_name,
        .transcript = transcript,
        .abundance = abundance_RNA,
        # .significance = !!pvalue_column_name,
        .significance = .significance,
        .do_check = is_significant,
        percent_false_positive_genes = 5,
        how_many_negative_controls = Inf
    ) %>% 
    mutate(sample_wise_data = map(sample_wise_data, ~ {attr(.x, "fit") = NULL; .x})) %>%
    saveRDS(file = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/ASDC_ppcseq.rds"))

# test_create_summary.R
# read ppcseq data
counts_ppcseq <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/ASDC_ppcseq.rds")


outlier_data <- counts_ppcseq %>% left_join(
    counts %>% pivot_transcript(transcript) %>%
        arrange(.significance) %>% rowid_to_column(var = "rank")
)

outlier_summary_table <-
    tibble(cell_type = cell_name,
           num_genes_with_outliers = sum(outlier_data$tot_deleterious_outliers > 0),
           total_num_genes = nrow(outlier_data),
           num_genes_with_outliers_in_top_PValue = sum(
               outlier_data$tot_deleterious_outliers > 0 &
                   outlier_data$rank <= 10
           ), 
           method = !!method
    )
outlier_summary_table %>% saveRDS(file = output_file)





# counts %>% group_by(sample) %>% 
#     dplyr::summarise(count = n_distinct(transcript))
#     
# aggregate(data=counts, transcript ~ sample, function(x) length(unique(x)))

