GSE120575_melanoma_bulk %>%
  nest(data_cell_type = -cell_type) %>%
  slice(1) %>%
  mutate(data_cell_type = map(
    data_cell_type,
    ~ .x %>%  keep_abundant(
      .sample = sample,
      .transcript = transcript,
      .abundance = abundance_RNA,
      factor_of_interest = condition
    ) %>%
      scale_abundance(
        .sample = sample,
        .transcript = transcript,
        .abundance = abundance_RNA
      )
  )) %>%
  mutate(
    plot_density = map(
      data_cell_type,
      ~ .x %>% pivot_longer(
        cols = c("abundance_RNA", "abundance_RNA_scaled"),
        names_to = "source",
        values_to = "abundance"
      ) %>%
        
        # Plotting
        ggplot(aes(x = abundance + 1, color = sample)) +
        geom_density() +
        facet_wrap( ~ source) +
        scale_x_log10() +
        custom_theme
    )
  ) %>%
  pull(plot_density) %>% .[[1]]
