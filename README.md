# single_cell_outliers
## Data directories: 

```R
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/GSE120575_melanoma/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/GSE129788_aging_mouse_brain/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/SCP1288_renal_cell_carcinoma/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/Overall_summary/
```

The name of these 6 directories indicate the specific dataset ID and the studied aim. Except the directory, **Overall_summary**, the rest 5 directories two sub-directories: *data and plot*. The **plot** directory includes all the plots generated during my process including PCA plots, volcano plots and several summary bar plots.The name of the plot are specifically nominated which should be easily understood. For the directory **Overall_summary**, it only contains the main diagrams for the paper and several GSEA plot. The names of these plots are specifically labelled.  

The data directory include the raw data that downloaded from the public repository and all the processed data that I produced during my project. For example in the **./COVID_19/data/** directory, it includes: 
```R
./s41587-020-0602-4_COVID_19.rds
./all_raw_data
./de_data
./ppcseq_data
./summary_table
./final_sum_table
```
The *s41587-020-0602-4_COVID_19.rds* represents the raw data downloaded from the GEO. Except this file, all the rest directories include the processed data separated by the cell types. All the data are labelled clearly with the cell type name and the function name. For example, the file


## R script direcotries

```R
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/COVID_19/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/covid_atlas/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/GSE120575_melanoma/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/GSE129788_aging_mouse_brain/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/SCP1288_renal_cell_carcinoma/
/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/R_scripts/Overall_summary/
```
Same as data directory, there are 6 R script directories that are used in each dataset. Except for the **Overall_summary** directory, there are 10 files in each directory used for data cleaning, data processing, Seurat object generation, pseudo-bulk samples generation, differential gene-transcript abundance analysis, outliers identification and plots generation. Basically, these 10 files are identical used in each dataset, except for some minor modification. I added comments during my process which should be understandable. For example, in the folder, **COVID_19**, the R scripts include:
```R 
./COVID_19_generate_seurat.R     # Generate seurat object
./save_tidy_data.R.      # Generate pseudo-bulk samples for each cell type
./run_DE.R     # Perform differential gene-transcript abundance analysis by 3 methods, including DESeq2, edgeR robust and edgeR quasi-likelihood 
./run_ppcseq.R   # Perform outlier identification by using ppcseq
./create_summary.R  # Generate a summary dataframe about the outlier distributions for each cell type 
./create_summary_table_for_all_cell_types.R   # Summarize the outliers of all the cell types
./create_5_summary_plots.R    # Generate 5 plots for all the cell types
./COVID_19_generate_makeflow.R  # Generate the makeflow script for running all the processes in each cell type
./COVID_19.makeflow  # The output makeflow script ready for makeflow execution
./COVID_19.makeflow.makeflowlog  # The output of the makeflow 
```


