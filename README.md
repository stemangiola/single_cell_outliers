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
The *s41587-020-0602-4_COVID_19.rds* represents the raw data downloaded from the GEO. Except this file, all the rest directories include the processed data separated by the cell types. All the data are labelled clearly with the cell type name and the function name. 

## Rscript direcotries





