# project1_DESeq2
Differential gene expression analysis with DESeq2 on TCGA-KIRC data

This repository performs differential gene expression analysis with [DESEq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) on normal (N=3) and tumor (N=3) samples from the [TCGA-KIRC](https://portal.gdc.cancer.gov/projects/TCGA-KIRC) (Kidney Renal Clear Cell Carcinoma) project. 

The repository contains the following scripts and files:

1. 01_data_download.R - obtaining the data and prepare the SummarizedExperimentObject
2. 02_inspecting_SE_object.R - this is an optional script to inspect the structure of our data as a SummarizedExperiment object
3. 03_DGE_DESeq2.R - performing the differential gene expression analysis
4. 04_exploring_results.R - exploring results with summary tables and figures
5. analysis.Rmd - R markdown to create the final html report
6. analysis.html - final report in html format

Remark:
The repository does not contain the data but those can be obtianed by running the script '01_data_download.R'.
The scripts and the Rmd file contain links to extra material and additional remarks.

Sources:

https://www.youtube.com/watch?v=UWXv9dUpxNE&t=2054s

https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf

https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html



 
