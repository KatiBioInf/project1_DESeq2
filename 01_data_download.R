################################################################################
# Author: Kati Buzasi

# Description: Differential gene expression analysis - DESeq2


# Sources and main inspirations:

# https://www.youtube.com/watch?v=UWXv9dUpxNE&t=2054s  (YT channel: Bioinformagician)
# https://portal.gdc.cancer.gov/

################################################################################
#                                 Obtaining data
################################################################################    

# necessary libraries

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(apeglm)
# library(annotables)
# library(biomaRt)
# library(AnnotationDbi)
# library(EnsDb.Hsapiens.v75)
library(dplyr)
library(RColorBrewer)
library(ggrepel)



################################################################################

# get projects present on the GDC portal

gdcprojects <- getGDCprojects()
View(gdcprojects)

# we are interested in TCGA-KIRC (Kidney Renal Clear Cell Carcinoma)

getProjectSummary('TCGA-KIRC')

# we are interested in RNA-Seq data

## we can inspect the samples in this project

query_tcga_kirc <- GDCquery(project = 'TCGA-KIRC',
                                     data.category = 'Transcriptome Profiling',
                                     access = 'open',
                                     experimental.strategy = 'RNA-Seq',
                                     workflow.type = 'STAR - Counts')

output_query <- getResults(query_tcga_kirc)

#


# we select 6 samples (just by eyeballing, no strategy), 
# the first three are primary tumor and the second three are normal sample 

# here we can find more on the TCGA barcodes: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

query_tcga_kirc_selected <- GDCquery(project = 'TCGA-KIRC',
                            data.category = 'Transcriptome Profiling',
                            access = 'open',
                            experimental.strategy = 'RNA-Seq',
                            workflow.type = 'STAR - Counts',
                            barcode = c('TCGA-CW-5584-01A-01R-1541-07',
                                        'TCGA-CZ-5984-01A-11R-1672-07',
                                        'TCGA-B0-4841-01A-01R-1277-07',
                                        'TCGA-CZ-5985-11A-01R-1672-07',
                                        'TCGA-CZ-5451-11A-01R-1503-07',
                                        'TCGA-B0-5691-11A-01R-1541-07'))

output_query <- getResults(query_tcga_kirc_selected)

# download data
# this will save the data in the project folder
# within subfolder: GDCdata
# this takes a little while

GDCdownload(query_tcga_kirc_selected)

# making summarized experiment object
## more on summarized experiments here: 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html


tcga_kirc_data <- GDCprepare(query_tcga_kirc_selected, summarizedExperiment = T)

# SUCCESS: data are downloaded and ready for preprocessing and analysis

save(tcga_kirc_data, file = "tcga_kirc_data.RData")

# clear environment
rm(list=ls())