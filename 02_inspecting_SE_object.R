################################################################################

# This part is optional
## the aim of this script is to get familiar with the structure and content
## of SummarizedExperiment objects

################################################################################

############################ ASSAY

## what types of assays we have in the summarized experiments?

### we have raw counts (e.g. unstranded) and normalized counts (e.g. fpkm-unstrand)

assayNames(tcga_kirc_data)

# [1] "unstranded"       "stranded_first"   "stranded_second"  "tpm_unstrand"     "fpkm_unstrand"   
# [6] "fpkm_uq_unstrand"

## retrieveing unstranded

assay(tcga_kirc_data, withDimnames = T, "unstranded")

# ## retrieving fpkm_unstrand
# 
# assay(tcga_kirc_data, withDimnames = T, "fpkm_unstrand")

############################ GENES (rows)

# what do we know about the genes? 
genedata <- rowData(tcga_kirc_data)

str(genedata)

# how many genes do we have?

length(genedata@rownames)
genedata@nrows

# gene names

head(genedata@listData$gene_name)

# gene types

head(genedata@listData$gene_type)

############################ SAMPLES (columns)

# what do we know about the samples

sampledata <- colData(tcga_kirc_data)

# definition

sampledata$definition

# tumor descriptor

sampledata$tumor_descriptor

# patient

sampledata$patient

# sample type id

sampledata$sample_type_id
