################################################################################
# the aim of this script is to perform:

## 1. data preprocessing
## 2. and differential gene expression analysis with DESeq2

################################################################################

# useful material/tutorials

# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#starting-from-summarizedexperiment
# https://nbisweden.github.io/workshop-RNAseq/2111/lab_preprocessing.html
# https://support.bioconductor.org/p/130210/
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# https://dputhier.github.io/ASG/practicals/rnaseq_diff_Snf2/rnaseq_diff_Snf2.html

# DESEq2 vignette

# vignette("DESeq2")

# https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf

# what do we know about how our data were processed

# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=rna_expression_workflow
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/


# on normalization

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3247-x

################################################################################
#                    STEP 1: creating DESeqDataSet object
################################################################################
# the DESeq2 package requires data in the so-called DESeqDataSet object form
# we create this object using the elements of the summarized experiment object

# we need three sorts of info: raw counts, samples and genes

########################### counts as data frame

counts <- as.data.frame(assay(tcga_kirc_data, withDimnames = T, "unstranded"))

# cleaning up sample names
## we keep only the first 16 characters of the column names

colnames(counts)[1:6] <- sapply(colnames(counts)[1:6], function(x){
  substr(x, 1,16)
})

# inspect the result

as_tibble(counts)

##################### information on samples as data frame

# extracting sample names
samples <- tcga_kirc_data$sample

# types of samples (tumor vs. normal)
types <- tcga_kirc_data$tissue_type

# making data frame
coldata <- as.data.frame(types, row.names = samples)

# setting factor levels (used as design later)
coldata$types <- factor(coldata$types, levels = c("Normal", "Tumor"))

# inspect
coldata

# IMPORTANT!!!: the order of samples should be the same in the count and coldata data frames

# checking if the order of rowdata and coldata are the same

rownames(coldata)

# [1] "TCGA-CW-5584-01A" "TCGA-CZ-5984-01A" "TCGA-B0-4841-01A" "TCGA-CZ-5985-11A" "TCGA-CZ-5451-11A"
# [6] "TCGA-B0-5691-11A"

colnames(counts)

# [1] "TCGA-CW-5584-01A" "TCGA-CZ-5984-01A" "TCGA-B0-4841-01A" "TCGA-CZ-5985-11A" "TCGA-CZ-5451-11A"
# [6] "TCGA-B0-5691-11A"

# do we have the samples in the same order in the two datasets? yes

# alternative

all(rownames(coldata) == colnames(counts))
# [1] TRUE

# if they do not align, use match()

################################ creating object

# putting it all together
# the design will be tumor vs. normal

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ types)

# inspect
dds

# it is not necessary but we can add extra info on genes, e.g. gene names, gene type

genename <- data.frame(gene_name=tcga_kirc_data@rowRanges$gene_name)
genetype <- data.frame(gene_type=tcga_kirc_data@rowRanges$gene_type)
mcols(dds) <- DataFrame(mcols(dds), genename, genetype)
mcols(dds)


# inspect again
dds

rowData(dds)
################################################################################
#                          STEP 2: preprocessing and visualisation
################################################################################

# # goal is to find out now if there is a problem with the samples and the experiment
# 
# # exploring the distribution of counts
# 
# ggplot(counts)+
#   geom_histogram(aes(x=`TCGA-CW-5584-01A`), stat = "bin", bins = 200)+
#   xlab("Raw expression counts")+
#   ylab("Number of genes")+
#   ggtitle("sample: TCGA-CW-5584-01A")
# 
# ggplot(counts)+
#   geom_histogram(aes(x=`TCGA-CZ-5984-01A`), stat = "bin", bins = 200)+
#   xlab("Raw expression counts")+
#   ylab("Number of genes")+
#   ggtitle("sample: TCGA-CZ-5984-01A")

######################## filtering

# we filter out genes with low counts

keep <- rowSums( counts(dds) >= 3 ) >= 3
dds <- dds[keep,]

# counts(dds)

# we are left with 25987 genes in place of 60660

# we can filter out genes with low variance as well, we do not do it now

######################## normalization for library size

## step 1: estimate size factors (median ratio method)
# extra info is added to colData 
dds <- estimateSizeFactors(dds)
dds

# inspect size factors

## method 1
dds$sizeFactor

## method 2
sizeFactors(dds)

## step 2: calculating normalized counts (this does not affect the dds object itself)

normalized_counts <- counts(dds, normalized = T)

############################### cluster analysis

## step 1: transformation (log transformation, variance stabilizing transformation - vst)

vst_dds <- vst(dds, blind=T)
vst_dds
assay(vst_dds)

## step 2: extract the vst matrix and compute pairwise correlation between samples

### extraction
vst_dds_mat <- assay(vst_dds)

### correlation

vst_dds_cor <- cor(vst_dds_mat)

## step 3: heatmap

pheatmap(vst_dds_cor, annotation = dplyr::select(coldata, types),
                     main = "Sample correlation and clustering")


# conclusion: TCGA-B0-4841-01A is a bit weird, it clusters relatively weakly with tumor samples

############################# pca

plotPCA(vst_dds, intgroup="types")+
  geom_text_repel(aes(label = rownames(coldata)))+
  theme_classic()+
  ggtitle("PCA plot")


# conclusion: good separation on pc1

# pca plots with more annotation possibilities
# https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_01_pca.html


################################################################################
#                            STEP 3: DGE analysis
################################################################################

# STEP 1: modeling
# estimate size factors and estimate variation across samples for each gene
# NB estimation

# we already have de DESeq2 object

# fitting the NB model

dds <- DESeq(dds)

# STEP 2: exploring dispersion estimates 

# how well do our data fit the model?
# log2 fold changes

plotDispEsts(dds, main = "Dispersion plot")


# we expect dispersion values to decrease with increasing mean
# when you have only a few replicates, the estimates are often inaccurate - use info across all genes 
# to estimate the most likely dispersion for a given mean expression value (red line)
# the original estimates are shrunken towards the red line to avoid false positives
# the shrinkage depends on distance to the curve and sample size

# STEP 3: contrasts testing log2 fold changes based on modelled counts

dge_model_results <- results(dds, alpha=0.05)
dge_model_results
# the design has already been defined (set as factor)

# MA plot

plotMA(dge_model_results, main = "MA plot - before shrinkage")


# each dot is a gene 
# x-axis: mean of normalized counts across samples
# y-axis: log fold change tumor vs. normal and significance is labeled with blue

# STEP 4: shrink log2 fold changes

dge_model_results <- lfcShrink(dds,
                               coef=2,
                                        res = dge_model_results)

dge_model_results

plotMA(dge_model_results, main = "MA plot - before shrinkage")

# tip: explore the columns of the result tables

# if I understand correctly with the above code we update the results
# compare the table (log fold change, SE and p-values)
# the number of differentially expressed genes does not change 

############################ STEP 5: explore results

mcols(dge_model_results)

# we use adjusted p-vales: BH method

summary(dge_model_results)

