# we explore significantly differentially expressed genes

# extracting dge results as dataframe
# we remove the gene names as rowid, it will be a column so that we can merge with additional data
# we also remove all characters after the . in the gene names

dge_model_results_sig <- data.frame(dge_model_results) %>% 
  arrange(padj) %>% 
  dplyr::filter(padj<0.05) %>% 
  rownames_to_column(var="ensid") 


### VISUALISATION 1: heatmap

# based on the differentially expressed genes
# heatmaps are based on normalized counts

norm_count_sig <- normalized_counts[dge_model_results_sig$ensid, ]

heat_colors <- brewer.pal(8, "YlOrRd")

pheatmap(norm_count_sig,
         color=heat_colors,
         cluster_rows=T,
         show_rownames = F,
         annotation=dplyr::select(coldata, types),
         scale="row")

### VISUALISATION 2: volcano plot

# data prep for volcano plot
dge_model_results_volcano <- data.frame(dge_model_results) %>% 
  rownames_to_column(var="ensid") %>% 
  dplyr::mutate(threshold = padj < 0.05)

# plot
ggplot(dge_model_results_volcano)+
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 color = threshold),
             size=0.8)+
  scale_color_manual(values=c("#CDB9B4", "#F43911"))+
  theme_classic()+
  xlab("log2 fold change")+
  ylab("-log10 adj. p-value")+
  theme(legend.position = "right",
        plot.title = element_text(size=rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  ylim(0,15)+
  ggtitle("Volcano plot")

### EXPLORING THE DIFFERENTIALLY EXPRESSED GENES based on type and function (only simple, no gsea)

# extracting dge results as dataframe
# we remove the gene names as rowid, it will be a column so that we can merge
# we also remove all characters after the . in the gene names

## adding gene symbol and gene type information to results data frame

df_results <- data.frame(dge_model_results) %>% 
  arrange(padj) %>% 
  dplyr::filter(padj<0.05) %>% 
  rownames_to_column(var="ensid") %>% 
  dplyr::mutate(ensid=gsub("\\..*","",ensid))


df_symbol <- data.frame(rowData(dds)[c('gene_name', 'gene_type')]) %>%
  rownames_to_column(var="ensid") %>% 
  dplyr::mutate(ensid=gsub("\\..*","",ensid))

df_results <- df_results %>% 
  left_join(df_symbol, by="ensid") 

# which types are the genes that are differentially expressed

# focus on the two major types: lncRNA and protein coding

df_results %>% 
  dplyr::select(ensid, gene_name, log2FoldChange, gene_type) %>% 
  dplyr::mutate(sign=case_when(log2FoldChange<0 ~ " neg.",
                               log2FoldChange>0 ~ " pos.")) %>% 
  ungroup() %>% 
  group_by(gene_type, sign) %>% 
  dplyr::summarise(N=n()) %>% 
  # focus only on the two largest groups: protein coding and lncRNA
  dplyr::filter(gene_type=="protein_coding" | gene_type=="lncRNA") %>% 
  ggplot(aes(fill=sign, y=N, x=gene_type)) + 
  scale_fill_manual(values=c("#CDB9B4", "#F43911"))+
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  ggtitle("The major types of sign. expressed genes")+
  xlab("gene type")


#################### visualisation 3: dot charts

### exploring certain genes: PGF, CALB1 

# data prep

df_symbol <- data.frame(rowData(dds)[c('gene_name', 'gene_type')]) %>%
  rownames_to_column(var="ensid")

df_dot_chart <- data.frame(normalized_counts) %>% 
  rownames_to_column(var="ensid") %>% 
  left_join(df_symbol, by="ensid") 

# PGF

df_coldata <- data.frame(coldata) %>% 
  rownames_to_column(var="sample")

df_pgf <- df_dot_chart %>% 
  dplyr::select(c(2:8)) %>% 
  dplyr::filter(gene_name=="PGF") %>% 
  dplyr::select(-gene_name) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="sample") 


colnames(df_pgf) <- c("sample", "norm_count")

df_pgf$sample <- gsub("\\.", "-", df_pgf$sample)

df_pgf <- df_pgf %>% 
    left_join(df_coldata, by="sample")

ggplot(df_pgf, aes(x=types, y=norm_count, fill=types)) + 
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values=c("#9FF3BD", "#F43911"))+
  theme_classic()+
  theme(legend.position = "none")+
  ggtitle("gene: PGF")+
  xlab("sample type")+
  ylab("normalized count")
  
# CALB1

df_calb1 <- df_dot_chart %>% 
  dplyr::select(c(2:8)) %>% 
  dplyr::filter(gene_name=="CALB1") %>% 
  dplyr::select(-gene_name) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="sample") 


colnames(df_calb1) <- c("sample", "norm_count")

df_calb1$sample <- gsub("\\.", "-", df_calb1$sample)

df_calb1 <- df_calb1 %>% 
  left_join(df_coldata, by="sample")

ggplot(df_calb1, aes(x=types, y=norm_count, fill=types)) + 
  geom_dotplot(binaxis='y', stackdir='center')+
  scale_fill_manual(values=c("#9FF3BD", "#F43911"))+
  theme_classic()+
  theme(legend.position = "none")+
  ggtitle("gene: CALB1")+
  xlab("sample type")+
  ylab("normalized count")
  

