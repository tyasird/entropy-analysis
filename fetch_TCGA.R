BiocManager::install("TCGAbiolinks")
BiocManager::install("PCAtools")

library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(biomaRt)
library(DESeq2)
library(PCAtools)
library(data.table)
library(dplyr)

project <- "TCGA-PAAD"  # Pancreatic adenocarcinoma

query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query)
data <- GDCprepare(query)

# counts and metadata
counts <- as.data.frame(assay(data))
clinical <- data.frame(data@colData)
data <- data[, rownames(clinical)]
all(rownames(clinical) == colnames(data))


# remove genes below 10 count
keep <- rowSums(assay(data)) >= 10
data <- data[keep,]
clinical$definition <-  gsub(" ", "_", clinical$definition)
clinical$definition <- as.factor(clinical$definition)
levels(clinical$definition)
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

# load as a DESeq object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(assay(data)),
                              colData = clinical,
                              design = ~ definition)


# data tranfromation
vsd <- vst(dds, blind=FALSE)
# p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)
# biplot(p, colby = "definition", lab = NULL, legendPosition = 'right')


# deg analysis, normalization
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)


# annotation
eset = as.data.frame(assay(vsd))
eset$ensembl_gene_id <- sub("\\..*", "", rownames(eset))
eset <- as.data.frame(eset)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = eset$ensembl_gene_id,
  mart = ensembl
)
gene_mapping <- gene_mapping[!duplicated(gene_mapping$ensembl_gene_id), ]



# merge annotation
eset <- data.frame(ensembl_gene_id = eset$ensembl_gene_id, eset, row.names = NULL)
merged_data <- merge(gene_mapping, eset, by = "ensembl_gene_id", all.y = TRUE)
merged_data <- merged_data[!is.na(merged_data$hgnc_symbol) & merged_data$hgnc_symbol != "", ]
merged_data$ensembl_gene_id = NULL

# save 
write.table(merged_data, paste0(project, ".csv"), sep=",", row.names=FALSE)
clinical_df = clinical_subset <- clinical %>% select(barcode, sample_type)
write.table(clinical_df, paste0(project, "_clinical.csv"), sep=",", row.names=FALSE)

