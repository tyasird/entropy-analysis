# Load necessary libraries
# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))
library(GEOquery)
library(limma)




# Define dataset ID
dataset_id <- "GSE32676"
gpl = "GPL570"
# Define group membership with 'X' as exclusion
gsms <- "11111111111111111111111110000000"

# Load the series and platform data from GEO
gset <- getGEO(dataset_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep(gpl, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make proper column names to match topTable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Log2 transformation if necessary
ex <- exprs(gset)
# Check if log transformation is needed
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}




# Split the 'gsms' string into individual group assignments
sml <- strsplit(gsms, split="")[[1]]

# Exclude samples labeled as 'X' and keep only '0' (control) and '1' (treatment)
valid_samples <- sml != "X"  # Logical vector to identify non-'X' samples

# Filter the expression data and metadata to keep only valid samples
gset_filtered <- gset[, valid_samples]  # Apply the filter to gset

# Reassign group labels (excluding 'X' samples)
sml_filtered <- sml[valid_samples]

# Assign filtered samples to groups (control and treatment)
gs <- factor(sml_filtered)
groups <- make.names(c("C", "T"))  # "C" for control, "T" for treatment
levels(gs) <- groups
gset_filtered$group <- gs  # Assign group information to the filtered gset object

# Design matrix (for later statistical analysis)
design <- model.matrix(~group + 0, gset_filtered)
colnames(design) <- levels(gs)


# Assuming gset has been filtered and group information has been assigned
# Access the group labels from the phenoData slot
group_labels <- pData(gset_filtered)$group  # 'group' is the column in phenoData where groups are stored
# Ensure group_labels is a character vector (in case it's a factor)
group_labels <- as.character(group_labels)

# Count the number of 'C' and 'T' samples
num_C <- sum(group_labels == "C")  # Number of control samples
num_T <- sum(group_labels == "T")  # Number of treatment samples


# Now sort the samples so that "C" samples come before "T" samples
gset_sorted <- gset[, order(group_labels)]  # Sort samples based on the group labels

# To check the sorted expression set:
# View the expression matrix for the sorted gset
sorted_expression_matrix <- exprs(gset_sorted)


# Extract the gene symbols and create an annotated expression set (sorted)
gene_symbols <- gset_sorted@featureData@data$Gene.symbol
annotated_expression_set <- data.frame(Gene.symbol = gene_symbols, exprs(gset_sorted))

# Define output filename with dataset ID
output_file <- paste0(dataset_id, "_C", num_C, "_T", num_T, ".csv")

# Save the annotated, sorted expression set as a CSV
write.table(annotated_expression_set, file = output_file, sep = ",", quote = FALSE, row.names = TRUE, col.names = NA)
message("Annotation and log2 transformation complete. Sorted annotated data saved to '", output_file, "'")
