
# Load required packages
library(ebGSEA)
library(dplyr)
library(GenomicRanges)
library(biomaRt)
library(readr)
library(ggplot2)
library(msigdb)
library(tibble)
library(parallel)
library(globaltest)
library(org.Mm.eg.db)
library(kpmt)

# Connect to Ensembl for mouse genome annotation
load("data/ensembl_mart.RData")

# Retrieve hallmark gene sets for mouse with Ensembl IDs
msigdb.mm <- getMsigdb(org = 'mm', id = 'EZID', version = '7.4')
hallmarks = subsetCollection(msigdb.mm, 'h')
print("Hallmark gene sets:")
print(hallmarks)

# Convert hallmarks to a list of gene symbols by hallmark name
hallmarks_list <- lapply(hallmarks, function(gs) {
  # Extract gene symbols and set name for each GeneSet object
  setNames(list(gs@geneIds), gs@setName)
})
print("Hallmarks list:")
print(head(hallmarks_list))

# Flatten the list structure to obtain a single list with hallmark names as keys
ptw.ls <- unlist(hallmarks_list, recursive = FALSE)
print("Flattened hallmarks list:")
print(head(ptw.ls))

# Load the percentage methylation data
percentage_meth <- read_tsv("data/WTvsG34R_CRX_10W.percentage_meth.tsv")
print("Percentage methylation data:")
print(head(percentage_meth))

percentage_meth <- percentage_meth %>%
  mutate(cpg_id = paste0(chr, "_", start, "_", end))
print("Percentage methylation data with cpg_id:")
print(head(percentage_meth))

# Display the dimensions and the first few rows of the subset
cat("Dimensions of the subset:", dim(percentage_meth), "\n")
print(head(percentage_meth))

# Construct chromosomal_region in the query dataset
percentage_meth <- percentage_meth %>%
  mutate(query_id = paste0(chr, "_", start, "_", end),
         chromosomal_region = paste0(chr, ":", start, "-", end))
print("Percentage methylation data with chromosomal regions:")
print(head(percentage_meth))

# Retrieve gene annotations for hypermethylated CpGs
annotations <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "entrezgene_id"),
  filters = "chromosomal_region",
  values = percentage_meth$chromosomal_region,
  mart = ensembl
)
print("Annotations:")
print(head(annotations))

# Add chromosomal_region column to the result
annotations <- annotations %>%
  mutate(chromosomal_region = paste0(chromosome_name, "_", start_position, "_", end_position))
print("Annotations with chromosomal regions:")
print(head(annotations))

# Filter out rows with NA in the external_gene_name column
annotations <- annotations %>%
  filter(!is.na(entrezgene_id) & entrezgene_id != "" & !grepl("^\\s*$", entrezgene_id))
print("Filtered annotations:")
print(head(annotations))

# Assume hyper_cpgs_subset has a 'cpg_id' column with the format "chr_start_end"
percentage_meth <- as.data.frame(percentage_meth)

# Extract chr, start, and end from cpg_id
percentage_meth <- percentage_meth %>%
  mutate(
    chr = sub("^(.*)_.*_.*$", "\\1", cpg_id),
    start = as.numeric(sub("^.*_(.*)_.*$", "\\1", cpg_id)),
    end = as.numeric(sub("^.*_.*_(.*)$", "\\1", cpg_id))
  )
print("Percentage methylation data with extracted chr, start, and end:")
print(head(percentage_meth))

# Parse 'chromosomal_region' in hyper_annotations to chr, start, and end
annotations <- annotations %>%
  mutate(
    chr = sub("^(.*)_.*_.*$", "\\1", chromosomal_region),
    start = as.numeric(sub("^.*_(.*)_.*$", "\\1", chromosomal_region)),
    end = as.numeric(sub("^.*_.*_(.*)$", "\\1", chromosomal_region))
  )
print("Parsed annotations:")
print(head(annotations))

# Convert CpG data to GRanges object
cpg_gr <- GRanges(
  seqnames = percentage_meth$chr,
  ranges = IRanges(start = percentage_meth$start, end = percentage_meth$end),
  cpg_id = percentage_meth$cpg_id
)
print("CpG GRanges object:")
print(cpg_gr)

# Convert gene annotation data to GRanges object
gene_gr <- GRanges(
  seqnames = annotations$chr,
  ranges = IRanges(start = annotations$start, end = annotations$end),
  entrezgene_id = annotations$entrezgene_id
)
print("Gene GRanges object:")
print(gene_gr)

# Find overlaps between CpG ranges and gene ranges
overlaps <- findOverlaps(cpg_gr, gene_gr)
print("Overlaps:")
print(overlaps)

# Map each CpG to all overlapping genes
mapped_genes <- data.frame(
  cpg_id = mcols(cpg_gr)$cpg_id[queryHits(overlaps)],
  entrezgene_id = mcols(gene_gr)$entrezgene_id[subjectHits(overlaps)]
)
print("Mapped genes:")
print(head(mapped_genes))
print(paste("Dimensions of mapped_genes:", dim(mapped_genes)))

percentage_meth <- as_tibble(percentage_meth)

mapped_genes <- mapped_genes %>%
  left_join(percentage_meth %>% dplyr::select(cpg_id), by = "cpg_id")
print("Mapped genes with p-values:")
print(head(mapped_genes))
print(paste("Dimensions of mapped_genes with p-values:", dim(mapped_genes)))

# Convert `mapped_genes` into a list where each CpG ID is associated with one or more gene names
custom_mapping <- mapped_genes %>%
  group_by(cpg_id) %>%
  summarise(genes = list(unique(entrezgene_id))) %>%
  deframe()
print("Custom mapping:")
print(head(custom_mapping, 15))

invalid_mappings <- which(sapply(custom_mapping, length) == 0)
if (length(invalid_mappings) > 0) {
  cat("Found empty mappings for CpGs:\n")
  print(names(custom_mapping)[invalid_mappings])
  custom_mapping <- custom_mapping[-invalid_mappings]
  cat("Updated custom_mapping size:", length(custom_mapping), "\n")
}

genes_in_mapping <- unique(unlist(custom_mapping))
genes_in_ptw <- unique(unlist(ptw.ls))

common_genes <- intersect(genes_in_mapping, genes_in_ptw)
cat("Number of common genes between custom_mapping and ptw.ls:", length(common_genes), "\n")

if (length(common_genes) == 0) {
  stop("Error: No overlap between custom_mapping and ptw.ls genes.")
}

# Ensure the data is a regular data frame
percentage_meth_df <- as.data.frame(percentage_meth)

percentage_meth_df <- percentage_meth_df[, !colnames(percentage_meth_df) %in% c("start", "end")]

# Identify the columns with numeric data (assuming they contain methylation percentages)
numeric_columns <- sapply(percentage_meth_df, is.numeric)

# Extract only the numeric columns (methylation percentages)
beta_values <- percentage_meth_df[, numeric_columns] / 100  # Convert percentages to beta values

# Convert the extracted values to a matrix
beta_matrix <- as.matrix(beta_values)

# Set row names as CpG site identifiers in the "chr_start_end" format
rownames(beta_matrix) <- percentage_meth_df$cpg_id

# Display the dimensions and a preview of the matrix
cat("Dimensions of the beta matrix:", dim(beta_matrix), "\n")
print(head(beta_matrix))

# Create a phenotype vector
pheno.v <- c("SK059_2427_WT_CRX", "SK060_2616_WT_CRX", "SK061_2617_G34R_CRX", "SK062_2428_G34R_CRX")

pheno.v_binary <- ifelse(grepl("WT", pheno.v), 0, 1)

# Display the phenotype vector
print(pheno.v_binary)

# Double-check formatting of CpG IDs in custom_mapping
names(custom_mapping) <- gsub(":", "_", names(custom_mapping))
names(custom_mapping) <- gsub("-", "_", names(custom_mapping))

# Ensure that the format matches with rownames of beta_matrix
common_cpgs <- intersect(rownames(beta_matrix), names(custom_mapping))
cat("Number of common CpGs:", length(common_cpgs), "\n")

# If overlap is zero, display a sample of CpG IDs
if (length(common_cpgs) == 0) {
  cat("Sample CpG IDs from beta_matrix:\n")
  print(head(rownames(beta_matrix)))
  
  cat("Sample CpG IDs from custom_mapping:\n")
  print(head(names(custom_mapping)))
}

# Check for missing values in beta_matrix
if (any(is.na(beta_matrix)) || any(is.nan(beta_matrix)) || any(is.infinite(beta_matrix))) {
  stop("Error: Beta matrix contains missing, NaN, or infinite values.")
}

# Check for empty gene lists in custom_mapping
empty_genes <- which(sapply(custom_mapping, length) == 0)
if (length(empty_genes) > 0) {
  cat("Found empty mappings for CpGs:\n")
  print(names(custom_mapping)[empty_genes])
  
  # Remove empty mappings
  custom_mapping <- custom_mapping[-empty_genes]
  cat("Updated custom_mapping size:", length(custom_mapping), "\n")
}

# Subset custom_mapping and beta_matrix to include only common CpGs
custom_mapping_valid <- custom_mapping[common_cpgs]
beta_matrix_valid <- beta_matrix[common_cpgs, , drop = FALSE]

# Check the dimensions and preview
cat("Dimensions of valid beta matrix:", dim(beta_matrix_valid), "\n")
print(head(beta_matrix_valid))

# Trim any leading/trailing whitespaces from both sets
rownames(beta_matrix_valid) <- trimws(rownames(beta_matrix_valid))
names(custom_mapping_valid) <- trimws(names(custom_mapping_valid))

# Recheck the intersection after trimming
common_cpgs <- intersect(rownames(beta_matrix_valid), names(custom_mapping_valid))
cat("Number of common CpGs:", length(common_cpgs), "\n")

# Print the first few CpG IDs in both datasets for inspection
print(head(rownames(beta_matrix_valid)))
print(head(names(custom_mapping_valid)))

# Check for NA, NaN, or Inf in beta_matrix
if (any(is.na(beta_matrix_valid)) || any(is.nan(beta_matrix_valid)) || any(is.infinite(beta_matrix_valid))) {
  cat("Error: Beta matrix contains NA, NaN, or Inf values.\n")
} else {
  cat("No NA/NaN/Inf values in beta matrix.\n")
}

beta_matrix_valid <- beta_matrix_valid[complete.cases(beta_matrix_valid) & !apply(beta_matrix_valid, 1, function(row) any(is.nan(row) | is.infinite(row))), ]

# Assuming mapped_genes has columns 'cpg_id' and 'external_gene_name'

# Step 1: Create a mapping between row names of beta_matrix_valid and cpg_id in mapped_genes
beta_matrix_valid_df <- as.data.frame(beta_matrix_valid)  # Convert matrix to data frame
beta_matrix_valid_df$cpg_id <- rownames(beta_matrix_valid_df)  # Add cpg_id as a new column

# Step 2: Merge the data frame with mapped_genes on cpg_id
merged_df <- merge(beta_matrix_valid_df, mapped_genes, by = "cpg_id", all.x = TRUE)
print("Merged data frame:")
print(head(merged_df))

# Step 3: Check for duplicates in external_gene_name
duplicate_genes <- merged_df[duplicated(merged_df$entrezgene_id), ]
if (nrow(duplicate_genes) > 0) {
  message("Warning: There are duplicate gene names in the mapping.")
  # You can decide how to handle duplicates here:
  # Option 1: Remove duplicates (if multiple rows correspond to the same gene, you might keep only one)
  merged_df <- merged_df[!duplicated(merged_df$entrezgene_id), ]
  # Option 2: You might want to keep all duplicates, but handle them separately.
}

# Step 4: Set row names of the merged data frame to external_gene_name
rownames(merged_df) <- merged_df$entrezgene_id

# Convert it back to a matrix if needed
beta_matrix_valid_with_genes <- as.matrix(merged_df)

beta_matrix_valid_with_genes <- beta_matrix_valid_with_genes[, !colnames(beta_matrix_valid_with_genes) %in% c("cpg_id", "entrezgene_id")]

# Display the final matrix
print("Beta matrix with genes:")
print(head(beta_matrix_valid_with_genes))

custom_mapping <- setNames(mapped_genes$cpg_id, mapped_genes$entrezgene_id)
print("Custom mapping:")
print(head(custom_mapping))

# Remove the cpg_id col from beta_matrix_valid_df
beta_matrix_valid_df <- beta_matrix_valid_df[, !colnames(beta_matrix_valid_df) %in% "cpg_id"]
print("Beta matrix valid data frame:")
print(head(beta_matrix_valid_df))

# Convert the cleaned data frame to a list grouped by entrezgene_id
mapped_genes_list <- split(mapped_genes$cpg_id, mapped_genes$entrezgene_id)

# Remove any entries with no associated gene (clean the list after creation)
mapped_genes_list <- mapped_genes_list[lengths(mapped_genes_list) > 0]

# Print the result
print("Mapped genes list:")
print(mapped_genes_list)

print("Phenotype vector:")
print(pheno.v_binary)
print("Beta matrix valid data frame:")
print(head(as.matrix(beta_matrix_valid_df)))
print("Mapped genes list:")
print(head(mapped_genes_list))

source("custom/doGT_custom.R")
source("custom/doGSEAwt_custom.R")

tryCatch({
  sgt.m <- doGT_custom(
    pheno.v = pheno.v_binary,
    data.m = as.matrix(beta_matrix_valid_df),
    model = c("linear"),
    array = "custom",
    ncores = 1,
    custom_map = mapped_genes_list
  )
  
  cat("Global Test completed successfully\n")
  print(head(sgt.m))
}, error = function(e) {
  cat("Error in doGT:", e$message, "\n")
})

# Save Global Test results to a file
write.csv(sgt.m, "sgt_results.csv", row.names = TRUE)

# Retrieve hallmark gene sets for mouse with Ensembl IDs
msigdb.mm.ez <- getMsigdb(org = 'mm', id = 'EZID', version = '7.4')

# Flatten the list structure to obtain a single list 
ptw.ls_ez <- unlist(msigdb.mm.ez, recursive = FALSE)

# Check the first few row names
print("Row names of sgt.m:")
print(head(rownames(sgt.m)))

# Find the intersection of genes between sgt.m and ptw.ls_ez
common_genes <- intersect(rownames(sgt.m), ptw.ls_ez)
print("Common genes between sgt.m and ptw.ls_ez:")
print(common_genes)

message("Structure of all_genes_ptw:")
str(all_genes_ptw)

# Iterate over the elements and check their class
for (i in seq_along(all_genes_ptw)) {
  message("Element ", i, " class: ", class(all_genes_ptw[[i]]))
}

common_genes <- as.character(unlist(common_genes))
all_genes_ptw <- lapply(all_genes_ptw, function(genes) {
  if (inherits(genes, "S4")) {
    as.character(genes@.Data)  # Extract the actual data from the S4 object
  } else if (is.vector(genes)) {
    as.character(genes)  # Ensure it's a character vector
  } else {
    stop("Unsupported data type in all_genes_ptw.")
  }
})

# Step 2: Filter sgt.m to include only common genes
sgt.m_filtered <- sgt.m[unlist(common_genes), , drop = FALSE]

# Step 3: Filter ptw.ls_ez to include only pathways with common genes
ptw.ls_filtered <- lapply(all_genes_ptw, function(genes) {
  intersect(genes, common_genes)
})

# Remove pathways that are now empty
ptw.ls_filtered <- ptw.ls_filtered[sapply(ptw.ls_filtered, length) > 0]

# Check dimensions and pathway counts after filtering
message("Filtered sgt.m dimensions: ", dim(sgt.m_filtered))
message("Number of pathways after filtering: ", length(ptw.ls_filtered))

topGSEA.lm <- doGSEAwt_custom(
  rankEID.m = as.matrix(sgt.m),
  ptw.ls = ptw.ls_ez,
  ncores = 1,
  minN = 3,
  adjPVth = 0.05,
  org.db = org.Mm.eg.db
)

# Save GSEA results to a file
write.csv(topGSEA.lm, "topGSEA_results.csv", row.names = TRUE)
