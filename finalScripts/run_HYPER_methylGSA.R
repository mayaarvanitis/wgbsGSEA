# Load required libraries
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationDbi)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(msigdb)
library(dplyr)
library(tibble)
library(methylGSA)
library(org.Mm.eg.db)
library(biomaRt)

# Read the data
diff_data <- readr::read_tsv("/Users/mayaarvanitis/Desktop/COMP401/data/WTvsG34R_CRX_10W.diff_meth.tsv", col_names = TRUE)

# Inspect the data
print("Initial data:")
print(head(diff_data))
print(paste("Dimensions of diff_data:", dim(diff_data)))

significant_cpgs <- diff_data %>% filter(pvalue < 0.05)
print("Significant CpGs:")
print(head(significant_cpgs))
print(paste("Dimensions of significant_cpgs:", dim(significant_cpgs)))

hyper_cpgs <- significant_cpgs %>% filter(meth.diff > 0)
print("Hypermethylated CpGs:")
print(head(hyper_cpgs))
print(paste("Dimensions of hyper_cpgs:", dim(hyper_cpgs)))

hypo_cpgs <- significant_cpgs %>% filter(meth.diff < 0)
print("Hypomethylated CpGs:")
print(head(hypo_cpgs))
print(paste("Dimensions of hypo_cpgs:", dim(hypo_cpgs)))


# Remove columns we don't need
hyper_cpgs <- hyper_cpgs[, c("chr", "start", "end", "pvalue")]
hyper_cpgs <- hyper_cpgs %>%
  mutate(cpg_id = paste0(chr, "_", start, "_", end))
print("Processed hypermethylated CpGs:")
print(head(hyper_cpgs))
print(paste("Dimensions of processed hyper_cpgs:", dim(hyper_cpgs)))


load("data/ensembl_mart.RData")

# Construct chromosomal_region in the query dataset
hyper_cpgs <- hyper_cpgs %>%
  mutate(query_id = paste0(chr, "_", start, "_", end),
         chromosomal_region = paste0(chr, ":", start, "-", end))
print("Hyper CpGs with chromosomal regions:")
print(head(hyper_cpgs))
print(paste("Dimensions of hyper_cpgs with chromosomal regions:", dim(hyper_cpgs)))


# Retrieve gene annotations for hypermethylated CpGs
hyper_annotations <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = "chromosomal_region",
  values = hyper_cpgs$chromosomal_region,
  mart = ensembl
)
print("Hyper annotations:")
print(head(hyper_annotations))
print(paste("Dimensions of hyper_annotations:", dim(hyper_annotations)))


# Add chromosomal_region column to the result
hyper_annotations <- hyper_annotations %>%
  mutate(chromosomal_region = paste0(chromosome_name, "_", start_position, "_", end_position))

# Filter out rows with NA in the external_gene_name column
hyper_annotations <- hyper_annotations %>%
  filter(!is.na(external_gene_name) & external_gene_name != "" & !grepl("^\\s*$", external_gene_name))
print("Filtered hyper annotations:")
print(head(hyper_annotations))
print(paste("Dimensions of filtered hyper_annotations:", dim(hyper_annotations)))

# Display the first few rows to confirm
head(hyper_annotations)

# Assume hyper_cpgs has a 'cpg_id' column with the format "chr_start_end"
hyper_cpgs <- as.data.frame(hyper_cpgs)
# Extract chr, start, and end from cpg_id
hyper_cpgs <- hyper_cpgs %>%
  mutate(
    chr = sub("^(.*)_.*_.*$", "\\1", cpg_id),
    start = as.numeric(sub("^.*_(.*)_.*$", "\\1", cpg_id)),
    end = as.numeric(sub("^.*_.*_(.*)$", "\\1", cpg_id))
  )
print("Hyper CpGs with extracted chr, start, and end:")
print(head(hyper_cpgs))
print(paste("Dimensions of hyper_cpgs with extracted chr, start, and end:", dim(hyper_cpgs)))


# Parse 'chromosomal_region' in hyper_annotations to chr, start, and end
hyper_annotations <- hyper_annotations %>%
  mutate(
    chr = sub("^(.*)_.*_.*$", "\\1", chromosomal_region),
    start = as.numeric(sub("^.*_(.*)_.*$", "\\1", chromosomal_region)),
    end = as.numeric(sub("^.*_.*_(.*)$", "\\1", chromosomal_region))
  )
print("Parsed hyper annotations:")
print(head(hyper_annotations))
print(paste("Dimensions of parsed hyper_annotations:", dim(hyper_annotations)))


# Convert CpG data to GRanges object
cpg_gr <- GRanges(
  seqnames = hyper_cpgs$chr,
  ranges = IRanges(start = hyper_cpgs$start, end = hyper_cpgs$end),
  cpg_id = hyper_cpgs$cpg_id
)
print("CpG GRanges object:")
print(cpg_gr)

# Convert gene annotation data to GRanges object
gene_gr <- GRanges(
  seqnames = hyper_annotations$chr,
  ranges = IRanges(start = hyper_annotations$start, end = hyper_annotations$end),
  external_gene_name = hyper_annotations$external_gene_name
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
  external_gene_name = mcols(gene_gr)$external_gene_name[subjectHits(overlaps)]
)
print("Mapped genes:")
print(head(mapped_genes))
print(paste("Dimensions of mapped_genes:", dim(mapped_genes)))


hyper_cpgs <- as_tibble(hyper_cpgs)

mapped_genes <- mapped_genes %>%
  left_join(hyper_cpgs %>% dplyr::select(cpg_id, pvalue), by = "cpg_id")
print("Mapped genes with p-values:")
print(head(mapped_genes))
print(paste("Dimensions of mapped_genes with p-values:", dim(mapped_genes)))


cpg_pval <- setNames(mapped_genes$pvalue, mapped_genes$cpg_id)
print("CpG p-values:")
print(head(cpg_pval))

# Remove duplicate rows based on cpg_id, keeping the first occurrence
# Use dplyr::select to avoid namespace conflicts
CpG2Gene <- mapped_genes %>%
  dplyr::select(cpg_id, external_gene_name) %>%
  distinct(cpg_id, .keep_all = TRUE)
print("CpG to Gene mapping:")
print(head(CpG2Gene))
print(paste("Dimensions of CpG2Gene:", dim(CpG2Gene)))


# Prepare the annotation with unique CpG to gene mappings
FullAnnot3 <- prepareAnnot(CpG2Gene, geneidtype = "SYMBOL")
print("Full annotation:")
print(head(FullAnnot3))
print(paste("Dimensions of FullAnnot3:", dim(FullAnnot3)))

# Retrieve hallmark gene sets for mouse with Ensembl IDs
msigdb.mm <- getMsigdb(org = 'mm', id = 'SYM', version = '7.4')
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
hallmarks_list <- unlist(hallmarks_list, recursive = FALSE)
print("Flattened hallmarks list:")
print(head(hallmarks_list))


# Run methylglm with custom settings
methylglm_results <- methylglm(
    cpg.pval = cpg_pval,              # Named vector of CpG p-values
    array.type = "custom",             # Custom array type
    FullAnnot = FullAnnot3,             # Data frame with CpG and gene annotations
    group = "all",                     # Consider all CpGs
    GS.list = hallmarks_list,      # List of hallmark gene sets
    GS.idtype = "SYM",             # Ensembl gene IDs in hallmark gene sets
    GS.type = "GO",                    # Type of gene sets; specify if needed
    minsize = 2,                       # Minimum gene set size
    maxsize = 800,                     # Maximum gene set size
    parallel = FALSE                   # Parallel processing (adjust as needed)
)
print("Methylglm results:")
print(head(methylglm_results, 50))
print(paste("Dimensions of methylglm_results:", dim(methylglm_results)))

# Save methylglm results to a file
write.csv(methylglm_results, "methylglm_results.csv", row.names = FALSE)

# ORA approach
res_ORA = methylRRA(
  cpg.pval = cpg_pval,     
  method = "ORA",           
  FullAnnot = FullAnnot3,   
  GS.list = hallmarks_list,           
  GS.idtype = "SYMBOL",       
  minsize = 2,              
  maxsize = 800            
)
print("ORA results:")
print(head(res_ORA, 15))
print(paste("Dimensions of res_ORA:", dim(res_ORA)))

# Save ORA results to a file
write.csv(res_ORA, "res_ORA.csv", row.names = FALSE)

# GSEA approach
res_GSEA <- methylRRA(
  cpg.pval = cpg_pval,         
  method = "GSEA",            
  FullAnnot = FullAnnot3,       
  GS.list = hallmarks_list,          
  GS.idtype = "SYMBOL",        
  minsize = 2,              
  maxsize = 800               
)
print("GSEA results:")
print(head(res_GSEA, 15))
print(paste("Dimensions of res_GSEA:", dim(res_GSEA)))

# Save GSEA results to a file
write.csv(res_GSEA, "res_GSEA.csv", row.names = FALSE)
