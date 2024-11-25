# Load required libraries
library(GenomicRanges)
library(AnnotationDbi)
library(methylGSA)
library(msigdb)
library(org.Mm.eg.db)
library(dplyr)
library(readr)
library(biomaRt)

# Step 1: Load the differential methylation data
diff_data <- readr::read_tsv("data_maya/WTvsG34R_CRX_10W.diff_meth.tsv.gz", col_names = TRUE)

# Step 2: Filter significant CpGs
significant_cpgs <- diff_data %>% filter(pvalue < 0.05)
hyper_cpgs <- significant_cpgs %>% filter(meth.diff > 0)
hypo_cpgs <- significant_cpgs %>% filter(meth.diff < 0)

# Step 3: Map CpGs to genes
# Helper function for gene annotation
map_cpgs_to_genes <- function(cpg_data, ensembl_mart) {
  cpg_data <- cpg_data %>%
    mutate(chromosomal_region = paste0(chr, ":", start, "-", end)) %>%
    mutate(cpg_id = paste0(chr, "_", start, "_", end))  # Create cpg_id for mapping
  
  
  # Retrieve annotations
  annotations <- getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
    filters = "chromosomal_region",
    values = cpg_data$chromosomal_region,
    mart = ensembl_mart
  )
  
  # Format annotations for methylGSA
  annotations <- annotations %>%
    filter(!is.na(external_gene_name) & external_gene_name != "" & !grepl("^\\s*$", external_gene_name)) %>%
    mutate(chromosomal_region = paste0(chromosome_name, "_", start_position, "_", end_position))
  
  # Create GRanges objects for overlaps
  cpg_gr <- GRanges(
    seqnames = cpg_data$chr,
    ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
    cpg_id = cpg_data$cpg_id
  )
  
  gene_gr <- GRanges(
    seqnames = annotations$chromosome_name,
    ranges = IRanges(start = annotations$start_position, end = annotations$end_position),
    external_gene_name = annotations$external_gene_name
  )
  
  overlaps <- findOverlaps(cpg_gr, gene_gr)
  mapped_genes <- data.frame(
    cpg_id = mcols(cpg_gr)$cpg_id[queryHits(overlaps)],
    external_gene_name = mcols(gene_gr)$external_gene_name[subjectHits(overlaps)]
  )
  
  # Merge pvalues from the original cpg_data
  mapped_genes <- mapped_genes %>%
    left_join(cpg_data %>% dplyr::select(cpg_id, pvalue), by = "cpg_id") %>%  # Retain pvalues
    distinct(cpg_id, .keep_all = TRUE)  # Remove duplicates
  
  return(mapped_genes)
}

# Step 4: Prepare annotation data
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")

# Annotate hypermethylated CpGs
hyper_annotations <- map_cpgs_to_genes(hyper_cpgs, ensembl)
hypo_annotations <- map_cpgs_to_genes(hypo_cpgs, ensembl)

# Step 5: Prepare methylGSA inputs
prepare_methylGSA_inputs <- function(mapped_genes) {
  
  # Generate CpG to gene mapping
  CpG2Gene <- mapped_genes %>%
    dplyr::select(cpg_id, external_gene_name) %>%
    distinct(cpg_id, .keep_all = TRUE)
  
  # Prepare annotation for methylGSA
  FullAnnot <- prepareAnnot(CpG2Gene, geneidtype = "SYMBOL")
  
  # Prepare CpG p-value vector
  cpg_pval <- setNames(mapped_genes$pvalue, mapped_genes$cpg_id)
  
  return(list(FullAnnot = FullAnnot, cpg_pval = cpg_pval))
}

hyper_inputs <- prepare_methylGSA_inputs(hyper_annotations)
hypo_inputs <- prepare_methylGSA_inputs(hypo_annotations)

# Step 6: Retrieve gene sets
msigdb.mm <- getMsigdb(org = 'mm', id = 'SYM', version = '7.4')
hallmarks <- subsetCollection(msigdb.mm, 'h')

# Convert hallmarks to a list for methylGSA
hallmarks_list <- lapply(hallmarks, function(gs) {
  setNames(list(gs@geneIds), gs@setName)
})
hallmarks_list <- unlist(hallmarks_list, recursive = FALSE)

# Step 7: Run methylGSA (methylglm, ORA, GSEA) for hyper and hypo CpGs
run_methylGSA <- function(inputs, hallmarks_list) {
  results <- list()
  
  # methylglm
  results$methylglm <- methylglm(
    cpg.pval = inputs$cpg_pval,
    array.type = "custom",
    FullAnnot = inputs$FullAnnot,
    group = "all",
    GS.list = hallmarks_list,
    GS.idtype = "SYM",
    minsize = 2,
    maxsize = 800,
    parallel = FALSE
  )
  
  # ORA
  results$ORA <- methylRRA(
    cpg.pval = inputs$cpg_pval,
    method = "ORA",
    FullAnnot = inputs$FullAnnot,
    GS.list = hallmarks_list,
    GS.idtype = "SYM",
    minsize = 2,
    maxsize = 800
  )
  
  # GSEA
  results$GSEA <- methylRRA(
    cpg.pval = inputs$cpg_pval,
    method = "GSEA",
    FullAnnot = inputs$FullAnnot,
    GS.list = hallmarks_list,
    GS.idtype = "SYM",
    minsize = 2,
    maxsize = 800
  )
  
  return(results)
}

# Run methylGSA for hypermethylated CpGs
hyper_results <- run_methylGSA(hyper_inputs, hallmarks_list)

# Run methylGSA for hypomethylated CpGs
hypo_results <- run_methylGSA(hypo_inputs, hallmarks_list)

# Step 8: Save results
saveRDS(hyper_results, "hyper_results.rds")
saveRDS(hypo_results, "hypo_results.rds")

# Display top results
print("Hyper Results - methylglm")
print(head(hyper_results$methylglm, 15))
print("Hypo Results - methylglm")
print(head(hypo_results$methylglm, 15))

