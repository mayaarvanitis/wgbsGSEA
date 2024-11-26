# Load necessary libraries
library(ebGSEA)
library(dplyr)
library(GenomicRanges)
library(biomaRt)
library(readr)
library(msigdb)
library(org.Mm.eg.db)
library(parallel)
library(globaltest)
library(ggplot2)

source("/project/kleinman/maya.arvanitis/from_hydra/COMP401/doGT_custom.R")
source("/project/kleinman/maya.arvanitis/from_hydra/COMP401/doGSEAwt_custom.R")

# Step 1: Connect to Ensembl for annotation
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://useast.ensembl.org")

# Step 2: Load data
diff_data <- read_tsv("data_maya/WTvsG34R_CRX_10W.percentage_meth.tsv.gz")
head(diff_data)

# Step 3: Prepare the phenotype vector
pheno.v <- c("SK059_2427_WT_CRX", "SK060_2616_WT_CRX", "SK061_2617_G34R_CRX", "SK062_2428_G34R_CRX")
pheno.v_binary <- ifelse(grepl("WT", pheno.v), 0, 1)

# Step 4: Prepare CpG-to-gene mapping
diff_data <- diff_data %>%
  mutate(chromosomal_region = paste0(chr, ":", start, "-", end)) %>%
  mutate(cpg_id = paste0(chr, "_", start, "_", end))

annotations <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "entrezgene_id"),
  filters = "chromosomal_region",
  values = diff_data$chromosomal_region,
  mart = ensembl
)

annotations <- annotations %>%
  mutate(cpg_id = paste0(chromosome_name, "_", start_position, "_", end_position)) %>%
  filter(!is.na(entrezgene_id))

cpg_gr <- GRanges(
  seqnames = diff_data$chr,
  ranges = IRanges(start = diff_data$start, end = diff_data$end),
  cpg_id = diff_data$cpg_id
)

gene_gr <- GRanges(
  seqnames = annotations$chromosome_name,
  ranges = IRanges(start = annotations$start_position, end = annotations$end_position),
  entrezgene_id = annotations$entrezgene_id
)

overlaps <- findOverlaps(cpg_gr, gene_gr)

mapped_genes <- data.frame(
  cpg_id = mcols(cpg_gr)$cpg_id[queryHits(overlaps)],
  entrezgene_id = mcols(gene_gr)$entrezgene_id[subjectHits(overlaps)]
)

mapped_genes <- mapped_genes %>%
  distinct(cpg_id, .keep_all = TRUE)

# Step 5: Prepare beta value matrix
beta_matrix <- as.matrix(diff_data[, grepl("SK", colnames(diff_data))]) / 100
rownames(beta_matrix) <- diff_data$cpg_id

# Filter beta matrix to include only CpGs in the mapping
common_cpgs <- intersect(rownames(beta_matrix), mapped_genes$cpg_id)
beta_matrix <- beta_matrix[common_cpgs, , drop = FALSE]

# Step 6: Map CpGs to genes
mapped_genes <- mapped_genes %>%
  filter(cpg_id %in% common_cpgs) %>%
  group_by(entrezgene_id) %>%
  summarise(cpgs = list(unique(cpg_id)))

custom_mapping <- mapped_genes %>% deframe()

# Step 7: Run ebGSEA
sgt.m <- doGT_custom(
  pheno.v = pheno.v_binary,
  data.m = beta_matrix,
  model = c("linear"),
  array = "custom",
  ncores = 1,
  custom_map = custom_mapping
)

# Step 8: Load pathways and perform GSEA
msigdb.mm <- getMsigdb(org = 'mm', id = 'EZID', version = '7.4')
ptw.ls <- unlist(msigdb.mm, recursive = FALSE)

topGSEA <- doGSEAwt_custom(
  rankEID.m = as.matrix(sgt.m),
  ptw.ls = ptw.ls,
  ncores = 1,
  minN = 3,
  adjPVth = 0.05,
  org.db = org.Mm.eg.db
)

# Step 9: Save results
saveRDS(sgt.m, "sgt_results.rds")
saveRDS(topGSEA, "gsea_results.rds")

# Step 10: Plot results
result_df <- as.data.frame(topGSEA[[1]])
ggplot(result_df, aes(x = AUC, y = -log10(adjP))) +
  geom_point() +
  labs(x = "Area Under Curve", y = "-log10(Adjusted p-value)", title = "Pathway Enrichment Analysis") +
  theme_minimal()
