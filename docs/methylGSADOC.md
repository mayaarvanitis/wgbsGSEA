# Gene Set Analysis with methylGSA

## Overview

### Key Steps
1. **Load and Preprocess Differential Methylation Data**
2. **Filter Significant CpGs**
3. **Annotate CpGs to Genes**
4. **Prepare Annotations**
5. **Retrieve Gene Sets**
6. **Run GSA Methods**

---

## Step-by-Step Process

### 1. Load and Preprocess Differential Methylation Data
- The pipeline begins by reading the differential methylation data file. CpGs are filtered based on their statistical significance (`pvalue < 0.05`), separating hyper- and hypo-methylated CpGs.

### 2. Filter Significant CpGs
- Significant CpGs are divided into hypermethylated and hypomethylated subsets. For debugging and testing purposes, a subset of hypermethylated CpGs is sampled.

```r
hyper_cpgs_subset <- hyper_cpgs %>% sample_n(1000)
```

---

### 3. Annotate CpGs to Genes
- CpGs are mapped to genes using the `biomaRt` package to query the Ensembl database. Overlaps between CpG coordinates and gene regions are identified using `GRanges`.

```r
overlaps <- findOverlaps(cpg_gr, gene_gr)
mapped_genes <- data.frame(
  cpg_id = mcols(cpg_gr)$cpg_id[queryHits(overlaps)],
  external_gene_name = mcols(gene_gr)$external_gene_name[subjectHits(overlaps)]
)
```

---

### 4. Prepare Annotations
- Unique CpG-to-gene mappings are prepared, ensuring no duplicates. The `prepareAnnot` function creates a compatible annotation format for methylGSA.

```r
CpG2Gene <- mapped_genes %>%
  dplyr::select(cpg_id, external_gene_name) %>%
  distinct(cpg_id, .keep_all = TRUE)

FullAnnot3 <- prepareAnnot(CpG2Gene, geneidtype = "SYMBOL")
```

---

### 5. Retrieve Gene Sets
- Gene sets for mouse (`mmusculus`) are retrieved from MSigDB. The hallmark gene sets are extracted and converted into a format compatible with methylGSA.

```r
msigdb.mm <- getMsigdb(org = 'mm', id = 'SYM', version = '7.4')
hallmarks <- subsetCollection(msigdb.mm, 'h')
```

---

### 6. Run GSA Methods
Three methods are implemented for analysis:

#### a. **methylglm**
- A logistic regression model is applied to identify gene set enrichment.

```r
methylglm_results <- methylglm(
    cpg.pval = cpg_pval,
    array.type = "custom",
    FullAnnot = FullAnnot3,
    group = "all",
    GS.list = hallmarks_list,
    GS.idtype = "SYM",
    minsize = 2,
    maxsize = 800
)
```

#### b. **Over-Representation Analysis (ORA)**
- The ORA approach identifies over-represented pathways.

```r
res_ORA <- methylRRA(
    cpg.pval = cpg_pval,
    method = "ORA",
    FullAnnot = FullAnnot3,
    GS.list = hallmarks_list,
    GS.idtype = "SYMBOL",
    minsize = 2,
    maxsize = 800
)
```

#### c. **Gene Set Enrichment Analysis (GSEA)**
- GSEA is performed to determine if gene sets are significantly enriched with CpGs.

```r
res_GSEA <- methylRRA(
    cpg.pval = cpg_pval,
    method = "GSEA",
    FullAnnot = FullAnnot3,
    GS.list = hallmarks_list,
    GS.idtype = "SYMBOL",
    minsize = 2,
    maxsize = 800
)
```

---

## Output

1. **`mapped_genes.rds`**: An RDS file containing the mapping of CpGs to genes.
2. **methylglm Results**: Results from logistic regression for gene set enrichment.
3. **ORA Results**: Pathways identified as over-represented.
4. **GSEA Results**: Enriched pathways based on ranked p-values.

---

## Notes
- Ensure all CpG IDs and gene symbols are properly formatted to avoid errors in annotation and analysis.
- Parallel processing can be enabled for large datasets by setting `parallel = TRUE` in methylGSA methods.
