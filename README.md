# Evaluating Gene Set Enrichment Analysis Methods (methylGSA & ebGSEA) for Whole Genome Bisulfite Sequencing Data

## Maya Arvanitis, McGill University  
Supervisor: Dr. Claudia Kleinman

### Overview
This project explores the application of gene set enrichment analysis (GSEA) methods to whole genome bisulfite sequencing (WGBS) data. WGBS provides single-base resolution on DNA methylation patterns across the genome, an invaluable tool for studying epigenetic changes and their impact on gene regulation. However, GSEA methodologies are traditionally developed for RNA-seq data, creating a gap in available tools for epigenomic data. This project aims to assess and adapt GSEA methods like methylGSA and ebGSEA to work effectively with WGBS data, with the potential of establishing a new pipeline for analyzing DNA methylation patterns and their biological implications.

---

### Biological Context
DNA methylation is a critical epigenetic modification involved in regulating gene expression and genome stability. Aberrant methylation patterns are frequently associated with diseases, particularly cancer and neurodevelopmental disorders. 

**Whole Genome Bisulfite Sequencing (WGBS)** captures DNA methylation across the genome, enabling detailed analyses of CpG sites. This level of coverage offers an advantage over array-based methods by uncovering methylation in previously unexplored genomic regions, including intergenic areas and non-CpG methylation.

However, WGBS data pose unique challenges for GSEA, such as sparse CpG coverage per gene and the need to align CpG sites with genes or functional gene sets. Current GSEA tools, including methylGSA and ebGSEA, are not inherently designed for WGBS, highlighting the need to evaluate and potentially adapt these methods.

---

### Objectives
1. **Evaluate the Effectiveness of methylGSA and ebGSEA on WGBS Data**  
   - Test the performance of methylGSA and ebGSEA on WGBS datasets, specifically looking at how well they identify known methylation patterns related to gene sets.
  
2. **Develop a Standardized WGBS-GSEA Analysis Pipeline**  
   - If methylGSA or ebGSEA demonstrates utility with WGBS data, establish a reproducible analysis workflow. This would include preprocessing, gene annotation, and running gene set enrichment analysis with optimized parameters.

---

### Workflow & Methodology
To evaluate and implement GSEA methods for WGBS data, this project will follow the outlined steps below. The workflow includes data preparation, gene annotation, GSEA execution, and result evaluation. 

#### 1. **Data Preparation**
   - **Input Data**: Preprocessed WGBS datasets that include percentage methylation per CpG, coverage, and differential methylation.
   - **Filtering**: Filter for CpGs with significant p-values. CpGs are separated into those with positive and negative methylation changes to differentiate hypermethylated (positive) and hypomethylated (negative) regions.
   - **Sampling**: A subset of the data is used for initial tests to optimize parameters and reduce processing time.

#### 2. **Gene Annotation**
   - **Objective**: Link CpG sites to genes to interpret methylation changes at the gene level.
   - **Methodology**:
     - CpG sites are mapped to genes using biomaRy
     - **Data Format**: Output includes mapped genes, CpG methylation values, and associated annotations for use in GSEA.

#### 3. **Running GSEA Methods on WGBS Data**
   - **methylGSA**: Designed to support DNA methylation data, this method will be assessed for its applicability to WGBS data.
   - **ebGSEA**: Originally developed for EWAS (epigenome-wide association studies), ebGSEA will be evaluated for its sensitivity to WGBS methylation data.
   - **Parameters**: Both methods will be optimized for WGBS, with specific parameters set for CpG gene coverage, p-value thresholds, and statistical corrections.

#### 4. **Evaluation of GSEA Results**
   - **Criteria**:
     - Sensitivity and specificity are assessed by comparing observed methylation patterns against known patterns for histone H3.3 G34R mutations, associated with differential methylation at neuronal and immune gene promoters.
   - **Validation**: Enrichment results are validated through comparison with established biological findings on H3.3 G34R mutation effects.

#### 5. **Pipeline Documentation**
   - A standardized pipeline will be created if these methods prove effective, including detailed steps for preprocessing, annotation, and GSEA execution with methylGSA and ebGSEA.

---
### Getting Started

# MethylGSA
---

# ebGSEA
---
## doGT_custom.R

This function is an extension of the doGT.R function applied in ebGSEA [2]. This custom implementation allows for custom arrays. 

The inputs for the function are the following when running with custom arrays. 

```{r}
sgt.m <- doGT_custom(
    pheno.v = pheno.v_binary,
    data.m = as.matrix(beta_matrix_valid_df),
    model = c("linear"),
    array = "custom",
    ncores = 1,
    custom_map = mapped_genes_list
  )
```

Where:
- `pheno.v` = a binary vector whose columns are in line with those of the sample columns in data.m matrix. E.g in this case we had two Wild Type samples and Two G34R Samples. Our respective pheno.v vector was [ 0, 0, 1, 1] where 0 = WT and 1 = G34R
- ## Example Structure of `data.m`

The `data.m` matrix should have CpG IDs as row names and sample beta values as columns. Below is an example structure:

| **CpG ID**          | **Sample 1** | **Sample 2** | **Sample 3** | **Sample 4** |
|----------------------|--------------|--------------|--------------|--------------|
| `cg00000029`         | 0.85         | 0.90         | 0.87         | 0.92         |
| `cg00000108`         | 0.76         | 0.81         | 0.78         | 0.80         |
| `cg00000289`         | 0.65         | 0.67         | 0.66         | 0.68         |
| `cg00000305`         | 0.50         | 0.55         | 0.53         | 0.52         |
| `cg00000409`         | 0.72         | 0.74         | 0.73         | 0.75         |

Note: In this project CpG IDs were formatted in chr_start_end formatting for consistency. 

### Explanation:
   - **Rows (`CpG ID`)**: Unique identifiers for CpG sites (e.g., `cg00000029`).
   - **Columns (`Samples`)**: Beta values for each sample, representing the methylation level (ranging from 0 to 1).

---
- `model` = as implemented by Dong et al [2]
- `array` = here we have the selection now of `"custom"` allowing for custom gene mappings. As well, the original `"450k"` and `"850k"` as implemented in ebGSEA.
- `nscores` = number of cores for parallel with default set to 4 [2]
- `custom_map` a list of genes in ENTREZ ID format with their associated cpg ID
  
## Example Structure of `custom_map`

The `custom_map` should be a named list where:
- **Names of the list**: Gene identifiers in ENTREZ ID format.
- **Values of the list**: Vectors containing the associated CpG IDs.

### Example

```R
custom_map <- list(
  "11418" = c("cg00000029", "cg00000108", "cg00000289"),
  "11569" = c("cg00000305", "cg00000409"),
  "11600" = c("cg00000532", "cg00000621"),
  "11608" = c("cg00000745", "cg00000856")
)
```
-------------------------
## `doGSEAwt_custom` Function Documentation

## Overview

The `doGSEAwt_custom` function performs Gene Set Enrichment Analysis (GSEA) using the Wilcoxon Rank Sum Test (WT) and the Known-Population Median Test (KPMT). This implementation builds upon the original `doGSEAwt` function by introducing enhanced input validation, flexible annotation support, improved debugging outputs, and pre-processing optimizations.

---

## Key Features and Enhancements

### 1. Robust Input Validation
- **Gene Validation**:
  - Validates both `rankEID.m` and `ptw.ls` against valid ENTREZ IDs from the provided annotation database (`org.db`).
  - Logs invalid IDs and filters them from the analysis.
- **Pathway Validation**:
  - Retains only pathways with at least `minN` valid genes before analysis.
- **Error Handling**:
  - Throws specific errors if no valid genes or pathways remain after filtering.

**Advantage**: Ensures data integrity and avoids errors from invalid inputs.

---

### 2. Detailed Debugging Messages
- Provides comprehensive feedback at each step:
  - Number of valid genes and pathways after filtering.
  - Overlap of genes between the ranked matrix and pathways.
  - Raw and adjusted p-values.
  - Ranks of genes inside and outside pathways during testing.
- Logs intermediate results to help users understand the progress and outcomes of the analysis.

**Advantage**: Simplifies troubleshooting and increases transparency.

---

### 3. Flexible Pathway Filtering
- Filters pathways using the `minN` threshold at two stages:
  - **Before analysis**: Removes pathways with insufficient valid genes.
  - **During testing**: Ensures robust results by excluding pathways with too few overlaps.

**Advantage**: Improves computational efficiency and focuses on biologically meaningful pathways.

---

### 4. Enhanced Gene Symbol Mapping
- Uses `AnnotationDbi::mapIds` to map ENTREZ IDs to gene symbols.
- Excludes genes without valid SYMBOL mappings, ensuring clean results.
- Ensures all results are human-readable and biologically interpretable.

**Advantage**: Reduces ambiguity and enhances the utility of output data.

---

### 5. Cross-Species Annotation Support
- Accepts any `org.db` annotation object, making the function compatible with various organisms (e.g., humans, mice).
- Automatically retrieves valid ENTREZ IDs for the selected species.

**Advantage**: Broadens applicability across different species and datasets.

---

### 6. Improved Error Handling
- Detects and stops execution if inputs are invalid or unsuitable for analysis.
- Alerts users with clear error messages when no significant pathways are identified.

**Advantage**: Prevents incomplete or misleading results and improves user experience.

---

## Comparison with Original Function

### Differences from `doGSEAwt`
| Feature                          | `doGSEAwt_custom`                     | `doGSEAwt`                          |
|----------------------------------|---------------------------------------|-------------------------------------|
| Input Validation                 | Validates and filters both `rankEID.m` and `ptw.ls`. | Assumes inputs are valid.          |
| Debugging Outputs                | Detailed messages for every step.     | Limited high-level messages.       |
| Pathway Filtering                | Pre-filters pathways by `minN`.       | Filters pathways only during testing. |
| Gene Symbol Mapping              | Uses `AnnotationDbi::mapIds`.         | Uses a custom `convertIDs` function. |
| Annotation Database Support      | Flexible, supports multiple species.  | Fixed to human datasets.           |
| Overlap Validation               | Logs overlap counts between genes and pathways. | No overlap validation.             |
| Error Handling                   | Stops with specific error messages.   | May return incomplete results.     |

---

## Parameters

| Parameter  | Description                                                                                                      | Default Value |
|------------|------------------------------------------------------------------------------------------------------------------|---------------|
| `rankEID.m`| Matrix of ranked genes (rows) by statistics from a global test. Rownames must be valid ENTREZ IDs.               | -             |
| `ptw.ls`   | List of pathways, where each entry contains a vector of ENTREZ IDs for the pathway.                             | -             |
| `ncores`   | Number of cores for parallel execution.                                                                         | `4`           |
| `minN`     | Minimum number of genes required in a pathway to be included in the analysis.                                    | `5`           |
| `adjPVth`  | Adjusted p-value threshold for determining significant pathways.                                                 | `0.05`        |
| `org.db`   | Annotation database object (e.g., `org.Mm.eg.db` for mice).                                                     | `org.Mm.eg.db`|

---

## Outputs

The function returns a list containing the following:
1. **Rank(P)**:
   - A matrix of pathways ranked by adjusted Wilcoxon p-values (`adjP`).
   - Columns include:
     - `nREP`: Number of genes mapped to the pathway.
     - `AUC`: Area Under Curve of the Wilcoxon test.
     - `P(WT)`: p-value of the Wilcoxon test.
     - `P(KPMT)`: p-value of the Known-Population Median Test.
     - `adjP`: Adjusted p-value of the Wilcoxon test.
2. **Rank(AUC)**:
   - A matrix of pathways ranked by AUC.
3. **Genestat**:
   - A list of genes in each enriched pathway, with their statistics and p-values from the global test.

---

## Example Usage

```r
library(org.Mm.eg.db)

# Simulate ranked gene matrix
sim_rankEID.m <- matrix(
    runif(500 * 2, min = 0, max = 1),
    nrow = 500,
    dimnames = list(sample(keys(org.Mm.eg.db, keytype = "ENTREZID"), 500), c("Score1", "Score2"))
)

# Simulate pathways
sim_ptw.ls <- list(
    Pathway1 = sample(rownames(sim_rankEID.m), 50),
    Pathway2 = sample(rownames(sim_rankEID.m), 30)
)

# Run GSEA
results <- doGSEAwt_custom(
    rankEID.m = sim_rankEID.m,
    ptw.ls = sim_ptw.ls,
    ncores = 2,
    minN = 10,
    adjPVth = 0.05
)
```
---
### References
[1] Xu Ren and Pei Fen Kuan. methylgsa documentation: Gene set analysis for methylation data. Available at: https://www.bioconductor.org/packages/release/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html.

[2] Dong Dong, Yuan Tian, Shu-Cheng Zheng, and Andrew E. Teschendorff. ebgsea: An improved gene set enrichment analysis method for epigenome-wide-association studies. Bioinformatics, 35(18):3514–3516, Sep 2019.

[3] Sima Khazaei, Carol CL Chen, Augusto Faria Andrade, et al. Single substitution in h3.3g34 alters dnmt3a recruitment to cause progressive neurodegeneration. Cell, 186(6):1162–1178, 2023.

[4]. **AnnotationDbi R Package**  
   Available at: [Bioconductor AnnotationDbi Documentation](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)  

[5]. **Wilcoxon Rank Sum Test**  
   Mann, H.B., Whitney, D.R.  
   *On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other*.  
   Annals of Mathematical Statistics, Vol. 18, No. 1, pp. 50-60, 1947.

[6]. **Known-Population Median Test (KPMT)**  
   Implemented using the **kpmt R package**. Documentation available at:  
   [kpmt GitHub Repository](https://github.com/your-link-to-repo)

[7]. **Molecular Signatures Database (MSigDB)**  
   Subramanian, A., Tamayo, P., et al.  
   *Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles*.  
   **Proceedings of the National Academy of Sciences** (2005), 102(43):15545-15550.  
   DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

