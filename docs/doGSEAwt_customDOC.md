# `doGSEAwt_custom` Function Documentation

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
## References

1. **Dong, D., Tian, Y., Zheng, S.C., Teschendorff, A.E.**  
   *ebGSEA: an improved Gene Set Enrichment Analysis method for Epigenome-Wide-Association Studies*.  
   **BMC Bioinformatics** (2019), 35(18):3514-3516.  
   DOI: [10.1093/bioinformatics/btz073](https://doi.org/10.1093/bioinformatics/btz073)

2. **AnnotationDbi R Package**  
   Available at: [Bioconductor AnnotationDbi Documentation](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)  

3. **Wilcoxon Rank Sum Test**  
   Mann, H.B., Whitney, D.R.  
   *On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other*.  
   Annals of Mathematical Statistics, Vol. 18, No. 1, pp. 50-60, 1947.

4. **Known-Population Median Test (KPMT)**  
   Implemented using the **kpmt R package**. Documentation available at:  
   [kpmt GitHub Repository](https://github.com/your-link-to-repo)

5. **Molecular Signatures Database (MSigDB)**  
   Subramanian, A., Tamayo, P., et al.  
   *Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles*.  
   **Proceedings of the National Academy of Sciences** (2005), 102(43):15545-15550.  
   DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
