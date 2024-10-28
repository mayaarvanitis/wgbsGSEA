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
---
### References
[1] Xu Ren and Pei Fen Kuan. methylgsa documentation: Gene set analysis for methylation data. Available at: https://www.bioconductor.org/packages/release/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html.
[2] Dong Dong, Yuan Tian, Shu-Cheng Zheng, and Andrew E. Teschendorff. ebgsea: An improved gene set enrichment analysis method for epigenome-wide-association studies. Bioinformatics, 35(18):3514–3516, Sep 2019.
[3] Sima Khazaei, Carol CL Chen, Augusto Faria Andrade, et al. Single substitution in h3.3g34 alters dnmt3a recruitment to cause progressive neurodegeneration. Cell, 186(6):1162–1178, 2023.

