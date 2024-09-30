# Evaluating Gene Set Enrichment Analysis Methods (methylGSA & ebGSEA) for Whole Genome Bisulfite Sequencing Data

## Maya Arvanitis, McGill University  
Supervisor: Dr. Claudia Kleinman

### Biological Context
Gene set enrichment analysis (GSEA) is an essential tool in genomics for identifying the enrichment of biological pathways from high-throughput data, such as gene expression and DNA methylation profiles. Whole genome bisulfite sequencing (WGBS) is a comprehensive method for assessing DNA methylation across the entire genome, providing valuable insights into epigenetic regulation. 

However, while GSEA methodologies have been well-established for RNA-seq data, there is a gap in the availability of effective tools specifically designed for WGBS data. Given the critical role of DNA methylation in gene regulation, the development of a robust GSEA framework for WGBS data would significantly advance our understanding of complex biological processes and their implications for disease.

### Objectives
1. **Evaluate the effectiveness of methylGSA and ebGSEA on WGBS data**: Test methylGSA [1] and ebGSEA [2] on WGBS datasets to assess their suitability for gene set enrichment analysis.
2. **Develop a standardized analysis pipeline if the methods prove effective**: If these methods are found to be effective, create a standardized pipeline for their application to WGBS data.

### Methodology
Pre-existing preprocessed WGBS datasets will be utilized to evaluate the compatibility and performance of methylGSA and ebGSEA. If these methods are not directly compatible, specific challenges will be identified, and potential solutions will be explored.

### Evaluation Methods
To evaluate the effectiveness of GSEA methods for whole genome bisulfite sequencing (WGBS) data, we will leverage the distinct methylation patterns associated with histone H3.3 mutations (G34R/W), which are implicated in cancer and neurodevelopmental disorders. The G34R mutation leads to hypermethylation at the promoters of neuronal genes and hypomethylation of immune genes, while the G34W mutation does not.

The evaluation will focus on each method’s ability to accurately detect these methylation patterns in WGBS data. Sensitivity and specificity will be assessed by comparing the results to established findings, ensuring both computational efficiency and biological relevance, as detailed in the following reference paper [3].

### References
1. Xu Ren and Pei Fen Kuan. methylgsa documentation: Gene set analysis for methylation data. Available at: https://www.bioconductor.org/packages/release/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html.
2. Dong Dong, Yuan Tian, Shu-Cheng Zheng, and Andrew E. Teschendorff. ebgsea: An improved gene set enrichment analysis method for epigenome-wide-association studies. *Bioinformatics*, 35(18):3514–3516, Sep 2019.
3. Sima Khazaei, Carol CL Chen, Augusto Faria Andrade, et al. Single substitution in h3.3g34 alters dnmt3a recruitment to cause progressive neurodegeneration. *Cell*, 186(6):1162–1178, 2023.
