## doGT_custom.R

This function is an extension of the doGT.R function applied in ebGSEA [1]. This custom implementation allows for custom arrays. 

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
- pheno.v = a binary vector whose columns are in line with those of the sample columns in data.m matrix. E.g in this case we had two Wild Type samples and Two G34R Samples. Our respective pheno.v vector was [ 0, 0, 1, 1] where 0 = WT and 1 = G34R
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

- model = as implemented by Dong et al [1]
- array = here we have the selection now of "custom" allowing for custom gene mappings. As well, the original "450k" and "850k" as implemented in ebGSEA.
- nscores = number of cores for parallel with default set to 4 [1]
- custom_map a list of genes in ENTREZ ID format with their associated cpg ID
  
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
