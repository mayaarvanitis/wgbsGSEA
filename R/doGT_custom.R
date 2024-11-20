doGT_custom <- function(pheno.v, data.m, model = c("linear"), array = c("450k", "850k", "custom"), ncores = 4, custom_map = NULL) {
  # Debug: Initial inputs
  message("Debug: Checking initial inputs...")
  message("Number of samples in pheno.v: ", length(pheno.v))
  message("Dimensions of data.m: ", paste(dim(data.m), collapse = " x "))
  message("Array type: ", array)
  
  # Check that phenotype vector is numeric and matches data matrix dimensions
  if (!is.numeric(pheno.v)) {
    stop("Error: Phenotype vector (pheno.v) must be numeric.")
  }
  if (length(pheno.v) != ncol(data.m)) {
    stop("Error: The length of phenotype vector does not match the number of rows in the data matrix.")
  }
  
  # Ensure valid values in phenotype and data matrix
  if (any(is.na(pheno.v)) || any(is.nan(pheno.v)) || any(is.infinite(pheno.v))) {
    stop("Error: Phenotype vector contains NA, NaN, or Inf values.")
  }
  if (any(is.na(data.m)) || any(is.nan(data.m)) || any(is.infinite(data.m))) {
    stop("Error: Data matrix contains NA, NaN, or Inf values.")
  }
  
  if (array == "custom" && is.null(custom_map)) {
    stop("Error: Custom mapping must be provided for custom array type.")
  }
  
  # Handle different array types or custom mapping
  if (array == "450k") {
    message("Mapping 450k probes to genes...")
    data("dualmap450kEID")
    subsets <- mclapply(mapEIDto450k.lv, intersect, rownames(data.m), mc.cores = ncores)
    message("Done")
  } else if (array == "850k") {
    message("Mapping EPIC probes to genes...")
    data("dualmap850kEID")
    subsets <- mclapply(mapEIDto850k.lv, intersect, rownames(data.m), mc.cores = ncores)
    message("Done")
  } else if (array == "custom") {
    message("Using custom mapping to genes...")
    message("Number of entries in custom_map: ", length(custom_map))
    message("First few CpG IDs in custom_map: ", paste(head(names(custom_map)), collapse = ", "))
    message("First few row names in data.m: ", paste(head(rownames(data.m)), collapse = ", "))
    
    # Check if there are empty CpG IDs in custom_map
    custom_map <- custom_map[!sapply(custom_map, function(x) length(x) == 0)]
    
    subsets <- mclapply(names(custom_map), intersect, rownames(data.m), mc.cores = ncores)
    message("Mapping complete. Checking subsets...")
    
    # Debug: Check subsets
    empty_subsets <- sum(sapply(subsets, length) == 0)
    message("Number of empty subsets: ", empty_subsets)
    if (empty_subsets > 0) {
      message("First few empty subset names: ", paste(head(names(subsets[sapply(subsets, length) == 0])), collapse = ", "))
    }
    message("Done")
  } else {
    stop("Error: Unsupported array type. Use '450k', '850k', or 'custom'.")
  }
  
  # Ensure subsets contain valid entries
  nrep.v <- unlist(lapply(subsets, length))
  message("Number of genes mapped to subsets: ", sum(nrep.v > 0))
  selG.idx <- which(nrep.v > 0)
  
  # Debug: Validate subsets
  if (length(selG.idx) == 0) {
    message("No valid genes found in the subsets after filtering. Subsets preview:")
    print(head(subsets))
    stop("Error: No valid genes found in the subsets after filtering.")
  }
  
  # Debug: Subsets information
  message("Number of valid subsets: ", length(selG.idx))
  message("First few subset names: ", paste(head(names(subsets[selG.idx])), collapse = ", "))
  message("Sizes of first few valid subsets: ", paste(head(nrep.v[selG.idx]), collapse = ", "))
  
  # Run the Global Test
  message("Running Global Test...")
  tryCatch({
    gt.o <- gt(response = pheno.v, alternative = t(data.m), model = model, directional = FALSE, 
               standardize = FALSE, permutations = 0, subsets = subsets[selG.idx], trace = FALSE)
  }, error = function(e) {
    message("Error during Global Test: ", e$message)
    stop(e)
  })
  message("Done")
  
  # Debug: Check Global Test object
  if (is.null(gt.o)) {
    stop("Error: Global Test object is NULL.")
  }
  
  # Extract and sort results
  message("Extracting and sorting results...")
  resGT.m <- as.matrix(result(gt.o))
  if (is.null(resGT.m) || nrow(resGT.m) == 0) {
    stop("Error: No results returned from Global Test.")
  }
  
  tmp.s <- sort(resGT.m[, 2], decreasing = TRUE, index.return = TRUE)
  sresGT.m <- resGT.m[tmp.s$ix, ]
  
  # Debug: Final results
  message("Global Test completed successfully. First few results:")
  print(head(sresGT.m))
  
  return(sresGT.m)
}
