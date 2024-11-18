
doGT_custom <- function(pheno.v, data.m, model = c("linear"), array = c("450k", "850k", "custom"), ncores = 4, custom_map = NULL) {
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
    if (is.null(custom_map)) {
      stop("Error: Custom mapping must be provided for custom array type.")
    }
    message("Using custom mapping to genes...")
    subsets <- mclapply(custom_map, intersect, rownames(data.m), mc.cores = ncores)
    message("Done")
  } else {
    stop("Error: Unsupported array type. Use '450k', '850k', or 'custom'.")
  }
  
  # Ensure subsets contain valid entries
  nrep.v <- unlist(lapply(subsets, length))
  selG.idx <- which(nrep.v > 0)
  
  # Check if subsets are empty after filtering
  if (length(selG.idx) == 0) {
    stop("Error: No valid genes found in the subsets after filtering.")
  }
  
  # Run the Global Test
  message("Running Global Test...")
  gt.o <- gt(response = pheno.v, alternative = t(data.m), model = model, directional = FALSE, 
             standardize = FALSE, permutations = 0, subsets = subsets[selG.idx], trace = FALSE)
  message("Done")
  
  # Extract and sort results
  resGT.m <- as.matrix(result(gt.o))
  tmp.s <- sort(resGT.m[, 2], decreasing = TRUE, index.return = TRUE)
  sresGT.m <- resGT.m[tmp.s$ix, ]
  
  return(sresGT.m)
}
