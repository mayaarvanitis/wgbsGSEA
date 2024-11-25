doGSEAwt_custom <- function(rankEID.m, ptw.ls, ncores = 4, minN = 5, adjPVth = 0.05, org.db = org.Mm.eg.db) {
  # Validate Inputs
  message("Validating inputs...")
  if (!is.matrix(rankEID.m)) stop("Error: 'rankEID.m' must be a matrix.")
  if (!is.list(ptw.ls)) stop("Error: 'ptw.ls' must be a list.")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Error: 'AnnotationDbi' package is required.")
  
  # Validate and filter rankEID.m rownames (ENTREZID)
  valid_entrez_ids <- AnnotationDbi::keys(org.db, keytype = "ENTREZID")
  invalid_rankEIDs <- setdiff(rownames(rankEID.m), valid_entrez_ids)
  if (length(invalid_rankEIDs) > 0) {
    warning("Invalid ENTREZIDs in rankEID.m: ", paste(invalid_rankEIDs, collapse = ", "))
  }
  rankEID.m <- rankEID.m[rownames(rankEID.m) %in% valid_entrez_ids, , drop = FALSE]
  message("Filtered rankEID.m to valid ENTREZIDs. Remaining rows: ", nrow(rankEID.m))
  
  # Validate and filter ptw.ls (ENTREZID)
  message("Validating and filtering pathways...")
  all_genes_in_ptw <- unique(unlist(ptw.ls))
  invalid_ptw_ids <- setdiff(all_genes_in_ptw, valid_entrez_ids)
  if (length(invalid_ptw_ids) > 0) {
    warning("Invalid ENTREZIDs in ptw.ls: ", paste(invalid_ptw_ids, collapse = ", "))
  }
  ptw.ls <- lapply(ptw.ls, function(genes) intersect(genes, valid_entrez_ids))
  ptw.ls <- ptw.ls[lengths(ptw.ls) > 0]  # Remove pathways with no valid genes
  message("Filtered pathways. Remaining pathways: ", length(ptw.ls))
  
  # Ensure there are still valid pathways and genes
  if (nrow(rankEID.m) == 0 || length(ptw.ls) == 0) {
    stop("Error: No valid ENTREZIDs left after filtering. Check your inputs.")
  }
  
  # Extract row names from rankEID.m
  rankEID.v <- rownames(rankEID.m)
  message("Number of valid ENTREZIDs: ", length(rankEID.v))
  
  # Run GSEA
  message("Running Wilcox Test and Known Population Median Test...")
  gseaWT.m <- matrix(
    unlist(parallel::mclapply(ptw.ls, gseaWTfn, rankEID.v, mc.cores = ncores, minN = minN)),
    ncol = 4, byrow = TRUE
  )
  message("GSEA completed.")
  colnames(gseaWT.m) <- c("nREP", "AUC", "P(WT)", "P(KPMT)")
  rownames(gseaWT.m) <- names(ptw.ls)
  
  # Adjust p-values
  message("Adjusting p-values using BH method...")
  tmp.s <- sort(gseaWT.m[, "P(WT)"], decreasing = FALSE, index.return = TRUE)
  sgseaWT.m <- gseaWT.m[tmp.s$ix, ]
  padj.v <- p.adjust(sgseaWT.m[, "P(WT)"], method = "BH")
  sel.idx <- which(padj.v <= adjPVth)
  message("Number of significant pathways: ", length(sel.idx))
  
  topGSEAwt.lm <- list()
  
  # Convert IDs and handle missing SYMBOLs
  message("Mapping ENTREZIDs to SYMBOLs...")
  sym.v <- AnnotationDbi::mapIds(org.db, keys = rankEID.v, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  del.idx <- which(is.na(sym.v))
  if (length(del.idx) > 0) {
    warning("Some ENTREZIDs could not be mapped to SYMBOLs and will be excluded.")
    sym.v <- sym.v[-del.idx]
    rankEID.m <- rankEID.m[-del.idx, ]
    rankEID.v <- rankEID.v[-del.idx]
    message("Remaining ENTREZIDs after SYMBOL mapping: ", length(rankEID.v))
  }
  
  # Process results
  if (length(sel.idx) > 1) {
    message("Processing results...")
    topGSEAwt.m <- cbind(sgseaWT.m[sel.idx, ], padj.v[sel.idx])
    colnames(topGSEAwt.m) <- c("nREP", "AUC", "P(WT)", "P(KPMT)", "adjP")
    topGSEAwt.lm[["Rank(P)"]] <- topGSEAwt.m
    tmp.s <- sort(topGSEAwt.m[, "AUC"], decreasing = TRUE, index.return = TRUE)
    topGSEAwt.lm[["Rank(AUC)"]] <- topGSEAwt.m[tmp.s$ix, ]
    topGSEAwt.lm[["Genestat"]] <- lapply(1:nrow(topGSEAwt.m), function(i) {
      pathway_name <- rownames(topGSEAwt.m)[i]
      EID.v <- intersect(ptw.ls[[pathway_name]], rankEID.v)
      pathgene.m <- matrix(NA, nrow = length(EID.v), ncol = 2)
      rownames(pathgene.m) <- sym.v[match(EID.v, rankEID.v)]
      colnames(pathgene.m) <- c("Pvalue", "Statistic")
      pathgene.m[, ] <- rankEID.m[match(EID.v, rankEID.v), ]
      pathgene.m
    })
    names(topGSEAwt.lm[["Genestat"]]) <- rownames(topGSEAwt.m)
  } else {
    stop("Error: No significant pathways identified after filtering.")
  }
  
  message("Function execution complete.")
  return(topGSEAwt.lm)
}

gseaWTfn <- function(termEID.v,rankEID.v,minN=5){
  commonEID.v <- intersect(termEID.v,rankEID.v);
  nrep <- length(commonEID.v);
  if(length(commonEID.v)>=minN){
    otherEID.v <- setdiff(rankEID.v,termEID.v);
    match(commonEID.v,rankEID.v) -> rank1.idx;
    match(otherEID.v,rankEID.v) -> rank2.idx;
    wt.o <- wilcox.test(rank1.idx,rank2.idx,alt="less");
    pv <- wt.o$p.value;
    n1 <- length(rank1.idx);
    n2 <- length(rank2.idx);
    auc <- 1 - wt.o$stat/(n1*n2);
    ### now do kpmt
    pop.v <- 1:length(rankEID.v);
    names(pop.v) <- rankEID.v;
    obs.v <- commonEID.v;
    pvKPMT <- kpmt(pop=pop.v,obs=obs.v,tail="lower")[[4]];
    out.v <- c(nrep,auc,pv,pvKPMT);
  }else {
    out.v <- c(nrep,0,1,1);
  }
  return(out.v);
}
