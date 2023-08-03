#' xCell2CleanGenes function
#'
#' This function clean and get shared genes between reference and mixture
#'
#' @param ref Reference gene expression matrix.
#' @param mix Mixture gene expression matrix.
#' @param gene_groups Gene group to clean c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
#' @return A list of ref and mix with the shared genes after cleaning
#' @export
xCell2CleanGenes <- function(ref, mix, gene_groups = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")){

  if (all(startsWith(rownames(ref), "ENSG"))) {
    message("Assuming genes are in Ensembl ID...")
    id <- 2
  }else{
    message("Assuming genes are in symbol ID...")
    id <- 3
  }

  # Remove gene groups
  ref <- ref[!rownames(ref) %in% hs.genelist[hs.genelist$gene_group %in% gene_groups, id],]

  # TODO: add option to choose only protein coding genes

  shared_genes <- intersect(rownames(ref), rownames(mix))
  message(paste0(shared_genes, " genes are shared between reference and mixture after cleaning."))

  ref <- ref[shared_genes,]
  mix <- mix[shared_genes,]

  return(list("ref" = ref, "mix" = mix))
}

