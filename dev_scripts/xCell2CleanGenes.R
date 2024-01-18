#' xCell2PrepData function
#'
#' This function ...
#'
#'@importFrom limma normalizeBetweenArrays
#'@importFrom Matrix colSums Diagonal
#' @param ref Reference gene expression matrix.
#' @param mix Mixture gene expression matrix.
#' @param ref_type .
#' @param val_type .
#' @param top_var_genes Select top variable genes using Seurat's FindVariableFeatures function (only for scRNA-seq data)
#' @param n_var_genes Number of top variable genes to select  (only for scRNA-seq data)
#' @return A list of ref and mix with the shared genes after cleaning
#' @export
xCell2PrepData <- function(ref, mix, ref_type, val_type, top_var_genes = FALSE, n_var_genes = 5000){

  if (!ref_type %in% c("rnaseq", "array", "sc")) {
    stop("ref_type should be 'rnaseq', 'array' or 'sc'.")
  }

  if (!val_type %in% c("rnaseq", "array")) {
    stop("val_type should be 'rnaseq' or 'array'.")
  }



  if (ref_type == "sc") {
    message("Normalizing scRNA-Seq reference to log2-CPM.")
    # TODO: For non 10X also do TPM
    lib_sizes <- Matrix::colSums(ref)
    norm_factor <- 1e6 / lib_sizes
    ref_norm <- ref %*% Matrix::Diagonal(x = norm_factor)
    colnames(ref_norm) <- colnames(ref)
    ref <- as.matrix(ref_norm) # TODO: Find a way to reduce memory usage by keeping matrix sparse
    ref <- log2(ref+1)
  }


  # TODO
  if (top_var_genes & ref_type == "sc") {

  }

  # Select shared genes
  shared_genes <- intersect(rownames(ref), rownames(mix))
  message(paste0(length(shared_genes), " genes are shared between reference and mixture."))
  ref <- ref[shared_genes,]
  mix <- mix[shared_genes,]


  # Data log-transformation
  if(max(ref) >= 50){
    ref <-log2(ref+1)
  }

  if(max(mix) >= 50){
    mix <-log2(mix+1)
  }


  y <- cbind(as.data.frame(ref), as.data.frame(mix))
  y <- limma::normalizeBetweenArrays(y)
  ref <- y[,1:ncol(ref)]
  mix <- y[,(ncol(ref)+1):ncol(y)]

  mix <- (2^mix)-1

  return(list("ref" = ref, "mix" = mix))
}
