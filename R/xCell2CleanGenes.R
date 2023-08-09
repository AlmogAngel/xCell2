#' xCell2CleanGenes function
#'
#' This function clean and get shared genes between reference and mixture
#'
#'@importFrom Seurat CreateSeuratObject FindVariableFeatures VariableFeatures
#' @param ref Reference gene expression matrix.
#' @param mix Mixture gene expression matrix.
#' @param top_var_genes Select top variable genes using Seurat's FindVariableFeatures function (only for scRNA-seq data)
#' @param gene_groups Gene group to clean c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
#' @param use_protein_coding Boolean - use only protein coding genes?
#' @param n_var_genes Number of top variable genes to select  (only for scRNA-seq data)
#' @return A list of ref and mix with the shared genes after cleaning
#' @export
xCell2CleanGenes <- function(ref, mix, gene_groups = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), use_protein_coding = TRUE, top_var_genes = FALSE, n_var_genes = 10000){

  if (all(startsWith(rownames(ref), "ENSG"))) {
    message("Assuming genes type are Ensembl.")
    id <- 2
  }else{
    message("Assuming genes type are Symbol.")
    id <- 3
  }

  # Remove gene groups
  ngenes_before <- nrow(ref)
  ref <- ref[!rownames(ref) %in% hs.genelist[hs.genelist$gene_group %in% gene_groups, id],]

  # Use protein coding genes
  if (use_protein_coding) {
    pc_genes <- gencode.v22.broad.category[gencode.v22.broad.category$genes_type == "protein_coding", id]
    ref <- ref[rownames(ref) %in% pc_genes,]
  }

  if (top_var_genes) {

    message(paste0("Selecting top ", n_var_genes, " variable genes - make sure you are using scRNA-Seq reference data!"))
    genes_names <- rownames(ref)
    ref.srt <- Seurat::CreateSeuratObject(counts = ref)
    ref.srt <- Seurat::FindVariableFeatures(ref.srt, selection.method = "vst", nfeatures = n_var_genes, verbose = FALSE)
    var_genes <- Seurat::VariableFeatures(ref.srt)

    # Check if Seurat changed genes names
    if (!all(var_genes %in% genes_names)) {
      # TODO: find a solution
      # rownames(ref.norm) <- genes_names # Because Seurat change genes names from "_" to "-"
      stop("Seurat genes name error")
    }

    ref <- ref[var_genes,]
  }

  ngenes_after <- nrow(ref)
  ngenes_removed <- ngenes_before - ngenes_after
  message(paste0(ngenes_removed, " total genes removed from reference."))


  shared_genes <- intersect(rownames(ref), rownames(mix))
  message(paste0(length(shared_genes), " genes are shared between reference and mixture after cleaning."))

  ref <- ref[shared_genes,]
  mix <- mix[shared_genes,]

  return(list("ref" = ref, "mix" = mix))
}

