#' @title Access Cell Type Signatures
#' @description Retrieve or assign the cell type-specific gene signatures for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of cell type-specific gene signatures (for the setter).
#' @return For `getSignatures`, a list of cell type-specific gene signatures. 
#'         For `setSignatures<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @aliases getSignatures setSignatures<-
#' @rdname signatures
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setSignatures(DICE_demo.xCell2Ref) <- list("T_cells" = c("GeneA", "GeneB"), "B_cells" = c("GeneC"))
#' getSignatures(DICE_demo.xCell2Ref)
#' @export
setMethod("getSignatures", "xCell2Object", function(object) object@signatures)

#' @rdname signatures
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setSignatures(DICE_demo.xCell2Ref) <- list("T_cells" = c("GeneA", "GeneB"))
#' @export
setReplaceMethod("setSignatures", "xCell2Object", function(object, value) {
  object@signatures <- value
  object
})

#' @title Access Cell Type Dependencies
#' @description Retrieve or assign the hierarchical dependencies between cell types for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of hierarchical dependencies (for the setter).
#' @return For `getDeps`, a list of hierarchical dependencies. 
#'         For `setDeps<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @aliases getDeps setDeps<-
#' @rdname dependencies
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setDeps(DICE_demo.xCell2Ref) <- list("Parent" = "ChildType1")
#' getDeps(DICE_demo.xCell2Ref)
#' @export
setMethod("getDeps", "xCell2Object", function(object) object@dependencies)

#' @rdname dependencies
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setDeps(DICE_demo.xCell2Ref) <- list("ParentType" = "ChildType1")
#' @export
setReplaceMethod("setDeps", "xCell2Object", function(object, value) {
  object@dependencies <- value
  object
})

#' @title Access Transformation Parameters
#' @description Retrieve or assign linear transformation parameters for cell types for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A data frame of transformation parameters (for the setter).
#' @return For `getParams`, a data frame of transformation parameters. 
#'         For `setParams<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @aliases getParams setParams<-
#' @rdname params
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setParams(DICE_demo.xCell2Ref) <- data.frame(celltype = "T_cells", a = 0.5, b = 2, m = 1.5)
#' getParams(DICE_demo.xCell2Ref)
#' @export
setMethod("getParams", "xCell2Object", function(object) object@params)

#' @rdname params
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setParams(DICE_demo.xCell2Ref) <- data.frame(celltype = "B_cells", a = 0.6, b = 1.8, m = 2.1)
#' @export
setReplaceMethod("setParams", "xCell2Object", function(object, value) {
  object@params <- value
  object
})

#' @title Access Spillover Matrix
#' @description Retrieve or assign the spillover correction matrix for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A matrix of spillover correction factors (for the setter).
#' @return For `getSpillMat`, a matrix of spillover correction factors. 
#'         For `setSpillMat<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @aliases getSpillMat setSpillMat<-
#' @rdname spillMat
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' spill_mat <- matrix(c(1, 0.1, 0.1, 1), nrow = 2, byrow = TRUE)
#' rownames(spill_mat) <- colnames(spill_mat) <- c("T_cells", "B_cells")
#' setSpillMat(DICE_demo.xCell2Ref) <- spill_mat
#' getSpillMat(DICE_demo.xCell2Ref)
#' @export
setMethod("getSpillMat", "xCell2Object", function(object) object@spill_mat)

#' @rdname spillMat
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' spill_mat <- matrix(c(1, 0.05, 0.05, 1), nrow = 2, byrow = TRUE)
#' rownames(spill_mat) <- colnames(spill_mat) <- c("T_cells", "B_cells")
#' setSpillMat(DICE_demo.xCell2Ref) <- spill_mat
#' @export
setReplaceMethod("setSpillMat", "xCell2Object", function(object, value) {
  object@spill_mat <- value
  object
})

#' @title Access Genes Used
#' @description Retrieve or assign the genes used in training the reference for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A character vector of genes (for the setter).
#' @return For `getGenesUsed`, a character vector of genes. 
#'         For `setGenesUsed<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @aliases getGenesUsed setGenesUsed<-
#' @rdname genesUsed
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setGenesUsed(DICE_demo.xCell2Ref) <- c("GeneA", "GeneB", "GeneC")
#' getGenesUsed(DICE_demo.xCell2Ref)
#' @export
setMethod("getGenesUsed", "xCell2Object", function(object) object@genes_used)

#' @rdname genesUsed
#' @examples
#' data(DICE_demo.xCell2Ref, package = "xCell2")
#' setGenesUsed(DICE_demo.xCell2Ref) <- c("GeneA", "GeneB", "GeneC")
#' @export
setReplaceMethod("setGenesUsed", "xCell2Object", function(object, value) {
  object@genes_used <- value
  object
})