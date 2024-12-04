#' @title Get Cell Type Signatures
#' @description Retrieve the cell type-specific gene signatures from an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @return A list of cell type-specific gene signatures.
#' @examples
#' obj <- new("xCell2Object")
#' getSignatures(obj)
#' @name getSignatures
#' @rdname getSignatures-xCell2Object-method
#' @export
setMethod("getSignatures", "xCell2Object", function(object) object@signatures)

#' @title Set Cell Type Signatures
#' @description Assign cell type-specific gene signatures to an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of cell type-specific gene signatures.
#' @return The modified \linkS4class{xCell2Object}.
#' @examples
#' obj <- new("xCell2Object")
#' setSignatures(obj) <- list("Signature1" = c("GeneA", "GeneB"))
#' @name setSignatures
#' @rdname setSignatures-xCell2Object-method
#' @export
setReplaceMethod("setSignatures", "xCell2Object", function(object, value) {
  object@signatures <- value
  object
})

#' @title Get Cell Type Dependencies
#' @description Retrieve hierarchical dependencies between cell types from an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @return A list of cell type dependencies.
#' @examples
#' obj <- new("xCell2Object")
#' getDeps(obj)
#' @name getDeps
#' @rdname getDeps-xCell2Object-method
#' @export
setMethod("getDeps", "xCell2Object", function(object) object@dependencies)

#' @title Set Cell Type Dependencies
#' @description Assign hierarchical dependencies between cell types to an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of cell type dependencies.
#' @return The modified \linkS4class{xCell2Object}.
#' @examples
#' obj <- new("xCell2Object")
#' setDeps(obj) <- list("Dependency1" = "ParentType")
#' @name setDeps
#' @rdname setDeps-xCell2Object-method
#' @export
setReplaceMethod("setDeps", "xCell2Object", function(object, value) {
  object@dependencies <- value
  object
})

#' @title Get Transformation Parameters
#' @description Retrieve linear transformation parameters from an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @return A data frame of transformation parameters.
#' @examples
#' obj <- new("xCell2Object")
#' getParams(obj)
#' @name getParams
#' @rdname getParams-xCell2Object-method
#' @export
setMethod("getParams", "xCell2Object", function(object) object@params)

#' @title Set Transformation Parameters
#' @description Assign linear transformation parameters to an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A data frame of transformation parameters.
#' @return The modified \linkS4class{xCell2Object}.
#' @examples
#' obj <- new("xCell2Object")
#' setParams(obj) <- data.frame(param1 = 1, param2 = 2)
#' @name setParams
#' @rdname setParams-xCell2Object-method
#' @export
setReplaceMethod("setParams", "xCell2Object", function(object, value) {
  object@params <- value
  object
})

#' @title Get Spillover Matrix
#' @description Retrieve the spillover correction matrix from an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @return A matrix of spillover correction factors.
#' @examples
#' obj <- new("xCell2Object")
#' getSpillMat(obj)
#' @name getSpillMat
#' @rdname getSpillMat-xCell2Object-method
#' @export
setMethod("getSpillMat", "xCell2Object", function(object) object@spill_mat)

#' @title Set Spillover Matrix
#' @description Assign the spillover correction matrix to an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A matrix of spillover correction factors.
#' @return The modified \linkS4class{xCell2Object}.
#' @examples
#' obj <- new("xCell2Object")
#' setSpillMat(obj) <- matrix(0, nrow = 2, ncol = 2)
#' @name setSpillMat
#' @rdname setSpillMat-xCell2Object-method
#' @export
setReplaceMethod("setSpillMat", "xCell2Object", function(object, value) {
  object@spill_mat <- value
  object
})

#' @title Get Genes Used
#' @description Retrieve the genes used in training the reference from an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @return A character vector of genes.
#' @examples
#' obj <- new("xCell2Object")
#' getGenesUsed(obj)
#' @name getGenesUsed
#' @rdname getGenesUsed-xCell2Object-method
#' @export
setMethod("getGenesUsed", "xCell2Object", function(object) object@genes_used)

#' @title Set Genes Used
#' @description Assign the genes used in training the reference to an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A character vector of genes.
#' @return The modified \linkS4class{xCell2Object}.
#' @examples
#' obj <- new("xCell2Object")
#' setGenesUsed(obj) <- c("GeneA", "GeneB", "GeneC")
#' @name setGenesUsed
#' @rdname setGenesUsed-xCell2Object-method
#' @export
setReplaceMethod("setGenesUsed", "xCell2Object", function(object, value) {
  object@genes_used <- value
  object
})