#' @title Access Cell Type Signatures
#' @description Retrieve or assign the cell type-specific gene signatures for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of cell type-specific gene signatures (for the setter).
#' @return For `getSignatures`, a list of cell type-specific gene signatures. 
#'         For `setSignatures<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @name signatures
#' @rdname signatures
#' @export
setGeneric("getSignatures", function(object) standardGeneric("getSignatures"))

#' @rdname signatures
#' @export
setGeneric("setSignatures<-", function(object, value) standardGeneric("setSignatures<-"))

#' @title Access Cell Type Dependencies
#' @description Retrieve or assign hierarchical dependencies between cell types for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A list of hierarchical dependencies (for the setter).
#' @return For `getDeps`, a list of hierarchical dependencies. 
#'         For `setDeps<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @name dependencies
#' @rdname dependencies
#' @export
setGeneric("getDeps", function(object) standardGeneric("getDeps"))

#' @rdname dependencies
#' @export
setGeneric("setDeps<-", function(object, value) standardGeneric("setDeps<-"))

#' @title Access Transformation Parameters
#' @description Retrieve or assign linear transformation parameters for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A data frame of transformation parameters (for the setter).
#' @return For `getParams`, a data frame of transformation parameters. 
#'         For `setParams<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @name params
#' @rdname params
#' @export
setGeneric("getParams", function(object) standardGeneric("getParams"))

#' @rdname params
#' @export
setGeneric("setParams<-", function(object, value) standardGeneric("setParams<-"))

#' @title Access Spillover Matrix
#' @description Retrieve or assign the spillover correction matrix for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A matrix of spillover correction factors (for the setter).
#' @return For `getSpillMat`, a matrix of spillover correction factors. 
#'         For `setSpillMat<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @name spillMat
#' @rdname spillMat
#' @export
setGeneric("getSpillMat", function(object) standardGeneric("getSpillMat"))

#' @rdname spillMat
#' @export
setGeneric("setSpillMat<-", function(object, value) standardGeneric("setSpillMat<-"))

#' @title Access Genes Used
#' @description Retrieve or assign the genes used in training the reference for an \linkS4class{xCell2Object}.
#' @param object An \linkS4class{xCell2Object}.
#' @param value A character vector of genes (for the setter).
#' @return For `getGenesUsed`, a character vector of genes. 
#'         For `setGenesUsed<-`, the updated \linkS4class{xCell2Object}.
#' @seealso \link{xCell2Object-class}
#' @name genesUsed
#' @rdname genesUsed
#' @export
setGeneric("getGenesUsed", function(object) standardGeneric("getGenesUsed"))

#' @rdname genesUsed
#' @export
setGeneric("setGenesUsed<-", function(object, value) standardGeneric("setGenesUsed<-"))