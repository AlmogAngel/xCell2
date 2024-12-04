#' @title Get Signatures
#' @description Retrieve the cell type-specific gene signatures from an \linkS4class{xCell2Object}.
#' @export
setGeneric("getSignatures", function(object) standardGeneric("getSignatures"))

#' @title Set Signatures
#' @description Assign cell type-specific gene signatures to an \linkS4class{xCell2Object}.
#' @export
setGeneric("setSignatures<-", function(object, value) standardGeneric("setSignatures<-"))

#' @title Get Dependencies
#' @description Retrieve the hierarchical dependencies between cell types from an \linkS4class{xCell2Object}.
#' @export
setGeneric("getDeps", function(object) standardGeneric("getDeps"))

#' @title Set Dependencies
#' @description Assign the hierarchical dependencies between cell types to an \linkS4class{xCell2Object}.
#' @export
setGeneric("setDeps<-", function(object, value) standardGeneric("setDeps<-"))

#' @title Get Parameters
#' @description Retrieve the linear transformation parameters for cell types from an \linkS4class{xCell2Object}.
#' @export
setGeneric("getParams", function(object) standardGeneric("getParams"))

#' @title Set Parameters
#' @description Assign the linear transformation parameters for cell types to an \linkS4class{xCell2Object}.
#' @export
setGeneric("setParams<-", function(object, value) standardGeneric("setParams<-"))

#' @title Get Spillover Matrix
#' @description Retrieve the spillover correction matrix from an \linkS4class{xCell2Object}.
#' @export
setGeneric("getSpillMat", function(object) standardGeneric("getSpillMat"))

#' @title Set Spillover Matrix
#' @description Assign the spillover correction matrix to an \linkS4class{xCell2Object}.
#' @export
setGeneric("setSpillMat<-", function(object, value) standardGeneric("setSpillMat<-"))

#' @title Get Genes Used
#' @description Retrieve the genes used in training the reference from an \linkS4class{xCell2Object}.
#' @export
setGeneric("getGenesUsed", function(object) standardGeneric("getGenesUsed"))

#' @title Set Genes Used
#' @description Assign the genes used in training the reference to an \linkS4class{xCell2Object}.
#' @export
setGeneric("setGenesUsed<-", function(object, value) standardGeneric("setGenesUsed<-"))