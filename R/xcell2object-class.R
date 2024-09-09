#' xCell2Object Class
#'
#' An S4 class to represent the xCell2 reference object.
#'
#' @slot signatures List of xCell2 signatures.
#' @slot dependencies List of cell type dependencies.
#' @slot params Tibble of cell type linear transformation parameters.
#' @slot spill_mat Matrix of cell types spillover correction factors.
#' @slot genes_used Character vector of genes used to train the xCell2 reference object.
#' @importFrom methods new
#' @exportClass xCell2Object
setClass("xCell2Object", slots = list(
  signatures = "list",
  dependencies = "list",
  params = "data.frame",
  spill_mat = "matrix",
  genes_used = "character"
))

#' Get Signatures
#'
#' Access the \code{signatures} slot from an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @return A list containing the xCell2 signatures.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' getSignatures(DICE_demo.xCell2Ref)
setGeneric("getSignatures", function(object) standardGeneric("getSignatures"))

setMethod("getSignatures", "xCell2Object", function(object) {
  return(object@signatures)
})

#' Set Signatures
#'
#' Set the \code{signatures} slot in an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @param value A list to set in the \code{signatures} slot.
#' @return The modified \code{xCell2Object}.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' setSignatures(DICE_demo.xCell2Ref) <- list("B cell" = c("MS4A1", "IGHM", "FCRL1", "IGHD"))
setGeneric("setSignatures<-", function(object, value) standardGeneric("setSignatures<-"))

setMethod("setSignatures<-", "xCell2Object", function(object, value) {
  object@signatures <- value
  return(object)
})

#' Get Dependencies
#'
#' Access the \code{dependencies} slot from an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @return A list containing the cell type dependencies.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' getDeps(DICE_demo.xCell2Ref)
setGeneric("getDeps", function(object) standardGeneric("getDeps"))

setMethod("getDeps", "xCell2Object", function(object) {
  return(object@dependencies)
})

#' Set Dependencies
#'
#' Set the \code{dependencies} slot in an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @param value A list to set in the \code{dependencies} slot.
#' @return The modified \code{xCell2Object}.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' setDeps(DICE_demo.xCell2Ref) <-
#'  list("T cells, CD4+" = list("descendants" = "T cells, CD4+, memory",  "ancestors" = "T cells"))
setGeneric("setDeps<-", function(object, value) standardGeneric("setDeps<-"))

setMethod("setDeps<-", "xCell2Object", function(object, value) {
  object@dependencies <- value
  return(object)
})

#' Get Params
#'
#' Access the \code{params} slot from an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @return A tibble containing the linear transformation parameters.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' getParams(DICE_demo.xCell2Ref)
setGeneric("getParams", function(object) standardGeneric("getParams"))

setMethod("getParams", "xCell2Object", function(object) {
  return(object@params)
})

#' Set Params
#'
#' Set the \code{params} slot in an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @param value A data frame/tibble to set in the \code{params} slot.
#' @return The modified \code{xCell2Object}.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' setParams(DICE_demo.xCell2Ref) <-
#'  data.frame(celltype = "B cells", a = 0.856, b = 0.197, m = 1.85, n = 0.00185)
setGeneric("setParams<-", function(object, value) standardGeneric("setParams<-"))

setMethod("setParams<-", "xCell2Object", function(object, value) {
  object@params <- value
  return(object)
})

#' Get Spillover Matrix
#'
#' Access the \code{spill_mat} slot from an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @return A matrix representing the spillover correction factors.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' getSpillMat(DICE_demo.xCell2Ref)
setGeneric("getSpillMat", function(object) standardGeneric("getSpillMat"))

setMethod("getSpillMat", "xCell2Object", function(object) {
  return(object@spill_mat)
})

#' Set Spillover Matrix
#'
#' Set the \code{spill_mat} slot in an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @param value A matrix to set in the \code{spill_mat} slot.
#' @return The modified \code{xCell2Object}.
#' @export
setGeneric("setSpillMat<-", function(object, value) standardGeneric("setSpillMat<-"))

setMethod("setSpillMat<-", "xCell2Object", function(object, value) {
  object@spill_mat <- value
  return(object)
})

#' Get Genes Used
#'
#' Access the \code{genes_used} slot from an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @return A character vector of gene names used to train the xCell2 reference object
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' getGenesUsed(DICE_demo.xCell2Ref)
setGeneric("getGenesUsed", function(object) standardGeneric("getGenesUsed"))

setMethod("getGenesUsed", "xCell2Object", function(object) {
  return(object@genes_used)
})

#' Set Genes Used
#'
#' Set the \code{genes_used} slot in an \code{xCell2Object}.
#'
#' @param object An \code{xCell2Object}.
#' @param value A character vector to set in the \code{genes_used} slot.
#' @return The modified \code{xCell2Object}.
#' @export
#' @examples
#' data("DICE_demo.xCell2Ref")
#' setGenesUsed(DICE_demo.xCell2Ref) <- c("GeneA", "GeneB", "GeneC")
setGeneric("setGenesUsed<-", function(object, value) standardGeneric("setGenesUsed<-"))

setMethod("setGenesUsed<-", "xCell2Object", function(object, value) {
  object@genes_used <- value
  return(object)
})