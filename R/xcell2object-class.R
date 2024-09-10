#' xCell2Object Class
#'
#' An S4 class to represent the xCell2 reference object.
#'
#' @slot signatures List of xCell2 signatures.
#' @slot dependencies List of cell type dependencies.
#' @slot params Data frame of cell type linear transformation parameters.
#' @slot spill_mat Matrix of cell types spillover correction factors.
#' @slot genes_used Character vector of genes used to train the xCell2 reference object.
#' @importFrom methods setClass
#' @name xCell2Object-class
#' @rdname xCell2Object-class
#' @concept objects
#' @exportClass xCell2Object
setClass("xCell2Object", slots = list(
  signatures = "list",
  dependencies = "list",
  params = "data.frame",
  spill_mat = "matrix",
  genes_used = "character"
))
