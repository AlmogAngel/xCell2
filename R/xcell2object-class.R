#' @title xCell2Object Class
#' @description
#' An S4 class to represent the xCell2 reference object. This object contains
#' cell type-specific gene signatures, hierarchical dependencies, linear transformation parameters,
#' spillover correction matrices, and genes used for training.
#'
#' @slot signatures A list of cell type-specific gene signatures.
#' @slot dependencies A list of hierarchical dependencies between cell types.
#' @slot params A data frame containing linear transformation parameters for cell types.
#' @slot spill_mat A matrix containing spillover correction factors for cell types.
#' @slot genes_used A character vector of genes used for training the xCell2 reference object.
#' @seealso
#' \link[xCell2]{xcell2object-methods}
#'
#' @name xCell2Object-class
#' @rdname xCell2Object-class
#' @docType class
#' @exportClass xCell2Object
setClass(
  "xCell2Object",
  slots = list(
    signatures = "list",
    dependencies = "list",
    params = "data.frame",
    spill_mat = "matrix",
    genes_used = "character"
  ),
  prototype = list(
    signatures = list(),
    dependencies = list(),
    params = data.frame(),
    spill_mat = matrix(),
    genes_used = character()
  )
)