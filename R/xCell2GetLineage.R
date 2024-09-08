#' xCell2GetLineage function
#'
#' This function identify the descendants and ancestors of each cell type based on the cell ontology tree.
#' If no output file is specified, the function returns a list of cell type dependencies.
#' If an output file is specified, the function writes the cell type dependencies to a TSV file.
#'
#' @importFrom ontoProc getCellOnto
#' @importFrom ontologyIndex get_descendants get_ancestors
#' @importFrom dplyr select mutate rowwise
#' @importFrom tibble as_tibble
#' @importFrom readr write_tsv
#' @param labels A data frame with four columns:
#'   "ont": the cell type ontology as a character (i.e., "CL:0000545" or NA if there is no ontology).
#'   "label": the cell type name as a character (i.e., "T-helper 1 cell").
#'   "sample": the cell type sample that match the column name in ref.
#'   "dataset": sample's source dataset or subject (can be the same for all sample if no such information).
#' @param outFile An optional output file name to write the cell type dependencies.
#' @return A list of cell type dependencies, or a TSV file containing the cell type dependencies if an output file is specified.
#' @export
xCell2GetLineage <- function(labels, outFile = NULL) {

  labelsUniq <- labels %>%
    tibble::as_tibble() %>%
    dplyr::select("ont", "label") %>%
    unique()

  if (all(is.na(labelsUniq[, 1]))) {
    message("Cannot find cell types dependencies - no ontologies provided")
    lineageOut <- labelsUniq %>%
      dplyr::mutate(descendants = "", ancestors = "")
  } else {
    cl <- ontoProc::getCellOnto()
    labelsUniq$descendants <- NA
    labelsUniq$ancestors <- NA
    for (i in 1:nrow(labelsUniq)) {
      ont <- labelsUniq$ont[i]
      descendants <- ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)
      ancestors <- ontologyIndex::get_ancestors(cl, terms = ont)
      ancestors <- ancestors[ancestors != ont]
      descendants <- paste(dplyr::pull(labelsUniq[dplyr::pull(labelsUniq[, 1]) %in% descendants, 2]), collapse = ";")
      ancestors <- paste(dplyr::pull(labelsUniq[dplyr::pull(labelsUniq[, 1]) %in% ancestors, 2]), collapse = ";")
      labelsUniq$descendants[i] <- descendants
      labelsUniq$ancestors[i] <- ancestors
    }
    lineageOut <- labelsUniq
  }

  # If no output is provided, xCell2GetLineage will return dependencies right away
  if (is.null(outFile)) {
    celltypes <- gsub("_", "-", dplyr::pull(lineageOut[, 2]))
    depList <- vector(mode = "list", length = length(celltypes))
    names(depList) <- celltypes

    for (i in seq_len(nrow(lineageOut))) {
      descendants <- gsub("_", "-", strsplit(dplyr::pull(lineageOut[i, 3]), ";")[[1]])
      descendants <- descendants[!is.na(descendants)]

      ancestors <- gsub("_", "-", strsplit(dplyr::pull(lineageOut[i, 4]), ";")[[1]])
      ancestors <- ancestors[!is.na(ancestors)]

      depList[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)
    }

    return(depList)
  } else {
    readr::write_tsv(lineageOut, outFile)
    warning("It is recommended that you manually check the cell type lineage file: ", outFile)
  }
}
