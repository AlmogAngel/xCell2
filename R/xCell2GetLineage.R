#' xCell2GetLineage function
#'
#' This function uses the `ontoProc` and `ontologyIndex` packages to identify the descendants and ancestors of each cell type based on the cell ontology tree. If no output file is specified, the function returns a list of cell type dependencies as a list. If an output file is specified, the function writes the cell type dependencies to a TSV file.
#'
#' @importFrom ontoProc
#' @importFrom ontologyIndex get_descendants
#' @importFrom ontologyIndex get_ancestors
#' @import ontoProc
#' @import dplyr
#' @import tibble
#' @param labels A data frame containing two columns: first column named "ont" for cell type onthology (character) and second column named "label" for cell type labels (character)
#' @param out_file An optional output file name to write the cell type dependencies to.
#' @return A list of cell type dependencies, or a TSV file containing the cell type dependencies if an output file is specified.
#' @export
xCell2GetLineage <- function(labels, out_file){

  cl <- ontoProc::getCellOnto()

  labels_uniq <- labels %>%
    as_tibble() %>%
    select(ont, label) %>%
    unique()

  lineage.out <- labels_uniq %>%
    rowwise() %>%
    mutate(descendants = list(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)),
           ancestors = list(ontologyIndex::get_ancestors(cl, terms = ont))) %>%
    mutate(ancestors = list(ancestors[ancestors != ont])) %>%
    mutate(descendants = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% descendants, 2]), collapse = ";"),
           ancestors = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% ancestors, 2]), collapse = ";"))

  # If not output provided xCell2GetLineage will return dependencies right away
  if (is.null(out_file)) {
    celltypes <- pull(lineage.out[,2])
    celltypes <- gsub("_", "-", celltypes)
    dep_list <- vector(mode = "list", length = length(celltypes))
    names(dep_list) <- celltypes

    for (i in 1:nrow(lineage.out)) {
      descendants <-  gsub("_", "-", strsplit(pull(lineage.out[i,3]), ";")[[1]])
      descendants <- descendants[!is.na(descendants)]

      ancestors <-  gsub("_", "-", strsplit(pull(lineage.out[i,4]), ";")[[1]])
      ancestors <- ancestors[!is.na(ancestors)]

      dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)
    }

    return(dep_list)
  }else{
    write_tsv(lineage.out, out_file)
    warning("It is recommended that you manually check the cell-type ontology file: ", out_file)
  }

}
