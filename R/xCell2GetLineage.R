#########################################################################################
# Generate cell-type lineage automatically for dependencies
#########################################################################################

# USAGE:
# (1) labels - a data frame with:
#   (a) first column for cell type onthology - named "ont" (character)
#   (b) second column for cell type name - named "label" (character)
# (2) out_file - path to cell types lineage file for manual check (if NULL dependencies list will be made automatically)

# DEPENDENCIES:
# tidyverse, ontoProc, ontologyIndex

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
      dep_cells <- c(strsplit(pull(lineage.out[i,3]), ";")[[1]], strsplit(pull(lineage.out[i,4]), ";")[[1]])
      dep_cells <- gsub("_", "-", dep_cells)
      dep_list[[i]] <- dep_cells[!is.na(dep_cells)]
    }

    return(dep_list)
  }else{
    write_tsv(lineage.out, out_file)
    warning("It is recommended that you manually check the cell-type ontology file: ", out_file)
  }

}
