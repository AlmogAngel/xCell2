#' Identify Cell Type Lineage Dependencies
#'
#' Identifies cell type dependencies based on the Cell Ontology, including both descendants and ancestors for each cell type.  
#' Enables manual inspection and refinement of lineage relationships to improve biological accuracy in \code{xCell2} analyses.
#'
#' @importFrom magrittr %>%
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom ontologyIndex get_descendants get_ancestors
#' @importFrom dplyr select mutate rowwise
#' @importFrom tibble as_tibble
#' @importFrom readr write_tsv
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom SingleCellExperiment colData
#'
#' @param labels A data frame with the following required columns:
#'   \itemize{
#'     \item \code{"ont"}: Cell type ontology ID (e.g., \code{"CL:0000545"}). Use \code{NA} if unavailable.  
#'       Ontologies can be accessed via \href{https://www.ebi.ac.uk/ols4/ontologies/cl}{EBI Ontology Lookup Service (OLS)}  
#'       or the \link[ontologyIndex]{ontologyIndex} package.
#'     \item \code{"label"}: Cell type name (e.g., \code{"T-helper 1 cell"}).
#'     \item \code{"sample"}: Sample or cell identifier matching column names in the gene expression matrix.
#'     \item \code{"dataset"}: Dataset or subject source. Use a constant value if not applicable.
#'   }
#' @param outFile Optional. Output file name for saving dependencies as a TSV file.  
#'   The file includes columns for \code{"ont"}, \code{"label"}, \code{"descendants"}, and \code{"ancestors"}.  
#'   Suitable for manual inspection and refinement before use in downstream analyses.
#'
#' @return 
#'   If \code{outFile} is:
#'   \itemize{
#'     \item \code{NULL}: Returns a list of dependencies for each cell type, with descendants and ancestors as components.
#'     \item Specified: Writes a TSV file and warns the user to inspect and validate results manually.
#'   }
#'
#' @details
#' The \code{xCell2GetLineage} function generates lineage relationships for cell types based on the Cell Ontology.  
#' These relationships refine lineage-based dependencies, improving the biological relevance of gene signatures.  
#' Users can:
#' \itemize{
#'   \item Use the generated TSV file for manual adjustments before training custom references via \code{\link{xCell2Train}}.
#'   \item Skip this step entirely, allowing \code{xCell2Train} to infer dependencies automatically.
#' }
#'
#' If no ontology IDs (\code{"ont"}) are provided, the function outputs empty dependencies with a message for user guidance.
#'
#' \strong{Relationship with Other Functions:}
#' \itemize{
#'   \item \code{\link{xCell2Train}}: Incorporates lineage relationships during reference training.
#'   \item \code{\link{xCell2Analysis}}: Uses trained references for enrichment analysis.
#' }
#'
#' @examples
#' # For detailed examples, see the xCell2 vignette.
#'
#' library(xCell2)
#'
#' # Load demo reference object
#' data(dice_demo_ref, package = "xCell2")
#'
#' # Prepare labels data frame
#' dice_labels <- SummarizedExperiment::colData(dice_demo_ref)
#' dice_labels <- as.data.frame(dice_labels)
#' dice_labels$ont <- NA
#' dice_labels$sample <- colnames(dice_demo_ref)
#' dice_labels$dataset <- "DICE"
#'
#' # Assign ontology IDs
#' dice_labels[dice_labels$label == "B cells", ]$ont <- "CL:0000236"
#' dice_labels[dice_labels$label == "Monocytes", ]$ont <- "CL:0000576"
#' dice_labels[dice_labels$label == "NK cells", ]$ont <- "CL:0000623"
#' dice_labels[dice_labels$label == "T cells, CD8+", ]$ont <- "CL:0000625"
#' dice_labels[dice_labels$label == "T cells, CD4+", ]$ont <- "CL:0000624"
#' dice_labels[dice_labels$label == "T cells, CD4+, memory", ]$ont <- "CL:0000897"
#'
#' # Generate cell type lineage dependencies
#' xCell2::xCell2GetLineage(labels = dice_labels)
#'
#' # Manually inspect and adjust saved dependencies for refined lineage relationships
#' # Use the adjusted file as input to xCell2Train via the `lineageFile` parameter.
#'
#' @seealso 
#' \code{\link{xCell2Train}} for training custom references with lineage data.  
#' \code{\link{xCell2Analysis}} for enrichment analysis using trained references.  
#' \code{\link{AnnotationHub}} to access ontology data.  
#' \code{\link{ontologyIndex}} to programmatically explore ontologies.
#'
#' @author Almog Angel and Dvir Aran
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
    ah <- AnnotationHub::AnnotationHub()
    # AH111554 cellOnto_2023.02.15 (2023)
    cl <- ah[["AH111554"]]
    labelsUniq$descendants <- NA
    labelsUniq$ancestors <- NA
    for (i in seq_len(nrow(labelsUniq))) {
      ont <- labelsUniq$ont[i]
      descendants <- ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)
      ancestors <- ontologyIndex::get_ancestors(cl, terms = ont)
      ancestors <- ancestors[ancestors != ont]
      descendants <- paste(dplyr::pull(labelsUniq[dplyr::pull(
        labelsUniq[, 1]
      ) %in% descendants, 2]), collapse = ";")
      ancestors <- paste(dplyr::pull(labelsUniq[dplyr::pull(
        labelsUniq[, 1]
      ) %in% ancestors, 2]), collapse = ";")
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
