library(tidyverse)

# 1. References: make a list that store all human/mouse references  -------

# >>> Each reference prepared individually in prep_ref.R

# Load human cell types references
human.refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))

human.refsRDSList <- lapply(human.refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})
saveRDS(human.refsRDSList, "/bigdata/almogangel/xCell2_dev/paper_figures/data/human_references.rds")

# Load mouse cell types references
mouse.refList <- list(rna_seq = c(mixed = "mouse_rnaseq_data"),
                array = c("igd"),
                sc = c(mixed = "tm_blood"))

mouse.refsRDSList <- lapply(mouse.refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})
saveRDS(mouse.refsRDSList, "/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_references.rds")


# 2. Validations: make RDS files that store all human/mouse validation datasets -------

# >>> Each validation dataset prepared individually in prep_vals.R

loadVals <- function(valList,
                     valMixDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/",
                     valTruthDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/",
                     conversion_mat = celltype_conversion){
  
  # load mixtures
  valMixtures <- lapply(valList, function(tissue){
    mixtures <- lapply(tissue, function(val){
      as.matrix(read.table(paste0(valMixDir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE, sep = "\t"))
    })
    names(mixtures) <- tissue
    mixtures
  })
  
  # Load truth
  valTruth <- lapply(valList, function(tissue){
    truths <- lapply(tissue, function(val){
      truth <- as.matrix(read.table(paste0(valTruthDir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
      
      # Change cell types labels
      rownames(truth) <- plyr::mapvalues(rownames(truth), conversion_mat$all_labels, conversion_mat$xCell2_labels, warn_missing = FALSE)
      rows <- rownames(truth)
      truth <- suppressWarnings(apply(truth, 2, as.numeric))
      rownames(truth) <- rows
      
      # Remove (10e9/l) values
      truth <- truth[!grepl("/l)", rownames(truth)),]
      
      # Remove values with all NAs
      truth <- truth[!apply(truth, 1, function(x){all(is.na(x))}),]
      
      # Remove values > 100
      truth <- truth[!apply(truth, 1, function(x){max(x, na.rm = TRUE) > 100}),]
      
      # Scale fractions
      if (max(truth, na.rm = TRUE) > 1) {
        truth <- truth/100
      }

      truth
    })
    names(truths) <- tissue
    truths
  })
  
  return(list(mixtures = valMixtures, truth = valTruth))
}

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_dev/paper_figures/data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Load human validation datasets
human.vals.list <- list(blood = c("BG_blood", "GSE107011", "GSE107572", "GSE127813",
                                 "SDY311", "SDY420", "SDY67", "GSE64385", "GSE65133", "GSE106898", "GSE64655", "GSE59654",
                                 "GSE107990", "GSE20300", "GSE65135", "GSE77343", "GSE77344", "GSE130824"),
                       tumor = c("ccRCC_cytof_CD45+", "NSCLC_cytof", "WU_ccRCC_RCCTC"),
                       other = c("GSE120444", "GSE115823", "GSE93722"))

human.vals <- loadVals(valList = human.vals.list)
saveRDS(human.vals, "/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds")

# Load mouse validation datasets
mouse.vals.list <- list(other = c("chen", "petitprez"))
mouse.vals <- loadVals(valList = mouse.vals.list,
                            valMixDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/mixtures/mouse/",
                            valTruthDir = "/bigdata/almogangel/xCell2_data/benchmarking_data/validations/truth/mouse/")
saveRDS(mouse.vals, "/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_validation.rds")
# 3. Match reference-validation pairs by shared cell types -------

combineRefVal <- function(vals, val2type, refList, splitDependencies){
  
  # Tibble of validation datasets
  vals.tbl <- enframe(vals$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(vals$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))
  
  # Tibble of cell type references
  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")
  
  types <- unique(vals.tbl$val_type)
  
  # Generate a tibble with validations-references combinations
  combined.tbl <- lapply(types, function(type){
    val_tmp <- vals.tbl[vals.tbl$val_type %in% type,]
    type <- if(type == "other") "mixed" else c(type, "mixed")
    ref_tmp <- refs.tbl[refs.tbl$ref_tissue %in% type,]
    crossing(val_tmp, ref_tmp)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::select(-mixture, -truth)
  rm(vals.tbl, refs.tbl)
  
  combined.tbl <- combined.tbl %>%
    mutate(val_data_type = val2type[val_dataset])
  
  
  # Load refs dependencies
  refs_dir <- "/bigdata/almogangel/xCell2_data/benchmarking_data/"
  getDependencies <- function(lineage_file_checked){
    ont <- readr::read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
      mutate_all(as.character)
    
    celltypes <- pull(ont[,2])
    celltypes <- gsub("_", "-", celltypes)
    dep_list <- vector(mode = "list", length = length(celltypes))
    names(dep_list) <- celltypes
    
    for (i in 1:nrow(ont)) {
      descendants <-  gsub("_", "-", strsplit(pull(ont[i,3]), ";")[[1]])
      descendants <- descendants[!is.na(descendants)]
      
      ancestors <-  gsub("_", "-", strsplit(pull(ont[i,4]), ";")[[1]])
      ancestors <- ancestors[!is.na(ancestors)]
      
      dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)
      
    }
    
    return(dep_list)
  }
  
  refsDepList <- lapply(refList, function(ref_type){
    refsDep <- lapply(ref_type, function(ref){
      
      # Load reference
      ref <- readRDS(paste0(refs_dir, "references/", ref, "_ref.rds"))
      
      # Get cell types dependencies
      dep_list <- getDependencies(ref$lineage_file)
      
      dep_cts <- sapply(names(dep_list), function(d){
        deps <- unname(unlist(dep_list[[d]]))
      })
      dep_cts[order(sapply(dep_cts, length), decreasing = T)]
    })
    names(refsDep) <- ref_type
    refsDep
  })

  # Split references with dependent cell types for deconvolution methods
  if (!splitDependencies) {
    
    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      mutate(refsDeps = NA) %>%
      # In non deconvolution method celltype_classes = shared_celltypes
      mutate(celltype_classes = list(shared_celltypes)) %>%
      return()
    
  }else{
    
    splitDepCellTypes <- function(refsDeps){
      
      refsDeps <- refsDeps[order(sapply(refsDeps, length), decreasing = T)]
      
      classes <- list("1" = c())
      for (ct in names(refsDeps)) {
        
        classes <- classes[order(sapply(classes, length), decreasing = F)]
        deps <- refsDeps[[ct]]
        
        for (class in 1:length(classes)) {
          if (!any(deps %in% classes[[class]])) {
            classes[[class]] <- c(classes[[class]], ct)
            break
          }
        }
        
        if (!ct %in% unname(unlist(classes))) {
          classes[[paste0(length(classes)+1)]] <- ct
        }
      }
      
      one_element_list <- unname(lapply(classes, length) == 1)
      if (sum(one_element_list) > 0) {
        classes[!one_element_list][[1]] <- unname(c(unlist(classes[!one_element_list][1]), unlist(classes[one_element_list])))
        classes <- classes[!one_element_list]
      }
      
      
      return(classes)
    }
    
    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
      mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
      # Split dependent cell types in each CIBERSORTx run
      mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
      unnest(celltype_classes) %>%
      return()
    
  }
  
}


# Human
human.vals.names <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "SDY311",
              "SDY420", "SDY67", "GSE64385", "GSE65133", "GSE106898", "GSE64655",
              "GSE59654", "GSE107990", "GSE20300", "GSE65135", "GSE77343", "GSE77344",
              "ccRCC_cytof_CD45+", "NSCLC_cytof", "WU_ccRCC_RCCTC", "GSE120444", "GSE115823",
              "GSE93722", "GSE130824")
human.vals.type <- c("rnaseq", "rnaseq", "rnaseq", "rnaseq", "array",
               "array", "rnaseq", "array", "array", "array", "rnaseq",
               "array", "array", "array", "array", "array", "array",
               "rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq")
val2type.human <- setNames(human.vals.type, human.vals.names)

refList.human <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))

refs.vals.matched.human <- combineRefVal(vals = human.vals, val2type = val2type.human, refList = refList.human, splitDependencies = FALSE)
saveRDS(refs.vals.matched.human, "/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.human.rds")


# Match between mouse references and validation datasets
mouse.vals.names <- c("chen", "petitprez")
mouse.vals.type <- c("rnaseq", "rnaseq")
val2type.mouse <- setNames(mouse.vals.type, mouse.vals.names)

refList.mouse <- list(rna_seq = c(mixed = "mouse_rnaseq_data"),
                      array = c(mixed = "igd"),
                      sc = c(mixed = "tm_blood"))

refs.vals.matched.mouse <- combineRefVal(vals = mouse.vals, val2type = val2type.mouse, refList = refList.mouse, splitDependencies = FALSE)
saveRDS(refs.vals.matched.mouse, "/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.mouse.rds")