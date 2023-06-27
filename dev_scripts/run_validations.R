library(tidyverse)
library(xCell2)
library(EPIC)
library(BayesPrism)

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")

# Load cell types labels conversion file
celltype_conversion <- read_tsv("celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Set references to use
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"))


# Function to load validation mixtures and truth values
loadVals <- function(valList, valMixDir, valTruthDir, conversion_mat = celltype_conversion){

  # load mixtures
  valMixtures <- lapply(valList, function(tissue){
    mixtures <- lapply(tissue, function(val){
      as.matrix(read.table(paste0(valMixDir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
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
      truth[is.na(truth)] <- 0
      truth[!duplicated(rows),]
    })
    names(truths) <- tissue
    truths
  })

  return(list(mixtures = valMixtures, truth = valTruth))
}


# Method functions -----------------------------------------------

# Get correlations using xCell2
getxCell2Cors <- function(valList, refList, celltype_conversion){

  print("Running xCell2 validations...")

  # Load xCell2 signatures
  print("Loading signatures...")
  sigsList <- pbapply::pbsapply(unname(unlist(refList)), function(ref){
    readRDS(paste0("references/xcell2_sigs/", ref, "_sigs.rds"))
  })

  # Generate a tibble with validations-references combinations
  vals.tbl <- enframe(valList$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(valList$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")

  vals.refs.tbl <- crossing(vals.tbl, refs.tbl)[1:3,] %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(truth)) %>%
    # Load xCell2 signatures
    mutate(xcell2sigs = list(sigsList[[ref_name]])) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(truth), unique(gsub("#.*", "", names(xcell2sigs@signatures)))))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Run xCell2
    mutate(xcell2res = list(xCell2Analysis(mixture, xcell2sigs)))

  # Calculate correlations
  getCors <- function(res_mat, truth_mat){
    samples <- intersect(colnames(res_mat), colnames(truth_mat))
    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))
    suppressWarnings(sapply(celltypes, function(ct){
      cor(res_mat[ct,samples], truth_mat[ct,samples], method = "spearman")
    }))
  }

  vals.refs.tbl %>%
    # In case some all cell types signatures are lost
    mutate(lost_celltypes = list(shared_celltypes[!shared_celltypes %in% rownames(xcell2res)])) %>%
    mutate(method = "xCell2", .before = everything()) %>%
    mutate(cor = list(getCors(xcell2res, truth))) %>%
    unnest_longer(cor, indices_to = "celltype") %>%
    return()

}

# Get correlations using CIBERSORTx
getCIBERSORTxCors <- function(valList, refList, celltype_conversion){

  print("Running CIBERSORTx validations...")

  runCIBERSORTx <- function(mix, refName, refType, celltypes2use , dir = "/bigdata/almogangel/CIBERSORTx_docker"){


    # Subset cell types from the signature matrix
    sigmat_tmp <- read.csv(paste0(dir, "/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
    sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
    write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


    # Make mixture file
    mix_tmp <- cbind("genes" = rownames(mix), mix)
    mix_file <- paste0(dir, "/mix-tmp.txt")
    write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


    # Run CIBERSORTx
    token <- "b72da36961922443b75a1b65beef27c0"

    # Make results directory
    results_dir <- paste0(dir, "/results")
    if (!dir.exists(results_dir)){
      dir.create(results_dir)
    }

    # Clean old results
    if(length(list.files(results_dir)) > 0){
      system(paste0("rm -f ", results_dir, "/*"))
    }

    single_cell <- ifelse(refType == "sc", TRUE, FALSE)

    cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                  token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell ," --rmbatchBmode ", !single_cell, " --rmbatchSmode ", single_cell,
                  " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

    # Run Docker via shell
    system(cmd, wait = TRUE)


    # Load results
    res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
    cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
    cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


    return(cibersortx_out)
  }
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


  # Load references cell types dependencies - to use when splitting dependent cell types before running CIBERSORTx
  print("Loading dependencies...")
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
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
  refsDepList <- pbapply::pblapply(refList, function(ref_type){
    refsDep <- lapply(ref_type, function(ref){

      # Load reference
      ref <- readRDS(paste0("references/", ref, "_ref.rds"))

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



  # Generate a tibble with validations-references combinations
  vals.tbl <- enframe(valList$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(valList$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")


  vals.refs.tbl <- crossing(vals.tbl, refs.tbl) %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(truth)) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(truth), names(refsDepList[[ref_type]][[ref_name]])))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least two cell types in the signature matrix
    filter(n_shared_celltypes > 1) %>%
    # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
    mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
    # Split dependent cell types in each CIBERSORTx run
    mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
    unnest(celltype_classes) %>%
    # Run CIBERSORTx
    rowwise() %>%
    mutate(cbrx_res = list(runCIBERSORTx(mix = mixture, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes)))



  # Calculate correlations
  getCors <- function(res_mat, truth_mat){
    samples <- intersect(colnames(res_mat), colnames(truth_mat))
    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))
    suppressWarnings(sapply(celltypes, function(ct){
      cor(res_mat[ct,samples], truth_mat[ct,samples], method = "spearman")
    }))
  }

  vals.refs.tbl %>%
    mutate(method = "CIBERSORTx", .before = everything()) %>%
    mutate(cor = list(getCors(cbrx_res, truth))) %>%
    unnest_longer(cor, indices_to = "celltype") %>%
    select(method:n_val_samples, n_shared_celltypes, celltype, cor) %>%
    return()

}

# Get correlations using EPIC
getEPICCors <- function(valList, refList, celltype_conversion){

  runEPIC <- function(mix, refsRDSList, refName, refType, celltypes2use, dir = "references"){

    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigGenes <- rownames(sigmat)

    # Subset GEP
    gep <- read.csv(paste0(dir, "/gep/", refName, "_gep.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    gep <- gep[,celltypes2use]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- ref.in$ref[,celltypeIndex]
    colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label


    # Generate reference for EPIC
    ref.raw.var <- sapply(unique(colnames(ref.raw)), function(ct){
      if(sum(colnames(ref.raw) == ct) > 1){
        apply(ref.raw[,colnames(ref.raw) == ct], 1, sd)
      }else{
        rep(0, length(ref.raw[,colnames(ref.raw) == ct]))
      }
    })
    ref.raw.var <- ref.raw.var[rownames(gep), colnames(gep)]

    epic_ref <- list("sigGenes" = sigGenes,
                     "refProfiles" = as.matrix(gep),
                     "refProfiles.var" = ref.raw.var)

    # Run EPIC
    res <- t(EPIC(bulk=mix, reference=epic_ref, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions)
    res <- res[rownames(res) != "otherCells",]

    return(res)
  }
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

  # Load references matrices
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  refsRDSList <- lapply(refList, function(ref_type){
    refs <- lapply(ref_type, function(ref){
      # Load reference
      ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
      ref.in
    })
    names(refs) <- ref_type
    refs
  })

  # Get cell types dependencies - to use when splitting dependent cell types before running EPIC
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
      ref <- readRDS(paste0("references/", ref, "_ref.rds"))

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


  # Generate a tibble with validations-references combinations
  vals.tbl <- enframe(valList$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(valList$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")


  vals.refs.tbl <- crossing(vals.tbl, refs.tbl) %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(truth)) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(truth), names(refsDepList[[ref_type]][[ref_name]])))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least two cell types in the signature matrix
    filter(n_shared_celltypes > 1) %>%
    # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
    mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
    # Split dependent cell types in each CIBERSORTx run
    mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
    unnest(celltype_classes) %>%
    # Run EPIC
    rowwise() %>%
    mutate(epic_res = list(runEPIC(mix = mixture, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes)))



  # Calculate correlations
  getCors <- function(res_mat, truth_mat){
    samples <- intersect(colnames(res_mat), colnames(truth_mat))
    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))
    suppressWarnings(sapply(celltypes, function(ct){
      cor(res_mat[ct,samples], truth_mat[ct,samples], method = "spearman")
    }))
  }

  vals.refs.tbl %>%
    mutate(method = "EPIC", .before = everything()) %>%
    mutate(cor = list(getCors(epic_res, truth))) %>%
    unnest_longer(cor, indices_to = "celltype") %>%
    select(method:n_val_samples, n_shared_celltypes, celltype, cor) %>%
    return()

}

# Get correlations using BayesPrism
getBayesPrismCors <- function(valList, refList, celltype_conversion){

  runBayesPrism <- function(mix, refsRDSList, refName, refType, celltypes2use, CPUs = 30, dir = "references/sigmats/"){


    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- t(ref.in$ref[,celltypeIndex])


    ref.filtered <- cleanup.genes(input=ref.raw,
                                     input.type="GEP",
                                     species="hs",
                                     gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"))


    ref.filtered.pc <-  select.gene.type(ref.filtered,
                                            gene.type = "protein_coding")


    # Find use signature genes/marker genes for bulk reference? ?

    labels <- ref.in$labels[celltypeIndex,]$label
    tumorKey <- NULL
    if ("malignant cell" %in% labels) {
      tumorKey <- "malignant cell"
    }

    myPrism <- new.prism(
      reference=ref.filtered.pc,
      mixture=t(mix),
      input.type="GEP",
      cell.type.labels = labels,
      cell.state.labels =labels,
      key=tumorKey)

    bp.res <- run.prism(prism = myPrism, n.cores = CPUs)

    res <- t(get.fraction (bp=bp.res,
                           which.theta="final",
                           state.or.type="type"))

    return(res)
  }
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

  # Load references matrices
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  refsRDSList <- lapply(refList, function(ref_type){
    refs <- lapply(ref_type, function(ref){
      # Load reference
      ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
      ref.in
    })
    names(refs) <- ref_type
    refs
  })

  # Get cell types dependencies - to use when splitting dependent cell types before running EPIC
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
  refsDepList <- pbapply::pblapply(refList, function(ref_type){
    refsDep <- lapply(ref_type, function(ref){

      # Load reference
      ref <- readRDS(paste0("references/", ref, "_ref.rds"))

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


  # Generate a tibble with validations-references combinations
  vals.tbl <- enframe(valList$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(valList$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")


  vals.refs.tbl <- crossing(vals.tbl, refs.tbl) %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(truth)) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(truth), names(refsDepList[[ref_type]][[ref_name]])))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least two cell types in the signature matrix
    filter(n_shared_celltypes > 1) %>%
    # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
    mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
    # Split dependent cell types in each CIBERSORTx run
    mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
    unnest(celltype_classes) %>%
    # Run BayesPrism
    rowwise() %>%
    mutate(bp_res = list(runBayesPrism(mix = mixture, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes)))



  # Calculate correlations
  getCors <- function(res_mat, truth_mat){
    samples <- intersect(colnames(res_mat), colnames(truth_mat))
    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))
    suppressWarnings(sapply(celltypes, function(ct){
      cor(res_mat[ct,samples], truth_mat[ct,samples], method = "spearman")
    }))
  }

  vals.refs.tbl %>%
    mutate(method = "BayesPrism", .before = everything()) %>%
    mutate(cor = list(getCors(bp_res, truth))) %>%
    unnest_longer(cor, indices_to = "celltype") %>%
    select(method:n_val_samples, n_shared_celltypes, celltype, cor) %>%
    return()

}






# -------------- Cytometry validations --------------

print("Running Cytometry Validations...")

cytoVals.list <- list(blood = c("BG_blood", "GSE107011", "GSE107572", "GSE127813"),
                  tumor = c("ccRCC_cytof_CD45+", "NSCLC_cytof"),
                  other = c("GSE120444"))

cytoVals <- loadVals(cytoVals.list, valMixDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/", valTruthDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/")



# xCell2 -----

print("Running xCell2...")


xcell2.cyto.val.res <- lapply(cyto_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]

    lapply(xcell2_sigs_paths, function(sig_path){
      print(sig_path)
      xcell2sigs <- readRDS(sig_path)

      # Check number of shared cell types
      sigs_celltypes <- unique(gsub("#.*", "", names(xcell2sigs@signatures)))
      shared_celltypes <- intersect(sigs_celltypes, rownames(truth))

      xcell2.out.mat <- xCell2Analysis(mix, xcell2sigs)
      samples <- intersect(colnames(truth), colnames(xcell2.out.mat))
      sapply(shared_celltypes, function(x){
        cor(xcell2.out.mat[x, samples], truth[x, samples], method = "spearman")
      })

    })


  })
})


xCell2.allCors.cytoVal.tidy <- enframe(xcell2.cyto.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "xCell2") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(xCell2.allCors.cytoVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/xCell2.allCors.cytoVal.tidy")
# xCell2.allCors.cytoVal.tidy <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/xCell2.allCors.cytoVal.tidy")


# CIBERSORTx  -----
runCIBERSORTx <- function(mix, sigMatPath, celltypes2use, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  # Subset cell types from the signature matrix
  sigmat_tmp <- read.csv(sigMatPath, sep = "\t", header = T, check.names = F, row.names = 1)
  sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
  sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
  write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Make mixture file
  mix_tmp <- cbind("genes" = rownames(mix), mix)
  mix_file <- paste0(dir, "/mix-tmp.txt")
  write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Run CIBERSORTx
  token <- "b72da36961922443b75a1b65beef27c0"

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }


  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell FALSE --QN FALSE --rmbatchBmode TRUE --verbose TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)


  # Load results
  res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
  cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}
splitDepCellTypes <- function(refName, truth_celltypes, dir = "/bigdata/almogangel/xCell2_data/dev_data/"){
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

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", refName, "_ref.rds"))

  dep_list <- getDependencies(ref$lineage_file)
  dep_list <- dep_list[names(dep_list) %in% truth_celltypes]

  dep_cts <- sapply(names(dep_list), function(d){
    deps <- unname(unlist(dep_list[[d]]))
    deps[deps %in% truth_celltypes]
  })
  dep_cts <- dep_cts[order(sapply(dep_cts, length), decreasing = T)]

  classes <- list("1" = c())
  for (ct in names(dep_cts)) {

    classes <- classes[order(sapply(classes, length), decreasing = F)]
    deps <- dep_cts[[ct]]

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

  return(classes)
}

print("Running CIBERSORTx...")


cbrx.cyto.val.res <- lapply(cyto_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]




    all_sigmats <- lapply(names(cbrx_sigmats_paths), function(ref){

      print(ref)
      truth_celltypes <- rownames(truth)
      sigmat_celltypes_classes <- splitDepCellTypes(ref, truth_celltypes)

      sigmat_path <- cbrx_sigmats_paths[[ref]]

      # Check number of shared cell types
      sigmat_celltypes <- unname(unlist(sigmat_celltypes_classes))
      shared_celltypes <- intersect(sigmat_celltypes, rownames(truth))

      all.classes.out <- lapply(sigmat_celltypes_classes, function(class){
        runCIBERSORTx(mix, sigmat_path, class)
        })
      all.classes.out <- do.call(rbind, all.classes.out)

      samples <- intersect(colnames(truth), colnames(all.classes.out))
      sapply(shared_celltypes, function(x){
        cor(all.classes.out[x, samples], truth[x, samples], method = "spearman")
      })

    })
    names(all_sigmats) <- names(cbrx_sigmats_paths)
    all_sigmats

   })


})

CBRx.allCors.cytoVal.tidy <- enframe(cbrx.cyto.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(CBRx.allCors.cytoVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/CBRx.allCors.cytoVal.tidy")
# CBRx.allCors.cytoVal.tidy <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/CBRx.allCors.cytoVal.tidy")


# Run other methods  -----



# Combine and plot all correlation results -----------------

allCors.tidy <- rbind(CBRx.allCors.tidy, xCell2.allCors.tidy)


allCors.tidy %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=validation, y=cor)) +
  geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  #geom_jitter(aes(col=validation), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  #scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
   #                           "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  #scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(x = "Cytometry validation dataset", title = NULL, y = NULL, colour = NULL, fill = NULL) +
  facet_wrap(~reference, scales = "free_x")


xCell2.allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_wrap(~dataset, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Sorted Cells Validation - Blood (xCell2)", x = NULL, colour = NULL, fill = NULL)


allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=method), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                              "#424242", "#00F5FF", "#FF3E96", "#54FF9F")) +
  facet_wrap(~dataset, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Sorted Cells Validation - Blood", x = NULL, colour = NULL, fill = NULL)









# -------------- Single-cell validations --------------

print("Running Single-cell Validations...")


val_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"

sc_vals <- list(blood = c("sc_pbmc - blood" = "sc_pbmc"),
                  tumor = c("SC_glioblastoma - tumor" = "SC_glioblastoma", "SC_GSE84133 - tumor" = "SC_GSE84133", "SC_HNSCC - tumor" = "SC_HNSCC", "SC_lymphomas - tumor" = "SC_lymphomas",
                            "SC_melanoma - tumor" = "SC_melanoma", "SC_NSCLC - tumor" = "SC_NSCLC"))


# xCell2 -----

print("Running xCell2...")


xcell2.sc.val.res <- lapply(sc_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]

    lapply(xcell2_sigs_paths, function(sig_path){
      print(sig_path)
      xcell2sigs <- readRDS(sig_path)

      # Check number of shared cell types
      sigs_celltypes <- unique(gsub("#.*", "", names(xcell2sigs@signatures)))
      shared_celltypes <- intersect(sigs_celltypes, rownames(truth))

      xcell2.out.mat <- xCell2Analysis(mix, xcell2sigs)
      samples <- intersect(colnames(truth), colnames(xcell2.out.mat))
      sapply(shared_celltypes, function(x){
        cor(xcell2.out.mat[x, samples], truth[x, samples], method = "spearman")
      })

    })


  })
})


xCell2.allCors.scVal.tidy <- enframe(xcell2.sc.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "xCell2") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(xCell2.allCors.scVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/xCell2.allCors.scVal.tidy")

# CIBERSORTx  -----
runCIBERSORTx <- function(mix, sigMatPath, celltypes2use, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  # Subset cell types from the signature matrix
  sigmat_tmp <- read.csv(sigMatPath, sep = "\t", header = T, check.names = F, row.names = 1)
  sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
  sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
  write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Make mixture file
  mix_tmp <- cbind("genes" = rownames(mix), mix)
  mix_file <- paste0(dir, "/mix-tmp.txt")
  write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Run CIBERSORTx
  token <- "b72da36961922443b75a1b65beef27c0"

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }


  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell FALSE --QN FALSE --rmbatchBmode TRUE --verbose TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)


  # Load results
  res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
  cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}
splitDepCellTypes <- function(refName, truth_celltypes, dir = "/bigdata/almogangel/xCell2_data/dev_data/"){
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

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", refName, "_ref.rds"))

  dep_list <- getDependencies(ref$lineage_file)
  dep_list <- dep_list[names(dep_list) %in% truth_celltypes]

  dep_cts <- sapply(names(dep_list), function(d){
    deps <- unname(unlist(dep_list[[d]]))
    deps[deps %in% truth_celltypes]
  })
  dep_cts <- dep_cts[order(sapply(dep_cts, length), decreasing = T)]

  classes <- list("1" = c())
  for (ct in names(dep_cts)) {

    classes <- classes[order(sapply(classes, length), decreasing = F)]
    deps <- dep_cts[[ct]]

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

  return(classes)
}

print("Running CIBERSORTx...")


cbrx.sc.val.res <- lapply(sc_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]




    all_sigmats <- lapply(names(cbrx_sigmats_paths), function(ref){

      print(ref)
      truth_celltypes <- rownames(truth)
      sigmat_celltypes_classes <- splitDepCellTypes(ref, truth_celltypes)

      sigmat_path <- cbrx_sigmats_paths[[ref]]

      # Check number of shared cell types
      sigmat_celltypes <- unname(unlist(sigmat_celltypes_classes))
      shared_celltypes <- intersect(sigmat_celltypes, rownames(truth))

      all.classes.out <- lapply(sigmat_celltypes_classes, function(class){
        runCIBERSORTx(mix, sigmat_path, class)
      })
      all.classes.out <- do.call(rbind, all.classes.out)

      samples <- intersect(colnames(truth), colnames(all.classes.out))
      sapply(shared_celltypes, function(x){
        cor(all.classes.out[x, samples], truth[x, samples], method = "spearman")
      })

    })
    names(all_sigmats) <- names(cbrx_sigmats_paths)
    all_sigmats

  })


})

CBRx.allCors.scVal.tidy <- enframe(cbrx.sc.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(CBRx.allCors.scVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/CBRx.allCors.scVal.tidy")





# Run other methods  -----



# Combine and plot all correlation results -----------------

allCors.scVal.tidy <- rbind(CBRx.allCors.scVal.tidy, xCell2.allCors.scVal.tidy)


allCors.scVal.tidy %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=validation, y=cor)) +
  geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  #geom_jitter(aes(col=validation), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  #scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
  #                           "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  #scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(x = "Single cell validation dataset", title = NULL, y = NULL, colour = NULL, fill = NULL) +
  facet_wrap(~reference, scales = "free_x")








# -------------- Other validations --------------

print("Running Other Validations...")

val_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"

other_vals <- list(blood = c("Globin_Dep - blood" = "GSE53655", "Blood_Staining - blood" = "GSE60424"),
                tumor = c("Cell_Lines_Mix - tumor" = "GSE121127", "WU_ccRCC_RCCTC - tumor" = "WU_ccRCC_RCCTC"),
                other = c("Nasal_Asthma" = "GSE115823"))


# xCell2 -----

print("Running xCell2...")


xcell2.other.val.res <- lapply(other_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]

    lapply(xcell2_sigs_paths, function(sig_path){
      print(sig_path)
      xcell2sigs <- readRDS(sig_path)

      # Check number of shared cell types
      sigs_celltypes <- unique(gsub("#.*", "", names(xcell2sigs@signatures)))
      shared_celltypes <- intersect(sigs_celltypes, rownames(truth))

      xcell2.out.mat <- xCell2Analysis(mix, xcell2sigs)
      samples <- intersect(colnames(truth), colnames(xcell2.out.mat))
      sapply(shared_celltypes, function(x){
        cor(xcell2.out.mat[x, samples], truth[x, samples], method = "spearman")
      })

    })


  })
})



xCell2.allCors.otherVal.tidy <- enframe(xcell2.other.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "xCell2") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(xCell2.allCors.otherVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/xCell2.allCors.otherVal.tidy")

# CIBERSORTx  -----
runCIBERSORTx <- function(mix, sigMatPath, celltypes2use, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  # Subset cell types from the signature matrix
  sigmat_tmp <- read.csv(sigMatPath, sep = "\t", header = T, check.names = F, row.names = 1)
  sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
  sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
  write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Make mixture file
  mix_tmp <- cbind("genes" = rownames(mix), mix)
  mix_file <- paste0(dir, "/mix-tmp.txt")
  write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


  # Run CIBERSORTx
  token <- "b72da36961922443b75a1b65beef27c0"

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }


  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell FALSE --QN FALSE --rmbatchBmode TRUE --verbose TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)


  # Load results
  res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
  cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}
splitDepCellTypes <- function(refName, truth_celltypes, dir = "/bigdata/almogangel/xCell2_data/dev_data/"){
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

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/dev_data/", refName, "_ref.rds"))

  dep_list <- getDependencies(ref$lineage_file)
  dep_list <- dep_list[names(dep_list) %in% truth_celltypes]

  dep_cts <- sapply(names(dep_list), function(d){
    deps <- unname(unlist(dep_list[[d]]))
    deps[deps %in% truth_celltypes]
  })
  dep_cts <- dep_cts[order(sapply(dep_cts, length), decreasing = T)]

  classes <- list("1" = c())
  for (ct in names(dep_cts)) {

    classes <- classes[order(sapply(classes, length), decreasing = F)]
    deps <- dep_cts[[ct]]

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

  return(classes)
}

print("Running CIBERSORTx...")


cbrx.other.val.res <- lapply(other_vals, function(tissue){
  print(tissue)
  lapply(tissue, function(val){
    print(val)

    mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
    rows <- rownames(truth)
    truth <- apply(truth, 2, as.numeric)
    rownames(truth) <- rows
    truth[is.na(truth)] <- 0
    truth <- truth[!duplicated(rows),]




    all_sigmats <- lapply(names(cbrx_sigmats_paths), function(ref){

      print(ref)
      truth_celltypes <- rownames(truth)
      sigmat_celltypes_classes <- splitDepCellTypes(ref, truth_celltypes)

      sigmat_path <- cbrx_sigmats_paths[[ref]]

      # Check number of shared cell types
      sigmat_celltypes <- unname(unlist(sigmat_celltypes_classes))
      shared_celltypes <- intersect(sigmat_celltypes, rownames(truth))

      all.classes.out <- lapply(sigmat_celltypes_classes, function(class){
        runCIBERSORTx(mix, sigmat_path, class)
      })
      all.classes.out <- do.call(rbind, all.classes.out)

      samples <- intersect(colnames(truth), colnames(all.classes.out))
      sapply(shared_celltypes, function(x){
        cor(all.classes.out[x, samples], truth[x, samples], method = "spearman")
      })

    })
    names(all_sigmats) <- names(cbrx_sigmats_paths)
    all_sigmats

  })


})

CBRx.allCors.otherVal.tidy <- enframe(cbrx.other.val.res, value = "data", name = "tissue") %>%
  unnest_longer(data, indices_to = "validation") %>%
  unnest_longer(data, indices_to = "reference") %>%
  unnest_longer(data, indices_to = "celltype", values_to = "cor") %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, tissue, reference, validation, celltype, cor)

saveRDS(CBRx.allCors.otherVal.tidy, "/bigdata/almogangel/xCell2_data/dev_data/CBRx.allCors.otherVal.tidy")





# Run other methods  -----



# Combine and plot all correlation results -----------------

allCors.otherVal.tidy <- rbind(CBRx.allCors.otherVal.tidy, xCell2.allCors.otherVal.tidy)


allCors.otherVal.tidy %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=validation, y=cor)) +
  geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  #geom_jitter(aes(col=validation), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  #scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
  #                           "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  #scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(x = "Other validation dataset", title = NULL, y = NULL, colour = NULL, fill = NULL) +
  facet_wrap(~reference, scales = "free_x")








# -------------- Tabula Sapiens validation --------------
data("ts_labels_with_ontology")

# xCell2 -----
allData <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/")
allData <- allData[allData != "Liver.rds"] # Because contain only two cell types

xcell2.out.list <- lapply(allData, function(data){

  data.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", data))
  ref <- data.in$singleCellExpr
  labels <- tibble(ont = "NA",
                   label = data.in$singleCellLabels,
                   sample = colnames(ref),
                   dataset = data.in$singleCellSubjects)

  # Get ontology
  labels <- labels %>%
    mutate(ont = plyr::mapvalues(sample, ts_labels_with_ontology$sample, ts_labels_with_ontology$ont.fine, warn_missing = FALSE))
  labels <- as.data.frame(labels)

  # Run xCell2.0
  xcell2Sigs <- xCell2Train(ref = ref, labels = labels, data_type = "sc")
  xcell2.out.mat <- xCell2Analysis(bulk = data.in$bulk, xcell2sigs = xcell2Sigs)
  rownames(xcell2.out.mat) <- gsub("-", "_", rownames(xcell2.out.mat))


  # Calculate correlation
  truth <- data.in$bulkRatio
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(xcell2.out.mat[x,], truth[x,], method = "spearman")
  })
})
names(xcell2.out.list) <- allData

xCell2.allCors.tidy <- enframe(xcell2.out.list, value = "cor", name = "tissue") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0) %>%
  mutate(method = "xCell2") %>%
  select(method, tissue, cor, celltype)

saveRDS(xCell2.allCors.tidy, "/bigdata/almogangel/twelve_years_decon_paper/analysis/xcell2_correlations_080523.rds")

# Load this
xCell2.allCors.tidy <- readRDS("/bigdata/almogangel/twelve_years_decon_paper/analysis/xcell2_correlations_080523.rds")



# Load single-cell validation results for all other methods -----
allRes <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/results/accuracy/")
allRes <- allRes[allRes != "Liver.rds"] # Because contain only two cell types

allCors <- sapply(allRes, function(res){

  # Load results
  res.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/results/accuracy/", res))

  # Calculate correlation
  prop <- t(res.in$P)
  truth <- res.in$groundTruth
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(prop[x,], truth[x,], method = "spearman")
  })
})

allCors.tidy <- enframe(allCors, value = "cor") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  separate(name, into = c("method", "tissue"), sep="_", extra = "merge") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0)

# CIBERSORTx -----

# Make results - only run once
runCIBERSORTx <- function(mix, refsample, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  token <- "b72da36961922443b75a1b65beef27c0"

  mix_tmp <- cbind("genes" = rownames(mix), mix)
  mix_file <- paste0(dir, "/mix-tmp.txt")
  write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

  refsample_tmp <- cbind("genes" = rownames(refsample), refsample)
  refsample_file <- paste0(dir, "/refsample-tmp.txt")
  write.table(refsample_tmp, file = refsample_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

  # Make results directory
  results_dir <- paste0(dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }

  cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --single_cell TRUE --refsample ", refsample_file, " --mixture ", mix_file, " --rmbatchSmode TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  # Load results
  cibersortx_out <- t(read.table(paste0(results_dir,  "/CIBERSORTx_Adjusted.txt"), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}

allData <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/")
allData <- allData[allData != "Liver.rds"] # Because contain only two cell types

cbrx.out.list <- lapply(allData, function(data){

  data.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", data))
  mix <- data.in$bulk
  singleCellExpr <- data.in$singleCellExpr
  colnames(singleCellExpr) <- data.in$singleCellLabels
  refsample <- singleCellExpr

  cbrx.out.mat <- runCIBERSORTx(mix, refsample, dir = "/bigdata/almogangel/CIBERSORTx_docker")

  # Calculate correlation
  truth <- data.in$bulkRatio
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(cbrx.out.mat[x,], truth[x,], method = "spearman")
  })
})

cbrx.out.list <- all_validation_combined_results$cbrx.out.list
names(cbrx.out.list) <- allData

CBRx.allCors.tidy <- enframe(cbrx.out.list, value = "cor", name = "tissue") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0) %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, tissue, cor, celltype)

saveRDS(CBRx.allCors.tidy, "/bigdata/almogangel/twelve_years_decon_paper/analysis/cibersortx_correlations.rds")

# Load this
CBRx.allCors.tidy <- readRDS("/bigdata/almogangel/twelve_years_decon_paper/analysis/cibersortx_correlations.rds")

CBRx.allCors.tidy <- CBRx.allCors.tidy %>%
  filter(tissue != "Liver")

# Combine and plot all correlation results -----------------

allCors.tidy <- rbind(allCors.tidy, CBRx.allCors.tidy, xCell2.allCors.tidy)

allCors.tidy <- allCors.tidy %>%
  filter(tissue != "Liver")

allCors.tidy.median <- allCors.tidy %>%
  group_by(method, tissue) %>%
  summarise(cor = median(cor))

method_sorted <- allCors.tidy.median %>%
  group_by(method) %>%
  summarise(cor = median(cor)) %>%
  arrange(cor) %>%
  pull(method)

allCors.tidy.median$method <- factor(allCors.tidy.median$method, levels = method_sorted)

allCors.tidy.median %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=tissue), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median of all cell types)", title = "Tabula Sapiens Validation", x = NULL, colour = NULL, fill = NULL)


xCell2.allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=tissue), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_wrap(~tissue, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Tabula Sapiens Validation - xCell2", x = NULL, colour = NULL, fill = NULL)


rbind(allCors.tidy, CBRx.allCors.tidy) %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=method), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                              "#424242", "#00F5FF", "#FF3E96", "#54FF9F")) +
  facet_wrap(~tissue, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Tabula Sapiens Validation", x = NULL, colour = NULL, fill = NULL)
