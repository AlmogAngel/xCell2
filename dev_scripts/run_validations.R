library(tidyverse)
library(EPIC)
library(BayesPrism)
library(MCPcounter)
library(dtangle)
library(DeconRNASeq)

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data")

getDependencies <- function(lineage_file_checked){
  ont <- read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
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


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
# Read reference-validation pairs
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
# refval.tbl.nodeps <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val_nodeps.rds")
sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# sc.refval.tbl.nodeps <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val_nodeps.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Method functions

# Run CIBERSORTx
getCIBERSORTxRes <- function(ref_val_table, vals, celltype_conversion){


  runCIBERSORTx <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "/bigdata/almogangel/CIBERSORTx_docker", suffix = "7jan"){


    getResults <- function(cell_types, dir, refName, valName, valType, single_cell, useQN){

      # Subset cell types from the signature matrix
      if (refName == "lm22") {
        # For CIBERSORTx use non RMA normalized signature matrix as it preform quantile normalization
        sigmat_tmp <- read.csv("/bigdata/almogangel/CIBERSORTx_docker/lm22_sigmat.txt", sep = "\t", header = T, check.names = F, row.names = 1)
      }else{
        sigmat_tmp <- read.csv(paste0(dir, "/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      }

      sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,cell_types])
      sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
      write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

      if (single_cell) {
        # Subset reference file
        ref_tmp <- read.csv(paste0(dir, "/", refName, "_ref.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
        index <- colnames(ref_tmp) %in% cell_types
        names <- colnames(ref_tmp)[index]
        ref_tmp <- cbind(rownames(ref_tmp), ref_tmp[,index])
        colnames(ref_tmp) <- c("genes", names)
        ref_tmp_file <- paste0(dir, "/ref-tmp.txt")
        write.table(ref_tmp, file = ref_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
      }

      # Make mixture file
      mix <- vals$mixtures[[valType]][[valName]]

      # CIBERSORTx work on non-log space
      if(max(mix) < 50){
        mix <- (2^mix)+1
      }

      mix_tmp <- cbind("genes" = rownames(mix), mix)
      mix_file <- paste0(dir, "/mix-tmp.txt")
      write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


      # Run CIBERSORTx
      token <- readLines("/bigdata/almogangel/xCell2_data/CIBERSORTx_token.txt", n = 1)

      # Make results directory
      results_dir <- paste0(dir, "/results")
      if (!dir.exists(results_dir)){
        dir.create(results_dir)
      }

      # Clean old results
      if(length(list.files(results_dir)) > 0){
        system(paste0("rm -f ", results_dir, "/*"))
      }


      if (single_cell) {
        cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                      token,  " --mixture ", mix_file, " --single_cell ", single_cell , " --rmbatchSmode ", single_cell, " --QN ", useQN,
                      " --refsample ", ref_tmp_file, " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
      }else{
        cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                      token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell ," --rmbatchBmode ", !single_cell, " --QN ", useQN,
                      " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
      }



      # Run Docker via shell
      system(cmd, wait = TRUE)


      # Load results
      res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
      cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
      cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]

      return(cibersortx_out)
    }


    print(paste0(valName, "_", refName))
    single_cell <- ifelse(refType == "sc", TRUE, FALSE)
    useQN <- ifelse(refType == "array", TRUE, FALSE)


    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/CIBERSORTx#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      cibersortx_out <- readRDS(results_file)
      return(cibersortx_out)
    }

    # Get dependencies
    if (refName == "lm22") {
      # For CIBERSORTx use non RMA normalized reference as it preform quantile normalization
      ref.in <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")
    }else{
      ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
    }
    lineage_file <- ref.in$lineage_file
    rm(ref.in)
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)

    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))


    cibersortx_out <- getResults(leaf_cell_types, dir, refName, valName, valType, single_cell, useQN)


    if (all(celltypes2use %in% rownames(cibersortx_out))) {

      cibersortx_out <- cibersortx_out[celltypes2use,]
      saveRDS(cibersortx_out, results_file)
      return(cibersortx_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(cibersortx_out)]

      for (ancestor in cell_types_left) {

        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, dir, refName, valName, valType, single_cell, useQN)
        cibersortx_out <- rbind(cibersortx_out,  ancestor_out[ancestor,])
        rownames(cibersortx_out)[nrow(cibersortx_out)] <- ancestor

      }

      saveRDS(cibersortx_out, results_file)
      return(cibersortx_out)
    }

  }


  # i = 51
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]


# Run CIBERSORTx
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runCIBERSORTx(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "CIBERSORTx", .before = everything())


  return(vals.refs.res)

}

# Run EPIC
getEPICRes <- function(ref_val_table, vals, celltype_conversion){

  runEPIC <- function(vals, default_refs, valType, valName, refsRDSList, refName, refType, refTissue, celltypes2use, dir = "references", suffix = "7jan"){


    getResults <- function(cell_types, ref.in, refName, refType, dir){

      # Subset sigGenes from the signature matrix
      sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      sigGenes <- rownames(sigmat)

      # Subset GEP
      gep <- read.csv(paste0(dir, "/gep/", refName, "_gep.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      gep <- gep[,cell_types]

      # Subset cell types from the raw reference matrix
      celltypeIndex <- ref.in$labels$label %in% cell_types
      ref.raw <- as.matrix(ref.in$ref)[,celltypeIndex]
      colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label

      # EPIC does not recommend log space
      if(max(gep) < 50){
        gep <- (2^gep)+1
      }

      if(max(ref.raw) < 50){
        ref.raw <- (2^ref.raw)+1
      }


      # Generate reference for EPIC
      ref.raw.var <- sapply(unique(colnames(ref.raw)), function(ct){
        if(sum(colnames(ref.raw) == ct) > 1){
          apply(ref.raw[,colnames(ref.raw) == ct], 1, sd)
        }else{
          rep(0, length(ref.raw[,colnames(ref.raw) == ct]))
        }
      })
      rownames(ref.raw.var) <- rownames(ref.raw)
      ref.raw.var <- ref.raw.var[rownames(gep), colnames(gep)]


      epic_ref <- list("sigGenes" = sigGenes,
                       "refProfiles" = as.matrix(gep),
                       "refProfiles.var" = ref.raw.var)

      # Run EPIC
      res <- t(EPIC(bulk=mix, reference=epic_ref, withOtherCells=FALSE)$cellFractions)

      return(res)
    }

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # EPIC does not recommend log space
    if(max(mix) < 50){
      mix <- (2^mix)+1
    }


    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/EPIC#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      epic_out <- readRDS(results_file)
      return(epic_out)
    }


    ref.in <- refsRDSList[[refType]][[refName]]

    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)


    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))


    epic_out <- getResults(leaf_cell_types, ref.in, refName, refType, dir)


    if (all(celltypes2use %in% rownames(epic_out))) {


      epic_out <- epic_out[celltypes2use,]
      saveRDS(epic_out, results_file)
      return(epic_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(epic_out)]

      for (ancestor in cell_types_left) {

        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, ref.in, refName, refType, dir)
        epic_out <- rbind(epic_out,  ancestor_out[ancestor,])
        rownames(epic_out)[nrow(epic_out)] <- ancestor

      }

      epic_out <- epic_out[celltypes2use,]
      saveRDS(epic_out, results_file)
      return(epic_out)
    }

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


  # i = 51
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]


  # Run EPIC
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runEPIC(vals, default_refs, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, refTissue = ref_tissue, celltypes2use = celltype_classes))) %>%
    mutate(method = "EPIC", .before = everything())


  return(vals.refs.res)


}

# Run BayesPrism
getBayesPrismRes <- function(ref_val_table, vals, celltype_conversion){

  runBayesPrism <- function(vals, valType, valName, refsRDSList, refName, refType, celltypes2use, CPUs = 45, suffix = "7jan"){

    dir.create(tempdir())

    getResults <- function(cell_types, ref.in, refType, CPUs){

      celltypeIndex <- ref.in$labels$label %in% cell_types
      ref.raw <- as.matrix(t(ref.in$ref[,celltypeIndex]))

      # BayesPrism does not recommend log-space data
      if(max(ref.raw) < 50){
        ref.raw <- (2^ref.raw)+1
      }

      type <- ifelse(refType == "sc", "count.matrix", "GEP")
      ref.filtered <- cleanup.genes(input=ref.raw,
                                    input.type=type,
                                    species="hs",
                                    gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"))


      ref.filtered.pc <- select.gene.type(ref.filtered,
                                          gene.type = "protein_coding")

      labels <- ref.in$labels[celltypeIndex,]$label

      if (refType == "sc") {
        diff.exp.stat <- get.exp.stat(sc.dat=ref.raw[,colSums(ref.raw>0)>3],# filter genes to reduce memory use
                                      cell.type.labels=labels,
                                      cell.state.labels=labels,
                                      # psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                      cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                      n.cores=CPUs) #number of threads


        ref.filtered.pc <- select.marker(sc.dat=ref.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)
      }else{

        # Subset marker genes from the signature matrix
        markers <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
        shared_markers <- intersect(unique(markers$marker), colnames(ref.filtered.pc))
        ref.filtered.pc <- ref.filtered.pc[,shared_markers]
      }


      tumorKey <- NULL
      if ("malignant cell" %in% labels) {
        tumorKey <- "malignant cell"
      }

      myPrism <- new.prism(
        reference=ref.filtered.pc,
        mixture=t(mix),
        input.type=type,
        cell.type.labels = labels,
        cell.state.labels =labels,
        key=tumorKey)

      bp.res <- run.prism(prism = myPrism, n.cores = CPUs)


      res <- t(get.fraction (bp=bp.res,
                             which.theta="final",
                             state.or.type="type"))

      return(res)

    }


    print(paste0(valName, "_", refName))

    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/BayesPrism#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      bayes_out <- readRDS(results_file)
      return(bayes_out)
    }


    mix <- vals$mixtures[[valType]][[valName]]
    ref.in <- refsRDSList[[refType]][[refName]]

    # BayesPrism does not recommend log-space data
    if(max(mix) < 50){
      mix <- (2^mix)+1
    }


    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)
    all_celltypes <- names(dep_list)

    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))

    bayes_out <- getResults(leaf_cell_types, ref.in, refType, CPUs)


    if (all(celltypes2use %in% rownames(bayes_out))) {


      bayes_out <- bayes_out[celltypes2use,]
      saveRDS(bayes_out, results_file)
      return(bayes_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(bayes_out)]

      for (ancestor in cell_types_left) {

        dir.create(tempdir())
        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, ref.in, refType, CPUs)
        bayes_out <- rbind(bayes_out,  ancestor_out[ancestor,])
        rownames(bayes_out)[nrow(bayes_out)] <- ancestor

      }

      bayes_out <- bayes_out[celltypes2use,]
      saveRDS(bayes_out, results_file)
      return(bayes_out)
    }


  }

  # Load references matrices
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  refsRDSList <- lapply(refList, function(ref_type){
    refs <- lapply(ref_type, function(ref){
      # Load reference
      ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
      ref.in
    })
    names(refs) <- ref_type
    refs
  })

  # i = 51
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]
  #

  vals.refs.res <- ref_val_table %>%
    # Run BayesPrism
    rowwise() %>%
    mutate(res = list(runBayesPrism(vals, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "BayesPrism", .before = everything())



  return(vals.refs.res)

}

# Run MCPcounter
getMCPcounterRes <- function(ref_val_table, vals, celltype_conversion){

  runMCPcounter <- function(vals, valType, valName, markers, celltypes2use, refName){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # It is recommended to use log space with MCPcounter
    if(max(mix) >= 50){
      mix <-log2(mix+1)
    }

    markers_tmp <- markers[celltypes2use]

    markers_tmp <- enframe(markers_tmp, value = "HUGO symbols", name = "Cell population") %>%
      unnest(`HUGO symbols`) %>%
      dplyr::select(`HUGO symbols`, `Cell population`) %>%
      as.data.frame()

    res <- MCPcounter.estimate(expression = mix, featuresType = "HUGO_symbols", genes = markers_tmp)

    return(res)
  }

  # Get marker genes
  refs <- unique(ref_val_table$ref_name)
  refs_markers <- sapply(refs, function(refName){
    markers <- read.csv(paste0("references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
    split(markers$marker, markers$label)
  })


  vals.refs.res <- ref_val_table %>%
    # Add marker genes
    mutate(markers = list(refs_markers[[ref_name]])) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(markers)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    filter(n_shared_celltypes > 2) %>%
    # Run MCPcounter
    mutate(res = list(runMCPcounter(vals, valType = val_type, valName = val_dataset, markers, celltypes2use = shared_celltypes, refName = ref_name))) %>%
    mutate(method = "MCPcounter", .before = everything())


  return(vals.refs.res)

}

# Run dtangle
getdtangleRes <- function(ref_val_table, vals, celltype_conversion){

  rundtangle <- function(vals, valType, celltypes2use, refsRDSList, valName, valDataType, refName, refType, suffix = "7jan"){

    getResults <- function(cell_types, ref.in, refName, refType){


      # Subset cell types from reference
      celltypeIndex <- ref.in$labels$label %in% cell_types
      ref.raw <- as.matrix(ref.in$ref)[,celltypeIndex]
      colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label

      # Convert scRNA-Seq reference expression to log-CPM
      if (refType == "sc") {
        lib_sizes <- Matrix::colSums(ref.raw)
        norm_factor <- 1e6 / lib_sizes
        ref_norm <- ref.raw %*% Matrix::Diagonal(x = norm_factor)
        colnames(ref_norm) <- colnames(ref.raw)
        ref_norm <- log2(ref_norm+1)
        ref.raw <- ref_norm
      }

      shared_genes <- intersect(rownames(mix), rownames(ref.raw))
      mix <- mix[shared_genes,]
      ref.raw <- ref.raw[shared_genes,]
      ref.raw <- as.matrix(ref.raw)

      # Make sure reference in log space
      if (max(ref.raw) >= 50) {
        ref.raw <-log2(ref.raw+1)
      }


      if (refType == "sc") {
        y <- cbind(ref.raw, mix)
        y <- limma::normalizeBetweenArrays(y)
        ref.raw <- y[,1:ncol(ref.raw)]
        mix <- y[,(ncol(ref.raw)+1):ncol(y)]
      }


      celltypes <- unique(colnames(ref.raw))
      pure_samples_list <- lapply(celltypes, function(ct){
        which(colnames(ref.raw) == ct)
      })
      names(pure_samples_list) <- celltypes

      markerMethod <- ifelse(refType == "sc", "ratio", "p.value")
      dataType <- ifelse(refType == "array", "microarray-gene", "rna-seq")
      res <- dtangle(Y = t(mix), references = t(ref.raw), pure_samples = pure_samples_list, marker_method = markerMethod, data_type = dataType)
      res <- t(res$estimates)

      return(res)

    }


    print(paste0(valName, "_", refName))

    # Load mixture
    mix <- vals$mixtures[[valType]][[valName]]

    # dtangle is recommended in log space
    if(max(mix) >= 50){
      mix <-log2(mix+1)
    }

    # Check if results already exist
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/dtangle#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      dtangle_out <- readRDS(results_file)
      return(dtangle_out)
    }

    # Load reference data
    ref.in <- refsRDSList[[refType]][[refName]]

    # Get cell type dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)


    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))

    # Run dtangle for leaf cell types
    dtangle_out <- getResults(leaf_cell_types, ref.in, refName, refType)

    if (all(celltypes2use %in% rownames(dtangle_out))) {


      dtangle_out <- dtangle_out[celltypes2use,]
      saveRDS(dtangle_out, results_file)
      return(dtangle_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(dtangle_out)]

      for (ancestor in cell_types_left) {

        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, ref.in, refName, refType)
        dtangle_out <- rbind(dtangle_out,  ancestor_out[ancestor,])
        rownames(dtangle_out)[nrow(dtangle_out)] <- ancestor

      }

      dtangle_out <- dtangle_out[celltypes2use,]
      saveRDS(dtangle_out, results_file)
      return(dtangle_out)
    }

  }

  # Load references
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  refsRDSList <- lapply(refList, function(ref_type){
    refs <- lapply(ref_type, function(ref){
      ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
      ref.in
    })
    names(refs) <- ref_type
    refs
  })

  # i = 51
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]


  vals.refs.res <- ref_val_table %>%
    # Run dtangle
    rowwise() %>%
    mutate(res = list(rundtangle(vals, valType = val_type, valName = val_dataset, valDataType = val_data_type, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "dtangle", .before = everything())




  return(vals.refs.res)

}

# Run DeconRNASeq
getDeconRNASeqRes <- function(ref_val_table, vals, celltype_conversion){

  runDeconRNASeq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references", suffix = "7jan"){


    getResults <- function(cell_types, refName){

      sigmat <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      sigmat <- sigmat[,cell_types]

      if(max(sigmat) < 50){
        sigmat <- (2^sigmat)+1
      }

      # Run DeconRNASeq
      res <- DeconRNASeq(data.frame(mix), sigmat, checksig=FALSE, known.prop = FALSE,  fig = FALSE)
      res <- as.matrix(res$out.all)
      rownames(res) <- colnames(mix)
      res <- t(res)

      return(res)
    }


    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # DeconRNASeq recommended not in log space
    if(max(mix) < 50){
      mix <- (2^mix)+1
    }


    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/DeconRNASeq#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      decon_out <- readRDS(results_file)
      return(decon_out)
    }


    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]

    # if (refName == "lm22") {
    #   ref.in <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref_for_dtangle.rds")
    # }


    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)


    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))

    decon_out <- getResults(leaf_cell_types, refName)

    if (all(celltypes2use %in% rownames(decon_out))) {


      decon_out <- decon_out[celltypes2use,]
      saveRDS(decon_out, results_file)
      return(decon_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(decon_out)]

      for (ancestor in cell_types_left) {

        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, refName)
        decon_out <- rbind(decon_out,  ancestor_out[ancestor,])
        rownames(decon_out)[nrow(decon_out)] <- ancestor

      }

      decon_out <- decon_out[celltypes2use,]
      saveRDS(decon_out, results_file)
      return(decon_out)
    }

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

  # i = 4
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]


  vals.refs.res <- ref_val_table %>%
    # Run DeconRNASeq
    rowwise() %>%
    mutate(res = list(runDeconRNASeq(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "DeconRNASeq", .before = everything())


  return(vals.refs.res)

}

# Run quanTIseq
# Workaround from: https://github.com/icbi-lab/quanTIseq/tree/master/quantiseq/deconvolution
# as used in "Twelve Years of Cellular Deconvolution: Applications, Benchmark, Methodology, and Challenges"
# TODO: check this workaround
getquanTIseqRes <- function(ref_val_table, vals, celltype_conversion){

  source("/bigdata/almogangel/xCell2/dev_scripts/quantiseq_code.R")
  runquanTIseq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references", suffix = "7jan"){

    getResults <- function(cell_types, refName, refType, dir){



      # Subset sigGenes from the signature matrix
      sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      sigmat <- sigmat[,cell_types]


      if(max(sigmat) < 50){
        sigmat <- (2^sigmat)+1
      }

      # Run quanTIseq
      res <- t(quanTIseq(sigmat, mix, scaling=rep(1, ncol(sigmat)), method="lsei"))
      res <- res[rownames(res) != "Other",]

      return(res)
    }


    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    if(max(mix) < 50){
      mix <- (2^mix)+1
    }

    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/other_methods/quanTIseq#", valName, "#", refName, "#", suffix, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      quanti_out <- readRDS(results_file)
      return(quanti_out)
    }

    ref.in <- refsRDSList[[refType]][[refName]]

    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)


    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))



    quanti_out <- getResults(leaf_cell_types, refName, refType, dir)


    if (all(celltypes2use %in% rownames(quanti_out))) {


      quanti_out <- quanti_out[celltypes2use,]
      saveRDS(quanti_out, results_file)
      return(quanti_out)

    }else{

      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(quanti_out)]

      for (ancestor in cell_types_left) {

        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(ct2use, refName, refType, dir)
        quanti_out <- rbind(quanti_out,  ancestor_out[ancestor,])
        rownames(quanti_out)[nrow(quanti_out)] <- ancestor

      }

      quanti_out <- quanti_out[celltypes2use,]
      saveRDS(quanti_out, results_file)
      return(quanti_out)
    }



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


  # i = 4
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]


  vals.refs.res <- ref_val_table %>%
    # Run quanTIseq
    rowwise() %>%
    mutate(res = list(runquanTIseq(vals, valType =  val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "quanTIseq", .before = everything())


  return(vals.refs.res)

}


# Run cytometry/other validations -----------------------

print("Running CIBERSORTx...")
cbrx.cyto.res <- getCIBERSORTxRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(cbrx.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.cyto.res.rds")

print("Running EPIC...")
epic.cyto.res <- getEPICRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(epic.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.cyto.res.rds")

print("Running BayesPrism...")
bp.cyto.res <- getBayesPrismRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(bp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/bp.cyto.res.rds")

print("Running MCPcounter")
mcp.cyto.res <- getMCPcounterRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(mcp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mcp.cyto.res.rds")

print("Running dtangle")
dtan.cyto.res <- getdtangleRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(dtan.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.cyto.res.rds")

print("Running DeconRNASeq")
decon.cyto.res <- getDeconRNASeqRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(decon.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/decon.cyto.res.rds")

print("Running quanTIseq")
quanti.cyto.res <- getquanTIseqRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(quanti.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/quanti.cyto.res.rds")



# Run single-cell validations -----------------------

print("Running CIBERSORTx...")
#cbrx.sc.res <- getCIBERSORTxRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
#saveRDS(cbrx.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.sc.res.rds")

print("Running EPIC...")
#epic.sc.res <- getEPICRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
#saveRDS(epic.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.sc.res.rds")

print("Running BayesPrism...")
bp.sc.res <- getBayesPrismRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(bp.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/bp.sc.res.rds")

print("Running MCPcounter")
mcp.sc.res <- getMCPcounterRes(ref_val_table = sc.refval.tbl, vals = sc.vals, celltype_conversion)
saveRDS(mcp.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mcp.sc.res.rds")

print("Running dtangle")
dtan.sc.res <- getdtangleRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(dtan.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.sc.res.rds")

print("Running DeconRNASeq")
decon.sc.res <- getDeconRNASeqRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(decon.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/decon.sc.res.rds")

print("Running quanTIseq")
quanti.sc.res <- getquanTIseqRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(quanti.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/quanti.sc.res.rds")

