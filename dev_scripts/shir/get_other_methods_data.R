library(tidyverse)
library(EPIC)
library(BayesPrism)
library(MCPcounter)
library(dtangle)
library(DeconRNASeq)
library(omnideconv)


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



runCIBERSORTx <- function(mix, mixType, mixName, refName, refType, dir = "/bigdata/almogangel/CIBERSORTx_docker", suffix = "9jul"){


  getResults <- function(mix, cell_types, dir, refName, mixName, mixType, single_cell, useQN){

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


    # CIBERSORTx work on non-log space
    if(max(mix) < 50){
      mix <- (2^mix)-1
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


  print(paste0(mixName, "_", refName))
  single_cell <- ifelse(refType == "sc", TRUE, FALSE)
  useQN <- ifelse(refType == "array", TRUE, FALSE)


  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/CIBERSORTx#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
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
  if (!file.exists(lineage_file)) {
    stop("lineage file missing!")
  }
  dep_list <- getDependencies(lineage_file)

  # Get leaf cell types (no descendants)
  leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))

  celltypes2use <- unique(ref.in$labels$label)
  rm(ref.in)
  cibersortx_out <- getResults(mix, leaf_cell_types, dir, refName, mixName, mixType, single_cell, useQN)


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
      ancestor_out <- getResults(mix, ct2use, dir, refName, mixName, mixType, single_cell, useQN)
      cibersortx_out <- rbind(cibersortx_out,  ancestor_out[ancestor,])
      rownames(cibersortx_out)[nrow(cibersortx_out)] <- ancestor

    }

    saveRDS(cibersortx_out, results_file)
    return(cibersortx_out)
  }

}

runEPIC <- function(mix, mixType, mixName, refName, refType, dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul"){


  getResults <- function(mix, cell_types, ref.in, refName, refType, dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references"){

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
      gep <- (2^gep)-1
    }

    if(max(ref.raw) < 50){
      ref.raw <- (2^ref.raw)-1
    }


    genes2use <- intersect(rownames(mix), rownames(gep))

    # Generate reference for EPIC
    ref.raw <- ref.raw[genes2use, ]
    ref.raw.var <- sapply(unique(colnames(ref.raw)), function(ct){
      if(sum(colnames(ref.raw) == ct) > 1){
        apply(ref.raw[,colnames(ref.raw) == ct], 1, sd)
      }else{
        rep(0, length(ref.raw[,colnames(ref.raw) == ct]))
      }
    })
    rownames(ref.raw.var) <- rownames(ref.raw)
    ref.raw.var <- ref.raw.var[genes2use, colnames(gep)]


    sigGenes <- sigGenes[sigGenes %in% genes2use]

    epic_ref <- list("refProfiles" = as.matrix(gep)[genes2use,],
                     "sigGenes" = sigGenes,
                     "refProfiles.var" = ref.raw.var)

    # Run EPIC
    res <- t(EPIC(bulk=mix[genes2use,], reference=epic_ref, withOtherCells=FALSE)$cellFractions)

    return(res)
  }

  print(paste0(mixName, "_", refName))

  # EPIC does not recommend log space
  if(max(mix) < 50){
    mix <- (2^mix)-1
  }


  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/EPIC#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
    epic_out <- readRDS(results_file)
    return(epic_out)
  }


  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

  # Get dependencies
  lineage_file <- ref.in$lineage_file
  if (!file.exists(lineage_file)) {
    stop("lineage file missing!")
  }
  dep_list <- getDependencies(lineage_file)


  # Get leaf cell types (no descendants)
  leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))


  epic_out <- getResults(mix, leaf_cell_types, ref.in, refName, refType, dir)
  celltypes2use <- unique(ref.in$labels$label)


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
      ancestor_out <- getResults(mix, ct2use, ref.in, refName, refType, dir)
      epic_out <- rbind(epic_out,  ancestor_out[ancestor,])
      rownames(epic_out)[nrow(epic_out)] <- ancestor

    }

    epic_out <- epic_out[celltypes2use,]
    saveRDS(epic_out, results_file)
    return(epic_out)
  }

}

runBayesPrism <- function(mix, mixType, mixName, refName, refType, CPUs = 45, suffix = "9jul"){

  dir.create(tempdir())

  getResults <- function(mix, cell_types, ref.in, refType, CPUs){

    celltypeIndex <- ref.in$labels$label %in% cell_types
    ref.raw <- as.matrix(t(ref.in$ref[,celltypeIndex]))

    # BayesPrism does not recommend log-space data
    if(max(ref.raw) < 50){
      ref.raw <- (2^ref.raw)-1
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

  print(paste0(mixName, "_", refName))

  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/BayesPrism#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
    bayes_out <- readRDS(results_file)
    return(bayes_out)
  }


  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

  # BayesPrism does not recommend log-space data
  if(max(mix) < 50){
    mix <- (2^mix)-1
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

  bayes_out <- getResults(mix, leaf_cell_types, ref.in, refType, CPUs)
  celltypes2use <- unique(ref.in$labels$label)


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
      ancestor_out <- getResults(mix, ct2use, ref.in, refType, CPUs)
      bayes_out <- rbind(bayes_out,  ancestor_out[ancestor,])
      rownames(bayes_out)[nrow(bayes_out)] <- ancestor

    }

    bayes_out <- bayes_out[celltypes2use,]
    saveRDS(bayes_out, results_file)
    return(bayes_out)
  }


}

runMCPcounter <- function(mix, mixType, mixName, refName){

  print(paste0(mixName, "_", refName))

  # It is recommended to use log space with MCPcounter
  if(max(mix) >= 50){
    mix <-log2(mix+1)
  }

  # Get marker genes
  markers <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
  markers <- split(markers$marker, markers$label)


  markers_tmp <- markers

  markers_tmp <- enframe(markers_tmp, value = "HUGO symbols", name = "Cell population") %>%
    unnest(`HUGO symbols`) %>%
    dplyr::select(`HUGO symbols`, `Cell population`) %>%
    as.data.frame()

  res <- MCPcounter.estimate(expression = mix, featuresType = "HUGO_symbols", genes = markers_tmp)

  return(res)
}

rundtangle <- function(mix, mixType, mixName, refName, refType, suffix = "9jul"){

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


  print(paste0(mixName, "_", refName))

  # dtangle is recommended in log space
  if(max(mix) >= 50){
    mix <-log2(mix+1)
  }

  # Check if results already exist
  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/dtangle#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
    dtangle_out <- readRDS(results_file)
    return(dtangle_out)
  }

  # Load reference data
  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

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

  celltypes2use <- unique(ref.in$labels$label)

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

runDeconRNASeq <- function(mix, mixType, mixName, refName, refType, dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul"){


  getResults <- function(cell_types, refName){

    sigmat <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,cell_types]

    if(max(sigmat) < 50){
      sigmat <- (2^sigmat)-1
    }

    # Run DeconRNASeq
    res <- DeconRNASeq(data.frame(mix), sigmat, checksig=FALSE, known.prop = FALSE,  fig = FALSE)
    res <- as.matrix(res$out.all)
    rownames(res) <- colnames(mix)
    res <- t(res)

    return(res)
  }


  print(paste0(mixName, "_", refName))

  # DeconRNASeq recommended not in log space
  if(max(mix) < 50){
    mix <- (2^mix)-1
  }


  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/DeconRNASeq#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
    decon_out <- readRDS(results_file)
    return(decon_out)
  }


  # Subset cell types from the raw reference matrix
  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

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
  celltypes2use <- unique(ref.in$labels$label)


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

runBisqueRes <- function(mix, refName){

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

  ref.in <- ref$ref
  labels.in <- ref$labels$label
  ds.in <- ref$labels$dataset

  if (length(unique(ds.in)) == 1) {
    ds.in[sample(1:length(ds.in), round(length(ds.in)/2))] <- paste0(unique(ds.in), "_2")
  }

  bisque.out <- omnideconv::deconvolute_bisque(single_cell_object = as.matrix(ref.in),
                                               bulk_gene_expression = as.matrix(mix),
                                               cell_type_annotations = labels.in,
                                               batch_ids = ds.in
  )
  return(bisque.out$bulk.props)
}

runDWLSRes <- function(mix, refName){

  dwls.sigmat <- read.csv(paste0("/bigdata/almogangel/CIBERSORTx_docker/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)

  genes <- intersect(rownames(dwls.sigmat), rownames(mix))
  mix <- mix[genes, , drop = FALSE]
  dwls.sigmat <- dwls.sigmat[genes, , drop = FALSE]

  solutions_ols <- parallel::mclapply(1:ncol(mix), function(i){
    mix.in_i <- mix[, i]
    DWLS::solveSVR(as.matrix(dwls.sigmat), mix.in_i)

  }, mc.cores = 20)
  names(solutions_ols) <- colnames(mix)
  dwls.out <- do.call(cbind, solutions_ols)

  return(dwls.out)
}

runSCDCRes <- function(mix, refName){

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

  ref.in <- ref$ref
  labels.in <- ref$labels$label
  ds.in <- ref$labels$dataset

  scdc.out <- omnideconv::deconvolute_scdc(single_cell_object = as.matrix(ref.in),
                                           bulk_gene_expression = mix,
                                           cell_type_annotations = labels.in,
                                           batch_ids = ds.in
  )


  return(t(scdc.out$prop.est.mvw))
}

runScadenRes <- function(mix, refName){

  ref <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))

  ref.in <- ref$ref
  labels.in <- ref$labels$label
  ds.in <- ref$labels$dataset


  labels_tmp <- gsub("[^[:alnum:]]+", ".", labels.in)

  val_dataset <- sample(1:10000, 1)
  model_dir <- paste0("/bigdata/almogangel/scaden/icb/", refName, "_", val_dataset)

  if (!dir.exists(model_dir)) {
    dir.create(model_dir)
  }
  scaden.model <- omnideconv::build_model_scaden(single_cell_object = as.matrix(ref.in),
                                                 bulk_gene_expression = mix,
                                                 cell_type_annotations = labels_tmp,
                                                 model_path = model_dir)


  scaden.out <- omnideconv::deconvolute_scaden(
    scaden.model,
    bulk_gene_expression = as.matrix(mix),
    temp_dir = NULL,
    verbose = FALSE
  )


  find_closest_match <- function(source_vec, target_vec) {
    dist_matrix <- stringdist::stringdistmatrix(source_vec, target_vec, method = "jw")
    closest_match_indices <- apply(dist_matrix, 1, which.min)
    closest_matches <- target_vec[closest_match_indices]
    return(closest_matches)
  }
  colnames(scaden.out) <- find_closest_match(colnames(scaden.out), unique(labels.in))


  return(t(scaden.out))
}
