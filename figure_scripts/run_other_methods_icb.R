library(tidyverse)
library(EPIC)
library(BayesPrism)
library(MCPcounter)
library(dtangle)
library(DeconRNASeq)
library(omnideconv)

tme_c.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/kass_tumor_ref.rds")
sc_pan_cancer.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")
icb_expr <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_expression_tpm.rds")

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


# bayesprism - DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$bayesprism$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$bayesprism$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$bayesprism$sc_pan_cancer))

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
  
  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/BayesPrism#", mixName, "#", refName, "#", suffix, ".rds")
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
runBayesPrism(mix = icb_expr, mixType = "rnaseq", mixName = "icb", refName = "kass_tumor", refType = "rnaseq", CPUs = 45, suffix = "9jul")


# cbrx
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$cbrx$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$cbrx$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$cbrx$sc_pan_cancer))

# decon - DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$decon$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$decon$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$decon$sc_pan_cancer))

getDeconRNASeqRes <- function(mix = icb_expr, mixName, refName, celltype_conversion, suffix = "9jul"){
  

  print(paste0(mixName, "_", refName))
  
  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/DeconRNASeq#", mixName, "#", refName, "#", suffix, ".rds")
  if (file.exists(results_file)) {
    print(paste0(mixName, "_", refName, " exist -  skipping...."))
    bayes_out <- readRDS(results_file)
    return(bayes_out)
  }
  
  
  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
  
  
  
  getResults <- function(cell_types, refName){
    
    sigmat <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,cell_types]
    
    if(max(sigmat) < 50){
      sigmat <- (2^sigmat)-1
    }
    
    
    # sigmat already normalized by CIBERSORTx
    
    # Run DeconRNASeq
    res <- DeconRNASeq(data.frame(mix), sigmat, checksig=FALSE, known.prop = FALSE,  fig = FALSE)
    res <- as.matrix(res$out.all)
    rownames(res) <- colnames(mix)
    res <- t(res)
    
    return(res)
  }
  
  
  # DeconRNASeq recommended not in log space
  if(max(mix) < 50){
    mix <- (2^mix)-1
  }
  
  
  # Get dependencies
  lineage_file <- ref.in$lineage_file
  if (!file.exists(lineage_file)) {
    stop("lineage file missing!")
  }
  dep_list <- getDependencies(lineage_file)
  
  
  # Get leaf cell types (no descendants)
  leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
  
  decon_out <- getResults(cell_types = leaf_cell_types, refName)
  
  celltypes2use <- unique(ref.in$labels$label)
  
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
  
  
  
  
  runDeconRNASeq(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes)
  
  return(vals.refs.res)
  
}
getDeconRNASeqRes(mix = icb_expr, mixName = "icb", refName = "kass_tumor", suffix = "9jul")

# epic - DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$epic$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$epic$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$epic$sc_pan_cancer))

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
  
  
  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/EPIC#", mixName, "#", refName, "#", suffix, ".rds")
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
runEPIC(mix = icb_expr, mixType = "rnaseq", mixName = "icb", refName = "kass_tumor", refType = "rnaseq", dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul")

# mcpcounter - OK
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$mcpcounter$tme_c))
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$mcpcounter$sc_pan_cancer))

# dtangle - OK
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$dtangle$tme_c))
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$dtangle$sc_pan_cancer))

# Scaden - DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$scaden$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$scaden$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$scaden$sc_pan_cancer)) ############################

runScaden <- function(mix, valType, valName, refName, refType){
  
  
  getResults <- function(cell_types, ref.in, mix, ref_name, val_dataset, idx){
    
    ct_index <- ref.in$labels$label %in% cell_types
    labels <- ref.in$labels$label[ct_index]
    ds <- ref.in$labels$dataset[ct_index]
    ref <- ref.in$ref[,ct_index]
    
    
    labels_tmp <- gsub("[^[:alnum:]]+", ".", labels)
    
    model_dir <- paste0("/bigdata/almogangel/scaden/", ref_name, "_", val_dataset, "_", idx)
    
    if (!dir.exists(model_dir)) {
      dir.create(model_dir)
      
      # Run Scaden via omnideconv
      scaden.model <- omnideconv::build_model_scaden(single_cell_object = as.matrix(ref),
                                                     bulk_gene_expression = mix,
                                                     cell_type_annotations = labels,
                                                     model_path = model_dir)
      
    }else{
      scaden.model <- model_dir
    }
    
    
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
    colnames(scaden.out) <- find_closest_match(colnames(scaden.out), unique(labels))
    scaden.out <- t(scaden.out)
    
    res <- scaden.out[cell_types,]
    
    return(res)
  }
  
  results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/Scaden#", valName, "#", refName, ".rds")
  if (file.exists(results_file)) {
    print(paste0(valName, "_", refName, " exist -  skipping...."))
    decon_out <- readRDS(results_file)
    return(decon_out)
  }
  
  if(max(mix) < 50){
    mix <- (2^mix)-1
  }
  
  
  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
  
  if (refType == "sc") {
    ref.in$ref <- apply(as.matrix(ref.in$ref), 2, function(x){x*10^6/sum(x)}) # Convert to CPM
  }
  
  # Get dependencies
  lineage_file <- ref.in$lineage_file
  if (!file.exists(lineage_file)) {
    stop("lineage file missing!")
  }
  dep_list <- getDependencies(lineage_file)
  
  # Get leaf cell types (no descendants)
  leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
  
  idx <- 1
  scaden_out <- getResults(cell_types = leaf_cell_types, ref.in, mix, ref_name = refName, val_dataset = valName, idx)
  celltypes2use <- unique(ref.in$labels$label)
  
  if (all(celltypes2use %in% rownames(scaden_out))) {
    
    
    scaden_out <- scaden_out[celltypes2use,]
    saveRDS(scaden_out, results_file)
    return(scaden_out)
    
  }else{
    
    # Get all missing cell type ancestors
    cell_types_left <- celltypes2use[!celltypes2use %in% rownames(scaden_out)]
    
    for (ancestor in cell_types_left) {
      idx <- idx + 1
      descendants <- dep_list[[ancestor]]$descendants
      ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
      ancestor_out <- getResults(cell_types = ct2use, ref.in, mix, refName, valName, idx)
      scaden_out <- rbind(scaden_out,  ancestor_out[ancestor,])
      rownames(scaden_out)[nrow(scaden_out)] <- ancestor
      
    }
    
    scaden_out <- scaden_out[celltypes2use,]
    saveRDS(scaden_out, results_file)
    return(scaden_out)
  }
  
}
runScaden(mix = icb_expr, valType = "rnaseq", valName = "icb", refName = "kass_tumor", refType = "rnaseq")


# dwls- DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$dwls$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$dwls$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$dwls$sc_pan_cancer))

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
dwls.out <- runDWLSRes(mix = icb_expr, refName = "kass_tumor")
results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/dwls#icb#kass_tumor#9jul.rds")
saveRDS(dwls.out, results_file)

# bisque - DONE
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$bisque$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$bisque$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$bisque$sc_pan_cancer))

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
bisque.out <- runBisqueRes(mix = icb_expr, refName = "kass_tumor")
results_file <- paste0("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/bisque#icb#kass_tumor#9jul.rds")
saveRDS(bisque.out, results_file)

# scdc
all(unique(tme_c.in$labels$label) %in% rownames(all_methods_final$scdc$tme_c))
unique(tme_c.in$labels$label)[!unique(tme_c.in$labels$label) %in% rownames(all_methods_final$scdc$tme_c)] # Missing kass_tumor ref
all(unique(sc_pan_cancer.in$labels$label) %in% rownames(all_methods_final$scdc$sc_pan_cancer))

runSCDC <- function(mix, valType, valName, refName, refType){
  
  
  getResults <- function(cell_types, ref.in, mix){
    
    ct_index <- ref.in$labels$label %in% cell_types
    labels <- ref.in$labels$label[ct_index]
    ds <- ref.in$labels$dataset[ct_index]
    ref <- ref.in$ref[,ct_index]
    
    # Run SCDC via omnideconv
    scdc.out <- omnideconv::deconvolute_scdc(single_cell_object = as.matrix(ref),
                                             bulk_gene_expression = mix,
                                             cell_type_annotations = labels,
                                             batch_ids = ds
    )
    
    
    res <- t(scdc.out$prop.est.mvw)
    
    return(res)
  }
  
  
  # SCDC recommended not in log space
  if(max(mix) < 50){
    mix <- (2^mix)-1
  }
  
  
  results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/SCDC#", valName, "#", refName, ".rds")
  if (file.exists(results_file)) {
    print(paste0(valName, "_", refName, " exist -  skipping...."))
    decon_out <- readRDS(results_file)
    return(decon_out)
  }
  
  ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
  
  if (refType == "sc") {
    ref.in$ref <- apply(as.matrix(ref.in$ref), 2, function(x){x*10^6/sum(x)}) # Convert to CPM
  }
  
  # Get dependencies
  lineage_file <- ref.in$lineage_file
  if (!file.exists(lineage_file)) {
    stop("lineage file missing!")
  }
  dep_list <- getDependencies(lineage_file)
  
  # Get leaf cell types (no descendants)
  leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
  
  scdc_out <- getResults(cell_types = leaf_cell_types, ref.in, mix)
  
  celltypes2use <- unique(ref.in$labels$label)
  
  if (all(celltypes2use %in% rownames(scdc_out))) {
    
    
    scdc_out <- scdc_out[celltypes2use,]
    saveRDS(scdc_out, results_file)
    return(scdc_out)
    
  }else{
    
    # Get all missing cell type ancestors
    cell_types_left <- celltypes2use[!celltypes2use %in% rownames(scdc_out)]
    
    for (ancestor in cell_types_left) {
      
      descendants <- dep_list[[ancestor]]$descendants
      ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
      ancestor_out <- getResults(cell_types = ct2use, ref.in, mix)
      scdc_out <- rbind(scdc_out,  ancestor_out[ancestor,])
      rownames(scdc_out)[nrow(scdc_out)] <- ancestor
      
    }
    
    scdc_out <- scdc_out[celltypes2use,]
    saveRDS(scdc_out, results_file)
    return(scdc_out)
  }
  
}
runSCDC(mix = icb_expr, refName = "kass_tumor", valName = "icb", valType = "rnaseq", refType = "rnaseq")


if (FALSE) {
  parallel_param <- MulticoreParam(workers = params$num_threads)
  
  # kass_tumor
  ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/kass_tumor_ref.rds")
  xcell2_object <- xCell2Train(ref = ref.in$ref, mix = TPMs, labels = ref.in$labels, refType = "rnaseq",
                               lineageFile = ref.in$lineage_file, BPPARAM = parallel_param, 
                               useOntology = params$use_ontolog, returnSignatures = params$return_signatures,
                               returnAnalysis = params$return_analysis, useSpillover = params$use_sillover,
                               spilloverAlpha = params$spillover_alpha, minPbCells = params$min_pb_cells,
                               minPbSamples = params$min_pb_samples, minScGenes = params$min_sc_genes)
  xcell2_res <- xCell2Analysis(mix = TPMs, xcell2object = xcell2_object, rawScores = !params$use_sillover,
                               spillover = params$use_sillover, spilloverAlpha = params$spillover_alpha, BPPARAM = parallel_param)
  saveRDS(xcell2_res, paste0(training_data_dir, "xcell2#icb#kass_tumor#9jul.rds"))
  
  # sc_pan_cancer
  ref.in <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")
  xcell2_object <- xCell2Train(ref = ref.in$ref, mix = TPMs, labels = ref.in$labels, refType = "sc",
                               lineageFile = ref.in$lineage_file, BPPARAM = parallel_param, 
                               useOntology = params$use_ontolog, returnSignatures = params$return_signatures,
                               returnAnalysis = params$return_analysis, useSpillover = params$use_sillover,
                               spilloverAlpha = params$spillover_alpha, minPbCells = params$min_pb_cells,
                               minPbSamples = params$min_pb_samples, minScGenes = params$min_sc_genes)
  xcell2_res <- xCell2Analysis(mix = TPMs, xcell2object = xcell2_object, rawScores = !params$use_sillover,
                               spillover = params$use_sillover, spilloverAlpha = params$spillover_alpha, BPPARAM = parallel_param)
  saveRDS(xcell2_res, paste0(training_data_dir, "xcell2#icb#sc_pan_cancer#9jul.rds"))
}


# Load data
icb_expr <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_expression_tpm.rds")
metadata <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_metadata.rds")
training_data_dir <- "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/"

files <- list.files("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul", full.names = TRUE)
other_methods <- lapply(files, function(f){readRDS(f)})
names(other_methods) <- gsub("#","-", gsub(".rds", "", basename(files)))
other_methods$`BayesPrism-icb-kass_blood-9jul`

# BayesPrism, CIBERSORTx, DeconRNASeq, dtangle

files2 <- list.files("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/", full.names = TRUE, pattern = "*rds")
other_methods2 <- lapply(files2, function(f){readRDS(f)})
names(other_methods2) <- gsub("#","-", gsub(".rds", "", basename(files2)))
other_methods2$icb_sc_pan_cancer.7aug.scaden

method2use <- c("bayesprism", "cbrx", "decon", "epic", "mcpcounter", "dtangle", "scaden", "dwls", "bisque", "scdc")

all_methods_final <- list()
all_methods_final <- lapply(method2use, function(m){all_methods_final[[m]]})
names(all_methods_final) <- method2use

all_methods_final$bayesprism <- list()
all_methods_final$bayesprism[["tme_c"]] <- other_methods$`BayesPrism-icb-kass_tumor-9jul`
dim(all_methods_final$bayesprism$tme_c)
all_methods_final$bayesprism[["sc_pan_cancer"]] <- other_methods$`BayesPrism-icb-sc_pan_cancer-9jul`
dim(all_methods_final$bayesprism$sc_pan_cancer)

all_methods_final$cbrx <- list()
all_methods_final$cbrx[["tme_c"]] <- other_methods$`CIBERSORTx-icb-kass_tumor`
dim(all_methods_final$cbrx$tme_c)
all_methods_final$cbrx[["sc_pan_cancer"]] <- other_methods$`CIBERSORTx-icb-sc_pan_cancer-9jul`
dim(all_methods_final$cbrx$sc_pan_cancer)

all_methods_final$decon <- list()
all_methods_final$decon[["tme_c"]] <- other_methods$`DeconRNASeq-icb-kass_tumor-9jul`
dim(all_methods_final$decon$tme_c)
all_methods_final$decon[["sc_pan_cancer"]] <- other_methods$`DeconRNASeq-icb-sc_pan_cancer-9jul`
dim(all_methods_final$decon$sc_pan_cancer)

all_methods_final$epic <- list()
all_methods_final$epic[["tme_c"]] <- other_methods$`EPIC-icb-kass_tumor-9jul`
dim(all_methods_final$epic$tme_c)
all_methods_final$epic[["sc_pan_cancer"]] <- other_methods$`EPIC-icb-sc_pan_cancer-9jul`
dim(all_methods_final$epic$sc_pan_cancer)

# TODO: missing cell types
all_methods_final$mcpcounter <- list()
all_methods_final$mcpcounter[["tme_c"]] <- other_methods2$icb_kass_tumor.7aug.mcpcounter
dim(all_methods_final$mcpcounter$tme_c)
all_methods_final$mcpcounter[["sc_pan_cancer"]] <- other_methods2$icb_sc_pan_cancer.7aug.mcpcounter
dim(all_methods_final$mcpcounter$sc_pan_cancer)

all_methods_final$dtangle <- list()
all_methods_final$dtangle[["tme_c"]] <- other_methods$`dtangle-icb-kass_tumor-9jul`
dim(all_methods_final$dtangle$tme_c)
all_methods_final$dtangle[["sc_pan_cancer"]] <- other_methods$`dtangle-icb-sc_pan_cancer-9jul`
dim(all_methods_final$dtangle$sc_pan_cancer)

# TODO: too many cell types
all_methods_final$scaden <- list()
all_methods_final$scaden[["tme_c"]] <- other_methods$`Scaden-icb-kass_tumor`
dim(all_methods_final$scaden$tme_c)
all_methods_final$scaden[["sc_pan_cancer"]] <- other_methods$`Scaden-icb-sc_pan_cancer`
dim(all_methods_final$scaden$sc_pan_cancer)

all_methods_final$dwls <- list()
all_methods_final$dwls[["tme_c"]] <- other_methods$`dwls-icb-kass_tumor-9jul`
dim(all_methods_final$dwls$tme_c)
all_methods_final$dwls[["sc_pan_cancer"]] <- other_methods2$icb_sc_pan_cancer.7aug.dwls
dim(all_methods_final$dwls$sc_pan_cancer)

all_methods_final$bisque <- list()
all_methods_final$bisque[["tme_c"]] <- other_methods$`bisque-icb-kass_tumor-9jul`
dim(all_methods_final$bisque$tme_c)
all_methods_final$bisque[["sc_pan_cancer"]] <- other_methods$`bisque-icb-sc_pan_cancer-9jul`
dim(all_methods_final$bisque$sc_pan_cancer)

all_methods_final$scdc <- list()
all_methods_final$scdc[["tme_c"]] <- other_methods$`SCDC-icb-kass_tumor`
dim(all_methods_final$scdc$tme_c)
all_methods_final$scdc[["sc_pan_cancer"]] <- other_methods2$icb_sc_pan_cancer.7aug.scdc
dim(all_methods_final$scdc$sc_pan_cancer)

length(all_methods_final)

all_methods_final$xcell2 <- list()
all_methods_final$xcell2[["tme_c"]] <- other_methods$`xcell2-icb-kass_tumor-9jul`
dim(all_methods_final$xcell2$tme_c)
all_methods_final$xcell2[["sc_pan_cancer"]] <- other_methods$`xcell2-icb-sc_pan_cancer-9jul`
dim(all_methods_final$xcell2$sc_pan_cancer)



x <- sort(table(unlist(lapply(all_methods_final, function(m){rownames(m$tme_c)}))), decreasing = TRUE)
sum(x == 11) # 12/25

scpcan_cts <- unique(sc_pan_cancer.in$labels$label) 
x2 <- sort(table(unlist(lapply(all_methods_final, function(m){rownames(m$sc_pan_cancer)[rownames(m$sc_pan_cancer) %in% scpcan_cts]}))), decreasing = TRUE)
sum(x2 == 11)

xcell2_res_tmec <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/mar_25/icb_kass_tumor.xcell2.rds")
dim(xcell2_res_tmec)

lapply(all_methods_final, function(m){all(rownames(xcell2_res_tmec) %in% rownames(m$tme_c))})

rownames(all_methods_final$bayesprism$tme_c)[!rownames(all_methods_final$bayesprism$tme_c) %in% rownames(xcell2_res_tmec)]

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

for (m in 1:length(all_methods_final)) {
  rownames(all_methods_final[[m]]$tme_c) <- plyr::mapvalues(rownames(all_methods_final[[m]]$tme_c) , celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
  rownames(all_methods_final[[m]]$sc_pan_cancer) <- plyr::mapvalues(rownames(all_methods_final[[m]]$sc_pan_cancer) , celltype_conversion$all_labels, celltype_conversion$xCell2_labels, warn_missing = FALSE)
  
}

saveRDS(all_methods_final, "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig3/icb_all_methods_predictions.rds")



