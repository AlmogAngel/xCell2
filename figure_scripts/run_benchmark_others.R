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

# Run CIBERSORTx
getCIBERSORTxRes <- function(ref_val_table, vals, celltype_conversion){
  
  
  runCIBERSORTx <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "/bigdata/almogangel/CIBERSORTx_docker"){
    
    
    getResults <- function(cell_types, dir, refName, valName, valType, single_cell, useQN){
      
      # Load custom CIBERSORTx signature matrix made in: make_cbrx_sigmats.R
      sigmat_tmp <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,cell_types])
      sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
      
      # Make signature matrix file
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
      # (!) CIBERSORTx work on non-log space
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
      
      cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                    token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell ," --rmbatchBmode ", !single_cell, " --QN ", useQN,
                    " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
      
      # Run Docker via shell
      system(cmd, wait = TRUE)
      
      
      # Load results
      res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
      cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
      cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]
      
      return(cibersortx_out)
    }
    
    print(paste0(valName, "_", refName))
    
    # Run settings
    single_cell <- ifelse(refType == "sc", TRUE, FALSE)
    useQN <- ifelse(refType == "array", TRUE, FALSE)
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/CIBERSORTx#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      cibersortx_out <- readRDS(results_file)
      return(cibersortx_out)
    }
    
    # Get dependencies
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
    
    lineage_file <- ref.in$lineage_file
    rm(ref.in)
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)
    
    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
    
    # Run CIBERSORTx for the first time
    cibersortx_out <- getResults(cell_types = leaf_cell_types, dir, refName, valName, valType, single_cell, useQN)
    
    # If cell type dependencies exist, run CIBERSORTx again with ancestors cell types until all cell types have been estimated
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
  
  
  # Run CIBERSORTx
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runCIBERSORTx(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "CIBERSORTx", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run EPIC
getEPICRes <- function(ref_val_table, vals, celltype_conversion){
  
  runEPIC <- function(vals, default_refs, valType, valName, refsRDSList, refName, refType, refTissue, celltypes2use, dir = "references"){
    
    
    getResults <- function(mix, cell_types, ref.in, refName, refType){
      
      # Using CRIBERSORTx's signature matrix
      # Subset sigGenes from the signature matrix
      # sigmat already normalized by CIBERSORTx
      sigmat <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      sigGenes <- rownames(sigmat)
      
      # Using CIBERSORTx's gene expression profile
      # Subset GEP
      gep <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_gep/", refName, "_gep.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      gep <- gep[,unique(cell_types)]
      if (min(gep) == 1) {
        gep <- gep-1
      }
      
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
      
      
      gep <- apply(gep, 2, function(x){x*10^6/sum(x)}) # Convert to CPM
      
      
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
      mix <- (2^mix)-1
    }
    
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/EPIC#", valName, "#", refName, ".rds")
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
    
    # Run EPIC
    epic_out <- getResults(mix, cell_types = leaf_cell_types, ref.in, refName, refType)
    
    
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
        ancestor_out <- getResults(mix, ct2use, ref.in, refName, refType)
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
  
  
  # i = 2
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
# Remember: Change clean genes to mm for mouse datasets
getBayesPrismRes <- function(ref_val_table, vals, celltype_conversion){
  
  runBayesPrism <- function(vals, valType, valName, refsRDSList, refName, refType, celltypes2use, CPUs = 45){
    
    dir.create(tempdir())
    
    getResults <- function(cell_types, ref.in, refType, CPUs){
      
      celltypeIndex <- ref.in$labels$label %in% cell_types
      ref.raw <- as.matrix(t(ref.in$ref[,celltypeIndex]))
      
      
      type <- ifelse(refType == "sc", "count.matrix", "GEP")
      ref.filtered <- cleanup.genes(input=ref.raw,
                                    input.type=type,
                                    species="hs", # Change to mm for mouse datasets
                                    gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") # from the vignette
                                    )
      
      
      ref.filtered.pc <- select.gene.type(ref.filtered,
                                          gene.type = "protein_coding")
      
      labels <- ref.in$labels[celltypeIndex,]$label
      
      if (refType == "sc") {
        diff.exp.stat <- get.exp.stat(sc.dat=ref.raw[,colSums(ref.raw>0)>3],# filter genes to reduce memory use (from vignette)
                                      cell.type.labels=labels,
                                      cell.state.labels=labels,
                                      # psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                      cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                      n.cores=CPUs) #number of threads
        
        
        ref.filtered.pc <- select.marker(sc.dat=ref.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)
        markers <- colnames(ref.filtered.pc)
        
        ref.filtered.pc <- apply(ref.raw, 2, function(x){x*10^6/sum(x)}) # Convert to CPM
        ref.filtered.pc <- ref.filtered.pc[,markers]
        
      }else{
        
        # Subset marker genes from the signature matrix
        markers <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
        shared_markers <- intersect(unique(markers$marker), colnames(ref.filtered.pc))
        ref.filtered.pc <- ref.filtered.pc[,shared_markers]
      }
      
      # Mark cancer cells labels
      tumorKey <- NULL
      if ("malignant cell" %in% labels) {
        tumorKey <- "malignant cell"
      }
      
      myPrism <- new.prism(
        reference=ref.filtered.pc,
        mixture=t(mix),
        input.type="GEP", # Also for scRNA-Seq because reference is now CPM
        cell.type.labels = labels,
        cell.state.labels = labels,
        key=tumorKey)
      
      bp.res <- run.prism(prism = myPrism, n.cores = CPUs)
      
      
      res <- t(get.fraction (bp=bp.res,
                             which.theta="final",
                             state.or.type="type"))
      
      return(res)
      
    }
    
    
    print(paste0(valName, "_", refName))
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/BayesPrism#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      bayes_out <- readRDS(results_file)
      return(bayes_out)
    }
    
    
    mix <- vals$mixtures[[valType]][[valName]]
    ref.in <- refsRDSList[[refType]][[refName]]
    
    # BayesPrism does not recommend log-space data
    if(max(mix) < 50){
      mix <- (2^mix)-1
    }
    
    if(max(ref.in$ref) < 50){
      ref.in$ref <- (2^ref.in$ref)-1
    }
    
    
    # Get cell types dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)
    all_celltypes <- names(dep_list)
    
    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
    
    bayes_out <- getResults(cell_types = leaf_cell_types, ref.in, refType, CPUs)
    
    
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
    markers <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/dtangle_markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
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
  
  rundtangle <- function(vals, valType, celltypes2use, refsRDSList, valName, valDataType, refName, refType, suffix = "8jan"){
    
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
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/dtangle#", valName, "#", refName, "#", suffix, ".rds")
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
  
  
  vals.refs.res <- ref_val_table %>%
    # Run dtangle
    rowwise() %>%
    mutate(res = list(rundtangle(vals, valType = val_type, valName = val_dataset, valDataType = val_data_type, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "dtangle", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run DeconRNASeq
getDeconRNASeqRes <- function(ref_val_table, vals, celltype_conversion){
  
  runDeconRNASeq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){
    
    
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
    
    
    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]
    
    # DeconRNASeq recommended not in log space
    if(max(mix) < 50){
      mix <- (2^mix)-1
    }
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/DeconRNASeq#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      decon_out <- readRDS(results_file)
      return(decon_out)
    }
    
    
    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    
    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)
    
    
    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
    
    decon_out <- getResults(cell_types = leaf_cell_types, refName)
    
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
  

  vals.refs.res <- ref_val_table %>%
    # Run DeconRNASeq
    rowwise() %>%
    mutate(res = list(runDeconRNASeq(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "DeconRNASeq", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run Bisque
getBisqueRes <- function(ref_val_table, vals, celltype_conversion){
  
  runBisqueSeq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){
    
    
    getResults <- function(cell_types, refName, ref.in, mix){
      
      ct_index <- ref.in$labels$label %in% cell_types
      labels <- ref.in$labels$label[ct_index]
      ds <- ref.in$labels$dataset[ct_index]
      ref <- ref.in$ref[,ct_index]
      markers <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
      markers <- unique(markers$marker)
      
      # Fix for:
      # Caused by error in `BisqueRNA::ReferenceBasedDecomposition()`:
      #! Only one individual detected in single-cell data. At least two subjects are needed (three or more recommended).
      
      if (length(unique(ds)) == 1) {
        ds[sample(1:length(ds), round(length(ds)/2))] <- paste0(unique(ds), "_2")
      }
      
      bisque.out <- omnideconv::deconvolute_bisque(single_cell_object = ref,
                                                   bulk_gene_expression = mix,
                                                   cell_type_annotations = labels,
                                                   batch_ids = ds,
                                                   markers = markers
      )
      
      res <- bisque.out$bulk.props
      return(res)
    }
    
    
    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]
    
    # Bisque recommended not in log space
    if(max(mix) < 50){
      mix <- (2^mix)-1
    }
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/Bisque#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      decon_out <- readRDS(results_file)
      return(decon_out)
    }
    
    
    ref.in <- refsRDSList[[refType]][[refName]]
    
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
    
    bisque_out <- getResults(cell_types = leaf_cell_types, refName, ref.in, mix)
    
    if (all(celltypes2use %in% rownames(bisque_out))) {
      
      
      bisque_out <- bisque_out[celltypes2use,]
      saveRDS(bisque_out, results_file)
      return(bisque_out)
      
    }else{
      
      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(bisque_out)]
      
      for (ancestor in cell_types_left) {
        
        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(cell_types = ct2use, refName, ref.in, mix)
        bisque_out <- rbind(bisque_out,  ancestor_out[ancestor,])
        rownames(bisque_out)[nrow(bisque_out)] <- ancestor
        
      }
      
      bisque_out <- bisque_out[celltypes2use,]
      saveRDS(bisque_out, results_file)
      return(bisque_out)
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
  

  vals.refs.res <- ref_val_table %>%
    # Run DeconRNASeq
    rowwise() %>%
    mutate(res = list(runBisqueSeq(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "Bisque", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run DWLS
getDWLSRes <- function(ref_val_table, vals, celltype_conversion){
  
  runDWLS <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){
    
    
    getResults <- function(cell_types, dwls.sigmat, mix){
      
      dwls.sigmat_tmp <- dwls.sigmat[,cell_types]
      
      solutions_ols <- parallel::mclapply(1:ncol(mix), function(i){
        mix <- mix[, i]
        DWLS::solveSVR(as.matrix(dwls.sigmat_tmp), mix)
        
      }, mc.cores = 20)
      names(solutions_ols) <- colnames(mix)
      
      res <- do.call(cbind, solutions_ols)
      
      return(res)
    }
    
    
    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]
    
    # Bisque recommended not in log space
    if(max(mix) < 50){
      mix <- (2^mix)-1
    }
    
    dwls.sigmat <- read.csv(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    
    genes <- intersect(rownames(dwls.sigmat), rownames(mix))
    mix <- mix[genes, , drop = FALSE]
    dwls.sigmat <- dwls.sigmat[genes, , drop = FALSE]
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/DWLS#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      decon_out <- readRDS(results_file)
      return(decon_out)
    }
    
    ref.in <- refsRDSList[[refType]][[refName]]
    
    # Get dependencies
    lineage_file <- ref.in$lineage_file
    if (!file.exists(lineage_file)) {
      stop("lineage file missing!")
    }
    dep_list <- getDependencies(lineage_file)
    rm(ref.in)
    
    # Get leaf cell types (no descendants)
    leaf_cell_types <- names(Filter(function(x) length(x$descendants) == 0, dep_list))
    
    dwls_out <- getResults(cell_types = leaf_cell_types, dwls.sigmat, mix)
    
    if (all(celltypes2use %in% rownames(dwls_out))) {
      
      
      dwls_out <- dwls_out[celltypes2use,]
      saveRDS(dwls_out, results_file)
      return(dwls_out)
      
    }else{
      
      # Get all missing cell type ancestors
      cell_types_left <- celltypes2use[!celltypes2use %in% rownames(dwls_out)]
      
      for (ancestor in cell_types_left) {
        
        descendants <- dep_list[[ancestor]]$descendants
        ct2use <- c(ancestor, leaf_cell_types[!leaf_cell_types %in% descendants])
        ancestor_out <- getResults(cell_types = ct2use, dwls.sigmat, mix)
        dwls_out <- rbind(dwls_out,  ancestor_out[ancestor,])
        rownames(dwls_out)[nrow(dwls_out)] <- ancestor
        
      }
      
      dwls_out <- dwls_out[celltypes2use,]
      saveRDS(dwls_out, results_file)
      return(dwls_out)
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
  

  vals.refs.res <- ref_val_table %>%
    # Run DWLS
    rowwise() %>%
    mutate(res = list(runDWLS(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "DWLS", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run SCDC
getSCDCRes <- function(ref_val_table, vals, celltype_conversion){
  
  runSCDC <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){
    
    
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
    
    
    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]
    
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
    
    ref.in <- refsRDSList[[refType]][[refName]]
    
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
  
  
  vals.refs.res <- ref_val_table %>%
    # Run SCDC
    rowwise() %>%
    mutate(res = list(runSCDC(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "SCDC", .before = everything())
  
  
  return(vals.refs.res)
  
}

# Run Scaden
getScadenRes <- function(ref_val_table, vals, celltype_conversion){
  
  runScaden <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){
    
    
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
    
    
    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]
    
    # Scaden recommended not in log space
    if(max(mix) < 50){
      mix <- (2^mix)-1
    }
    
    
    results_file <- paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/tmp/Scaden#", valName, "#", refName, ".rds")
    if (file.exists(results_file)) {
      print(paste0(valName, "_", refName, " exist -  skipping...."))
      decon_out <- readRDS(results_file)
      return(decon_out)
    }
    
    ref.in <- refsRDSList[[refType]][[refName]]
    
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
  
  # i = 2
  # x = ref_val_table
  # valType = x[i,]$val_type[[1]]; valName = x[i,]$val_dataset[[1]]; refName = x[i,]$ref_name[[1]]; refType = x[i,]$ref_type[[1]]; celltypes2use = x[i,]$celltype_classes[[1]]
  
  
  vals.refs.res <- ref_val_table %>%
    # Run SCDC
    rowwise() %>%
    mutate(res = list(runScaden(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "Scaden", .before = everything())
  
  
  return(vals.refs.res)
  
}


# Run other methods (Human/Mouse) -----------------------

# Load reference-validation pairs (prep_benchmark_data.R)
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.human.rds") # Human
# refval.tbl <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.mouse.rds") # Mouse


# Load validation data (prep_benchmark_data.R)
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds") # Human
# cyto.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/mouse_validation.rds") # Mouse

# Set references
# Human
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))
# Mouse
# refList <- list(rna_seq = c(mixed = "mouse_rnaseq_data"),
#                 array = c("igd"),
#                 sc = c(mixed = "tm_blood"))


celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


print("Running CIBERSORTx...")
cbrx.cyto.res <- getCIBERSORTxRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(cbrx.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/cbrx.cyto.res.rds")

print("Running EPIC...")
epic.cyto.res <- getEPICRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(epic.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/epic.cyto.res.rds")

print("Running BayesPrism...")
bp.cyto.res <- getBayesPrismRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(bp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/bp.cyto.res.rds")

print("Running MCPcounter")
mcp.cyto.res <- getMCPcounterRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(mcp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/mcp.cyto.res.rds")
# saveRDS(mcp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/mouse/mcp.cyto.res.rds")

print("Running dtangle")
dtan.cyto.res <- getdtangleRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(dtan.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.cyto.res.rds")

print("Running DeconRNASeq")
decon.cyto.res <- getDeconRNASeqRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(decon.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/decon.cyto.res.rds")

print("DWLS")
dwls.cyto.res <- getDWLSRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(dwls.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/dwls.cyto.res.rds")

print("Running Bisque")
bisque.cyto.res <- getBisqueRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(bisque.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/bisque.cyto.res.rds")
# saveRDS(bisque.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/mouse/bisque.cyto.res.rds")

print("Running SCDC")
scdc.cyto.res <- getSCDCRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(scdc.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/scdc.cyto.res.rds")

print("Running Scaden")
scaden.cyto.res <- getScadenRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(scaden.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/proportions/scaden.cyto.res")
