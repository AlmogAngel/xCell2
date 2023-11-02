library(tidyverse)
library(EPIC)
library(BayesPrism)
library(MCPcounter)
library(dtangle)
library(DeconRNASeq)


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
# Read reference-validation pairs
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
refval.tbl.nodeps <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val_nodeps.rds")
sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
sc.refval.tbl.nodeps <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val_nodeps.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")


# Method functions

# Run CIBERSORTx
getCIBERSORTxRes <- function(ref_val_table, vals, celltype_conversion){


  runCIBERSORTx <- function(vals, valType, valName, refName, refType, celltypes2use , dir = "/bigdata/almogangel/CIBERSORTx_docker"){

    print(paste0(valName, "_", refName))
    single_cell <- ifelse(refType == "sc", TRUE, FALSE)


    # Subset cell types from the signature matrix
    sigmat_tmp <- read.csv(paste0(dir, "/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
    sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
    write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

    if (single_cell) {
      # Subset reference file
      ref_tmp <- read.csv(paste0(dir, "/", refName, "_ref.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      index <- colnames(ref_tmp) %in% celltypes2use
      names <- colnames(ref_tmp)[index]
      ref_tmp <- cbind(rownames(ref_tmp), ref_tmp[,index])
      colnames(ref_tmp) <- c("genes", names)
      ref_tmp_file <- paste0(dir, "/ref-tmp.txt")
      write.table(ref_tmp, file = ref_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
    }

    # Make mixture file
    mix <- vals$mixtures[[valType]][[valName]]
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


    if (single_cell) {
      cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                    token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell , " --rmbatchSmode ", single_cell,
                    " --refsample ", ref_tmp_file, " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
    }else{
      cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                    token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell ," --rmbatchBmode ", !single_cell,
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


  # Add CIBERSORT default reference (LM22)
  lm22 <- read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sigmats/CIBERSORT_LM22_sigmat.txt", header = T, sep = "\t", check.names = F)

  ref_val_table <- ref_val_table %>%
    dplyr::select(val_type, val_dataset, n_val_samples) %>%
    unique() %>%
    mutate(ref_type = "default",
           ref_name = "CIBERSORT_LM22",
           ref_tissue = "mixed") %>%
    rowwise() %>%
    dplyr::select(-n_val_samples, everything()) %>%
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), colnames(lm22)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least 3 cell types
    filter(n_shared_celltypes > 2) %>%
    mutate(refsDeps = NA) %>%
    mutate(celltype_classes = list(shared_celltypes)) %>%
    rbind(., ref_val_table)

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

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    if (refType == "default") {

      res <- t(EPIC(bulk=mix, reference=default_refs[[refTissue]], withOtherCells=FALSE, scaleExprs=FALSE)$cellFractions)
      return(res)

    }


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
    res <- t(EPIC(bulk=mix, reference=epic_ref, withOtherCells=FALSE, scaleExprs=FALSE)$cellFractions)
    return(res)
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


  # Add EPIC default references (BRef/TRef)
  bref <- EPIC::BRef
  tref <- EPIC::TRef
  colnames(bref$refProfiles) <- plyr::mapvalues(colnames(bref$refProfiles), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(bref$refProfiles.var) <- plyr::mapvalues(colnames(bref$refProfiles.var), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(tref$refProfiles) <- plyr::mapvalues(colnames(tref$refProfiles), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(tref$refProfiles.var) <- plyr::mapvalues(colnames(tref$refProfiles.var), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  default_refs <- list("tumor" = tref, "blood" = bref)

  ref_val_table <- ref_val_table %>%
    dplyr::select(val_type, val_dataset, n_val_samples) %>%
    unique() %>%
    mutate(ref_type = "default",
           ref_name = ifelse(val_type == "tumor", "TRef", "BRef"),
           ref_tissue = ifelse(val_type == "tumor", "tumor", "blood")) %>%
    rowwise() %>%
    dplyr::select(-n_val_samples, everything()) %>%
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), colnames(default_refs[[val_type]]$refProfiles)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least 3 cell types
    filter(n_shared_celltypes > 2) %>%
    mutate(refsDeps = NA) %>%
    mutate(celltype_classes = list(shared_celltypes)) %>%
    rbind(., ref_val_table)


  # Run EPIC
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runEPIC(vals, default_refs, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, refTissue = ref_tissue, celltypes2use = celltype_classes))) %>%
    mutate(method = "EPIC", .before = everything())


  return(vals.refs.res)


}

# Run BayesPrism
getBayesPrismRes <- function(ref_val_table, vals, celltype_conversion){

  runBayesPrism <- function(vals, valType, valName, refsRDSList, refName, refType, celltypes2use, CPUs = 30){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- as.matrix(t(ref.in$ref[,celltypeIndex]))

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
                                    psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                    cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                    n.cores=CPUs) #number of threads


      ref.filtered.pc <- select.marker(sc.dat=ref.filtered.pc,
                                       stat=diff.exp.stat,
                                       pval.max=0.01,
                                       lfc.min=0.1)
    }else{

      # Subset marker genes from the signature matrix
      markers <- read.csv(paste0("references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
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
# TODO: try to fix patch for lm22
getdtangleRes <- function(ref_val_table, vals, celltype_conversion){

  rundtangle <- function(vals, valType, celltypes2use, refsRDSList, valName, refName, refType){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]

    if (refName == "lm22") {
      ref.in <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref_for_dtangle.rds")
    }

    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- ref.in$ref[,celltypeIndex]
    colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label

    # Normalize to CPM
    if (refType == "sc") {
      seurat_object <- Seurat::CreateSeuratObject(counts = ref.raw)
      seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "RC", scale.factor = 1e6)
      ref.cpm <- seurat_object@assays$RNA@data
      colnames(ref.cpm) <- colnames(ref.raw)
      ref.raw <- ref.cpm
    }


    celltypes <- unique(colnames(ref.raw))
    pure_samples_list <- lapply(celltypes, function(ct){
      which(colnames(ref.raw) == ct)
    })
    names(pure_samples_list) <- celltypes

    shared_genes <- intersect(rownames(mix), rownames(ref.raw))
    mix <- mix[shared_genes,]
    ref.raw <- ref.raw[shared_genes,]
    ref.raw <- as.matrix(ref.raw)

    if (refType == "sc") {
      y <- cbind(ref.raw, mix)
      y <- limma::normalizeBetweenArrays(y)
      ref.raw <- y[,1:ncol(ref.raw)]
      mix <- y[,(ncol(ref.raw)+1):ncol(y)]
    }


    markerMethod <- ifelse(refType == "sc", "ratio", "p.value")
    dataType <- ifelse(refType == "array", "microarray-gene", "rna-seq")
    res <- dtangle(Y = t(mix), references = t(ref.raw), pure_samples = pure_samples_list, marker_method = markerMethod, data_type = dataType)
    res <- t(res$estimates)

    return(res)
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


  # valType = x[1,]$val_type[[1]]; valName = x[1,]$val_dataset[[1]]; refName = x[1,]$ref_name[[1]]; refType = x[1,]$ref_type[[1]]; celltypes2use = x[1,]$celltype_classes[[1]]

  vals.refs.res <- ref_val_table %>%
    # Run dtangle
    rowwise() %>%
    mutate(res = list(rundtangle(vals, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "dtangle", .before = everything())




  return(vals.refs.res)

}

# Run DeconRNASeq
getDeconRNASeqRes <- function(ref_val_table, vals, celltype_conversion){

  runDeconRNASeq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,celltypes2use]

    # Run DeconRNASeq
    res <- DeconRNASeq(data.frame(mix), sigmat, checksig=FALSE, known.prop = FALSE,  fig = FALSE)
    res <- as.matrix(res$out.all)
    rownames(res) <- colnames(mix)
    res <- t(res)

    return(res)
  }


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
  runquanTIseq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]


    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,celltypes2use]

    # Run quanTIseq
    res <- t(quanTIseq(sigmat, mix, scaling=rep(1, ncol(sigmat)), method="lsei"))
    res <- res[rownames(res) != "Other",]


    return(res)
  }


  vals.refs.res <- ref_val_table %>%
    # Run quanTIseq
    rowwise() %>%
    mutate(res = list(runquanTIseq(vals, valType =  val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "quanTIseq", .before = everything())


  return(vals.refs.res)

}


# Run cytometry/other validations

print("Running CIBERSORTx...")
cbrx.cyto.res <- getCIBERSORTxRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(cbrx.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.cyto.res.rds")

print("Running EPIC...")
epic.cyto.res <- getEPICRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(epic.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.cyto.res.rds")

print("Running BayesPrism...")
bp.cyto.res <- getBayesPrismRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(bp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/bp.cyto.res.rds")

print("Running MCPcounter")
mcp.cyto.res <- getMCPcounterRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(mcp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mcp.cyto.res.rds")

print("Running dtangle")
dtan.cyto.res <- getdtangleRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(dtan.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.cyto.res.rds")

print("Running DeconRNASeq")
decon.cyto.res <- getDeconRNASeqRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(decon.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/decon.cyto.res.rds")

print("Running quanTIseq")
quanti.cyto.res <- getquanTIseqRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(quanti.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/quanti.cyto.res.rds")



# Run single-cell validations

print("Running CIBERSORTx...")
cbrx.sc.res <- getCIBERSORTxRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(cbrx.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.sc.res.rds")

print("Running EPIC...")
epic.sc.res <- getEPICRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(epic.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.sc.res.rds")

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

