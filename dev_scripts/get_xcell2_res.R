library(tidyverse)


# library(xCell2)
# source("/bigdata/almogangel/xCell2/R/xCell2Analysis.R")

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# Load reference-validation data (prep_ref_val_pairs.R)
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")


# Load references
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))

refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})

vals.refs.res <- refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())



# xCell 2.0 settings
thisseed <- 123
cores2use <- 40
objects_dir <- "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun_nofilt/"


# Make xCell2 objects
if (FALSE) {
  lapply(1:nrow(vals.refs.res), function(i){

    print(paste0("-------------------- ", i, "/", nrow(vals.refs.res), " --------------------"))

    # Load data
    val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
    print(val_ref)
    mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
    ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
    labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
    lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
    refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
    valType <- vals.refs.res[i,]$val_type
    valDataset <- vals.refs.res[i,]$val_dataset
    refName <- vals.refs.res[i,]$ref_name
    valDataType <- vals.refs.res[i,]$val_data_type[[1]]


    file <-  paste0(objects_dir, val_ref, ".rds")


    if (file.exists(file)) {
      xcell2_object <- readRDS(file)
    }else{

      # Load filtering data
      if (valDataType == "array") {
        filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/array_filtering_data.rds")
      }else{
        filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/rnaseq_filtering_data.rds")
      }


      # Remove current validation from filtering data
      filt_datasets <- gsub(x = names(filtering_data$mixture), pattern = "#.*", replacement = "")
      filtering_data$mixture <- filtering_data$mixture[filt_datasets != valDataset]
      filtering_data$truth <- filtering_data$truth[filt_datasets != valDataset]


      xcell2_object <-  xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, filtering_data = filtering_data, ref_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                            return_sigs_filt = FALSE, sigsFile = NULL, nCores = cores2use,  return_analysis = FALSE, predict_res = FALSE)


      saveRDS(xcell2_object, file)


    }

  })

}


# Get results using xCell2 objects
return_raw_scores <- TRUE
predict_res = FALSE
use_sillover = FALSE
cores2use <- 15
spillAlpha = 0.8

if (TRUE) {
all_res <- parallel::mclapply(1:nrow(vals.refs.res), function(i){

    print(paste0("-------------------- ", i, "/", nrow(vals.refs.res), " --------------------"))

    # Load data
    val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
    print(val_ref)
    mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
    valType <- vals.refs.res[i,]$val_type
    valDataset <- vals.refs.res[i,]$val_dataset
    refName <- vals.refs.res[i,]$ref_name


    file <-  paste0(objects_dir, val_ref, ".rds")
    if (file.exists(file)) {
      xcell2_object <- readRDS(file)
    }else{
      errorCondition("Missing file!")
    }

    res <- xCell2::xCell2Analysis(mix.in, xcell2object = xcell2_object, raw_scores = return_raw_scores,
                                  spillover = use_sillover, spillover_alpha = spillAlpha, num_threads = cores2use)
    return(res)

  }, mc.cores = 15)
vals.refs.res$res <- all_res
}

# saveRDS(vals.refs.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_4apr/xcell2.cyto.predict.res.rds")
saveRDS(vals.refs.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun_nofilt/raw_scores_ontology.rds")
