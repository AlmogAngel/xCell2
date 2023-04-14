library(tidyverse)
source("R/utils.R")

makeCIBERSORTSigMat <- function(dir = "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/",
                                ref, pheno, method = "fractions", QN = FALSE){

  token <- "b72da36961922443b75a1b65beef27c0"

  ref_tmp <- paste0(dir, basename(ref), "-tmp")
  file.copy(ref, ref_tmp)

  pheno_tmp <- paste0(dir, basename(pheno), "-tmp")
  file.copy(pheno, pheno_tmp)

  # Make results directory
  results_dir <- paste0(dir, basename(ref), "-sigmat")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm ", results_dir, "/*"))
  }


  cmd <- paste0("docker run --platform linux/amd64 -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/",
                method, " --username almog.angel@campus.technion.ac.il --token ", token,
                " --refsample ", ref_tmp, " --phenoclasses ", pheno_tmp, " --QN ", QN,
                " 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  system(cmd, wait = TRUE)

  print("Done!")
}


runCIBERSORTx <- function(mix, sig_mat, dir = "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/",
                           method="fractions", single_cell=FALSE, QN=FALSE, absolute=TRUE, rmbatchBmode=TRUE, rmbatchSmode=FALSE){

  token <- "b72da36961922443b75a1b65beef27c0"


  mix_tmp <- paste0(dir, basename(mix), "-tmp")
  file.copy(mix, mix_tmp)

  # Make results directory
  results_dir <- paste0(dir, basename(mix), "-results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm ", results_dir, "/*"))
  }

  cmd <- paste0("docker run --platform linux/amd64 -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/",
                method, " --username almog.angel@campus.technion.ac.il --token ", token, " --single_cell ", single_cell,
                " --sigmatrix ", sig_mat, " --mixture ", mix_tmp, " --QN ", QN, " --absolute ",
                absolute, " --rmbatchBmode ", rmbatchBmode, " --rmbatchSmode ", rmbatchSmode,
                " 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  # Load results
  res_file <- ifelse(rmbatchBmode | rmbatchSmode, "/CIBERSORTx_Adjusted.txt", "/CIBERSORTx_Results.txt")
  #res_file <- "/CIBERSORTx_Results.txt"
  cibersortx_out <- t(read.table(paste0(results_dir, res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))

  unlink(mix_tmp)

  return(cibersortx_out)
}

getCIBERSORTxCor <- function(truth, cbrx_results){

  samples2use <- intersect(colnames(truth), colnames(cbrx_results))
  cbrx_results <- cbrx_results[,samples2use]
  truth <- truth[,samples2use]

  if (!all(colnames(truth) == colnames(cbrx_results))) {
    errorCondition(paste0("Error with dataset: ", file))
  }

  celltypes2use <- intersect(rownames(truth), rownames(cbrx_results))

  all_celltypes_cor <- sapply(celltypes2use, function(ctoi){

    scores_ctoi <- cbrx_results[ctoi,]
    truth_ctoi <- truth[ctoi, names(scores_ctoi)]
    truth_ctoi <- truth_ctoi[which(!is.na(truth_ctoi))]
    truth_ctoi <- truth_ctoi[which(truth_ctoi != "")]
    scores_ctoi <- scores_ctoi[names(truth_ctoi)]

    if (!all(names(truth_ctoi) == colnames(scores_ctoi))) {
      errorCondition(paste0("Error with dataset: ", file))
    }

    if (all(as.numeric(truth_ctoi) == 0) | all(as.numeric(scores_ctoi) == 0)) {
      NULL
    }else{
      cor(as.numeric(scores_ctoi), as.numeric(truth_ctoi), method = "spearman")
    }
  })

  out <- enframe(compact(all_celltypes_cor)) %>%
    arrange(-value) %>%
    rename(celltype = name , cor = value)

  return(out)

}

# Make Tabula Sapiens signature matrix for CIBERSORTx

ts_ref <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/ts_data.rds")
ts_labels <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/ts_labels.rds")

ts_ref.mat <-  cbind(rownames(ts_ref), as.matrix(ts_ref))
colnames(ts_ref.mat) <- c("GeneSymbol", ts_labels$label.main)


write_tsv(ts_ref.mat, "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_ts_ref.txt")



# # Make blood signature matrix File for CIBERSORTx ----
dep_list <- getDependencies(ontology_file_checked)

colnames(blood_ref) <- blood_labels$label
blood_ref %>%
  rownames_to_column(var = "genes") %>%
  write_tsv("/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_ref.txt")

all_celltypes <- unique(blood_labels$label)
pheno.df <- matrix(rep(1, ncol(blood_ref)*length(all_celltypes)), ncol = ncol(blood_ref))
pheno.df <- cbind(unique(blood_labels$label), pheno.df)

pheno.df <- t(sapply(all_celltypes, function(ct){
  dep_vec <- sapply(colnames(blood_ref), function(ref_label){
    ifelse(ct == ref_label, 1,
           ifelse(ref_label %in% dep_list[[ct]], 0, 2))
  })
  as.numeric(dep_vec)
}))
pheno.df %>%
  as.data.frame(.) %>%
  rownames_to_column() %>%
  write_tsv("/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_pheno_df.txt", col_names = FALSE)

makeCIBERSORTSigMat(ref = "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_ref.txt",
                    pheno = "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_pheno_df.txt")



# # Run CIBERSORTx with LM22 and BG_blood ----
# mix <- "/Users/almogang/Documents/xCell2.0/validation_data/expressions/BG_blood_expressions.tsv"
# sig_mat <- "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/LM22.txt"
#
#
# # Run CIBERSORTx
# bg_blood_lm22 <- runCIBERSORTx(mix = mix, sig_mat = sig_mat, absolute = FALSE, rmbatchBmode = TRUE)
# rownames(bg_blood_lm22) <- plyr::mapvalues(rownames(bg_blood_lm22), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
#
# # Combine scores to add broder cell types
# tmp <- t(data.frame("B-cells" = colSums(bg_blood_lm22[c("naive B-cells", "Memory B-cells", "Plasma cells"),]),
#              "CD4+ T-cells" = colSums(bg_blood_lm22[c("CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells"),]), check.names = FALSE))
# bg_blood_lm22 <- rbind(bg_blood_lm22, tmp)
#
# # Load truth
# truth_bg <- read.table(paste0("/Users/almogang/Documents/xCell2.0/validation_data/cell_values/BG_blood.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
# rownames(truth_bg) <- plyr::mapvalues(rownames(truth_bg), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
#
# getCIBERSORTxCor(truth = truth_bg, cbrx_results = bg_blood_lm22)
#
#
# # Run CIBERSORTx with Kassandra's blood ref and BG_blood ----
#
# mix <- "/Users/almogang/Documents/xCell2.0/validation_data/expressions/BG_blood_expressions.tsv"
# sig_mat <- "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/cibersortx_blood_ref.txt-sigmat/CIBERSORTx_cibersortx_blood_pheno_df.CIBERSORTx_cibersortx_blood_ref.bm.K999.txt"
#
#
# # Run CIBERSORTx
# bg_blood_kass <- runCIBERSORTx(mix = mix, sig_mat = sig_mat, absolute = FALSE, rmbatchBmode = TRUE)
# rownames(bg_blood_kass) <- plyr::mapvalues(rownames(bg_blood_kass), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
#
# # Load truth
# truth_bg <- read.table(paste0("/Users/almogang/Documents/xCell2.0/validation_data/cell_values/BG_blood.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
# rownames(truth_bg) <- plyr::mapvalues(rownames(truth_bg), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
#
# getCIBERSORTxCor(truth = truth_bg, cbrx_results = bg_blood_kass)
#
#
#
# # Run CIBERSORTx with all validation datasets -----
# mix_dir <- "/Users/almogang/Documents/xCell2.0/validation_data/expressions/"
# validation_ds <- gsub("_expressions.tsv", "", list.files(mix_dir))
# validation_ds_blood <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "GSE53655", "GSE60424", "sc_pbmc", "SDY67", "SDY420", "SDY311", "DREAM")
#
# cibersortx_res.list.lm22 <- list()
# for (mix_in in validation_ds_blood) {
#   use_qn <- ifelse(mix_in %in% c("SDY67", "SDY420", "SDY311"), TRUE, FALSE)
#   mix_file <- paste0(mix_dir, mix_in, "_expressions.tsv")
#   cbrx_out <- run_cibersortx(mix = mix_file, sig_mat = "/Users/almogang/Documents/xCell2.0/CIBERSORTx_docker/LM22.txt", QN=use_qn)
#   cibersortx_res.list.lm22[[mix_in]] <- cbrx_out
# }
#
# #saveRDS(cibersortx_res.list, "/Users/almogang/Documents/xCell2.0/cibersortx_res_list_kass_blood_ref.rds")
