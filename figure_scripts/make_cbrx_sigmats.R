library(tidyverse)

makeCIBERSORTxSigMat <- function(ref, labels, lineage_file, refName, single_cell, QN,
                                 sigmat_out_dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_sigmats",
                                 gep_out_dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/cbrx_gep",
                                 docker_dir = "/bigdata/almogangel/CIBERSORTx_docker"){

  # Make pheno file
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
  dep_list <- getDependencies(lineage_file)
  all_celltypes <- unique(labels$label)
  pheno.mat <- t(sapply(all_celltypes, function(ct){
      ifelse(ct == labels$label, 1, 2)
    }))
  colnames(pheno.mat) <- labels$label
  # Insert zero to dependent cell types
  for (ct in all_celltypes) {
      dep_cts <- unname(unlist(dep_list[[ct]]))
      dep_cts <- dep_cts[dep_cts != ct]
      pheno.mat[ct, dep_cts] <- 0
    }
  pheno.df <- data.frame(cbind(rownames(pheno.mat), pheno.mat), check.names = FALSE)
  pheno_file <- paste0(docker_dir, "/", refName, "_pheno.txt")
  write.table(pheno.df, file = pheno_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


  # Make reference file
  tmp <- ref
  if (single_cell) {
    tmp <- as.matrix(tmp)
  }
  colnames(tmp) <- labels$label
  tmp <- cbind("genes" = rownames(tmp), tmp)
  all(colnames(tmp)[-1] == colnames(pheno.df)[-1])
  refsample_file <- paste0(docker_dir, "/", refName, "_ref.txt")
  write.table(tmp, file = refsample_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


  # Run CIBERSORTx
  token <- readLines("/bigdata/almogangel/xCell2_data/CIBERSORTx_token.txt", n = 1)

  # Make results directory
  results_dir <- paste0(docker_dir, "/results")
  if (!dir.exists(results_dir)){
    dir.create(results_dir)
  }

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm -f ", results_dir, "/*"))
  }


  cmd <- paste0("docker run -v ", docker_dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                token, " --refsample ", refsample_file, " --phenoclasses ", pheno_file, " --single_cell ", single_cell, " --QN ", QN, " --rmbatchBmode ", !single_cell, " --rmbatchSmode ", single_cell, " --verbose TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  sigmat_file <- paste0(sigmat_out_dir, "/", refName, "_sigmat.txt")
  sigmat_output_file <- list.files(results_dir, pattern = ".bm.K999.txt", full.names = TRUE)

  out_tmp <- read.table(sigmat_output_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  if (min(out_tmp) == 1) {
    out_tmp <- out_tmp -1
    out_tmp <- cbind("NAME" = rownames(out_tmp), out_tmp)
    write.table(out_tmp, sigmat_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }else{
    file.copy(from = sigmat_output_file, to = sigmat_file, overwrite = TRUE)
  }


  # Save GEP for future use:
  gep_file <- paste0(gep_out_dir, "/", refName, "_gep.txt")
  gep_output_file <- list.files(results_dir, pattern = "_sourceGEP.txt", full.names = TRUE)
  file.copy(from = gep_output_file, to = gep_file, overwrite = TRUE)


  print("Done")

  return(sigmat_file)
}


# Human
print("Kassandra Blood Reference")
kass_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")
makeCIBERSORTxSigMat(ref = kass_blood_ref$ref, labels = kass_blood_ref$labels, lineage_file = kass_blood_ref$lineage_file, refName = "kass_blood", single_cell = FALSE, QN = FALSE)
print("Done")

print("Kassandra Tumor Reference")
kass_tumor_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_ref.rds")
makeCIBERSORTxSigMat(kass_tumor_ref$ref, kass_tumor_ref$labels, kass_tumor_ref$lineage_file, refName = "kass_tumor", single_cell = FALSE, QN = FALSE)
print("Done")

print("BlueprintEncode Reference")
bp_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")
makeCIBERSORTxSigMat(bp_ref$ref, bp_ref$labels, bp_ref$lineage_file, refName = "bp", single_cell = FALSE, QN = FALSE)
print("Done")

print("LM22 Reference")
# lm22_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")
# makeCIBERSORTxSigMat(ref = lm22_ref$ref, labels = lm22_ref$labels, lineage_file = lm22_ref$lineage_file, refName = "lm22", single_cell = FALSE, QN = TRUE)
# print("Done")
# We will just load LM22 as is: # "/bigdata/almogangel/CIBERSORTx_docker/CIBERSORT_LM22_sigmat.txt"

print("Pan Cancer References")
sc_pan_cancer_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sc_pan_cancer_ref.rds")
makeCIBERSORTxSigMat(ref = sc_pan_cancer_ref$ref, labels = sc_pan_cancer_ref$labels, lineage_file = sc_pan_cancer_ref$lineage_file, refName = "sc_pan_cancer", single_cell = TRUE, QN = FALSE)
print("Done")

print("TS Blood")
ts_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")
makeCIBERSORTxSigMat(ref = ts_blood_ref$ref, labels = ts_blood_ref$labels, lineage_file = ts_blood_ref$lineage_file, refName = "ts_blood", single_cell = TRUE, QN = FALSE)
print("Done")



# Mouse
igd_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/igd_ref.rds")
makeCIBERSORTxSigMat(ref = igd_ref$ref, labels = igd_ref$labels, lineage_file = igd_ref$lineage_file, refName = "igd", single_cell = FALSE, QN = TRUE)
print("Done")


tm_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/tm_blood_ref.rds")
makeCIBERSORTxSigMat(ref = tm_blood_ref$ref, labels = tm_blood_ref$labels, lineage_file = tm_blood_ref$lineage_file, refName = "tm_blood", single_cell = TRUE, QN = FALSE)
print("Done")

mouse_rnaseq_data_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/mouse_rnaseq_data_ref.rds")
makeCIBERSORTxSigMat(ref = mouse_rnaseq_data_ref$ref, labels = mouse_rnaseq_data_ref$labels, lineage_file = mouse_rnaseq_data_ref$lineage_file, refName = "mouse_rnaseq_data", single_cell = FALSE, QN = FALSE)
print("Done")
