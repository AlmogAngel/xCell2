library(tidyverse)

makeCIBERSORTxSigMat <- function(ref, labels, lineage_file, refName, single_cell, QN, dir = "/bigdata/almogangel/CIBERSORTx_docker"){

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
  pheno_file <- paste0(dir, "/", refName, "_pheno.txt")
  write.table(pheno.df, file = pheno_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


  # Make reference file
  tmp <- ref
  if (single_cell) {
    tmp <- as.matrix(tmp)
  }
  colnames(tmp) <- labels$label
  tmp <- cbind("genes" = rownames(tmp), tmp)
  all(colnames(tmp)[-1] == colnames(pheno.df)[-1])
  refsample_file <- paste0(dir, "/", refName, "_ref.txt")
  write.table(tmp, file = refsample_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


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
                token, " --refsample ", refsample_file, " --phenoclasses ", pheno_file, " --single_cell ", single_cell, " --QN ", QN, " --rmbatchBmode ", !single_cell, " --rmbatchSmode ", single_cell, " --verbose TRUE 1> ",
                results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  # Save signature matrix for future use:
  sigmat_file <- paste0(dir, "/", refName, "_sigmat.txt")
  file.copy(from = paste0(results_dir, "/CIBERSORTx_", refName, "_pheno.CIBERSORTx_", refName, "_ref.bm.K999.txt"), to = sigmat_file, overwrite = TRUE)

  # Save GEP for future use:
  gep_file <- paste0(dir, "/", refName, "_gep.txt")
  file.copy(from = paste0(results_dir, "/CIBERSORTx_cell_type_sourceGEP.txt"), to = gep_file, overwrite = TRUE)


  print("Done")

  return(sigmat_file)
}


# RNA-seq references
# TODO: Add Mahmoud's reference
print("Kassandra Blood Reference")
kass_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")
makeCIBERSORTxSigMat(kass_blood_ref$ref, kass_blood_ref$labels, kass_blood_ref$lineage_file, refName = "kass_blood", single_cell = FALSE, QN = FALSE)
print("Done")

print("Kassandra Tumor Reference")
kass_tumor_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_tumor_ref.rds")
makeCIBERSORTxSigMat(kass_tumor_ref$ref, kass_tumor_ref$labels, kass_tumor_ref$lineage_file, refName = "kass_tumor", single_cell = FALSE, QN = FALSE)
print("Done")

print("BlueprintEncode Reference")
bp_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/bp_ref.rds")
makeCIBERSORTxSigMat(bp_ref$ref, bp_ref$labels, bp_ref$lineage_file, refName = "bp", single_cell = FALSE, QN = FALSE)
print("Done")



# Array references
print("LM22 Reference")
lm22_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/lm22_ref.rds")
makeCIBERSORTxSigMat(lm22_ref$ref, lm22_ref$labels, lm22_ref$lineage_file, refName = "lm22_adjusted", single_cell = FALSE, QN = TRUE)
print("Done")



# scRNA-seq references
print("Tabula Sapiens Blood Reference")
ts_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/ts_blood_ref.rds")
makeCIBERSORTxSigMat(ref = ts_blood_ref$ref, labels = ts_blood_ref$labels, lineage_file = ts_blood_ref$lineage_file, refName = "ts_blood", single_cell = TRUE, QN = FALSE)
print("Done")

