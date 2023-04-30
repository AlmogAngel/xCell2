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




dataFiles <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim", full.names = T)



cbrx.out.list <- lapply(dataFiles, function(file){

  print(file)
  data <- readRDS(file)
  mix <- data$bulk
  singleCellExpr <- data$singleCellExpr
  colnames(singleCellExpr) <- data$singleCellLabels
  refsample <- singleCellExpr

  runCIBERSORTx(mix, refsample, dir = "/bigdata/almogangel/CIBERSORTx_docker")
})


saveRDS(cbrx.out.list, "/bigdata/almogangel/twelve_years_decon_paper/analysis/cibersort_res_list.rds")
