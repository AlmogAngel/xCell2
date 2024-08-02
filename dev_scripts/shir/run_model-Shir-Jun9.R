library(tidyverse)
# source("/home/sherlevy/shir/lightgbm_model_func.R")

number_of_iterations <- 50

# ------ load data and metadata -----
# >> Load tumors gene expression data
TPMs <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/TPMs.rds")
# metadata <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/totaldata_mini.rds") # Load metadata
# metadata <- metadata[, !(names(metadata) %in% c("ResponseNumeric", "isTop"))]
# categorical columns: Gender, treatment, Cancer_Type


#  --------------- Run xCell2 and other method (Almog Jun9) ---------------

source("/bigdata/almogangel/xCell2/dev_scripts/shir/get_other_methods_data.R")

# kass_tumor reference

# icb.dtangle <- rundtangle(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_tumor", refType = "rnaseq", suffix = "9jul")
# saveRDS(icb.dtangle, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.dtangle.rds")
# icb.kass_blood.dtangle <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.dtangle.rds")

# icb.mcpcounter <- runMCPcounter(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_tumor")
# saveRDS(icb.mcpcounter, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.mcpcounter.rds")
# icb.kass_blood.mcpcounter <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.mcpcounter.rds")

# print("Running DeconRNASeq...")
# icb.decon <- runDeconRNASeq(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_blood", refType = "rnaseq", dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul")
# saveRDS(icb.decon, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.decon.rds")
# icb.kass_blood.decon <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.decon.rds")

# print("Running BayesPrism...")
# icb.bayesprism <- runBayesPrism(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_blood", refType = "rnaseq", CPUs = 45, suffix = "9jul")
# saveRDS(icb.bayesprism, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.bayesprism.rds")
# icb.kass_blood.bayesprism <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.bayesprism.rds")

# print("Running EPIC...")
# icb.epic <- runEPIC(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_blood", refType = "rnaseq", dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul")
# saveRDS(icb.epic, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.epic.rds")

# print("Running CIBERSORTx...")
# icb.cbrx <- runCIBERSORTx(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "kass_blood", refType = "rnaseq", dir = "/bigdata/almogangel/CIBERSORTx_docker", suffix = "9jul")
# saveRDS(icb.cbrx, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.cbrx.rds")
# icb.kass_blood.cbrx <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_kass_tumor.9jul.cbrx.rds")



# sc_pan_cancer reference

# icb.dtangle <- rundtangle(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer", refType = "sc", suffix = "9jul")
# saveRDS(icb.dtangle, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.dtangle.rds")

# icb.mcpcounter <- runMCPcounter(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer")
# saveRDS(icb.mcpcounter, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.mcpcounter.rds")

# print("Running DeconRNASeq...")
# icb.decon <- runDeconRNASeq(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer", refType = "sc", dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul")
# saveRDS(icb.decon, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.decon.rds")

# print("Running BayesPrism...")
# icb.bayesprism <- runBayesPrism(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer", refType = "sc", CPUs = 25, suffix = "9jul")
# saveRDS(icb.bayesprism, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.bayesprism.rds")
# icb.sc_pan_cancer.bayesprism <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.bayesprism.rds")

# print("Running EPIC...")
# icb.epic <- runEPIC(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer", refType = "sc", dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/references", suffix = "9jul")
# saveRDS(icb.epic, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/icb_sc_pan_cancer.9jul.epic.rds")

print("Running CIBERSORTx...")
icb.cbrx <- runCIBERSORTx(mix = TPMs, mixType = "rnaseq", mixName = "icb", refName = "sc_pan_cancer", refType = "sc", dir = "/bigdata/almogangel/CIBERSORTx_docker", suffix = "9jul")
saveRDS(icb.cbrx, "/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/icb_sc_pan_cancer.9jul.cbrx.rds")
icb.cbrx <- readRDS("/bigdata/almogangel/xCell2_data/ICB_prediction/other_methods/9jul/icb_sc_pan_cancer.9jul.cbrx.rds")


# # ------ Run xCell2 ------
#
# xcell2_kass = readRDS('/bigdata/almogangel/Loaikathon/newRuns_0708_TPM_xcell2_kass_tumor.rds') %>% t() %>% as.data.frame()
# xcell2 = xcell2_kass
# colnames(xcell2) <- clean_names(colnames(xcell2))
# rownames(xcell2) <- clean_names(rownames(xcell2))
#
# # to use a different method replace xcellw with cibersort/immune_response_scores/quanTIseq_MINI/TIDE...
# method_scores = xcell2[rownames(metadata),] %>% as.data.frame()
#
# # order rows
# method_scores <- method_scores[rownames(methos_scores) %in% rownames(metadata), ]
# metadata <- metadata[rownames(metadata) %in% rownames(methos_scores), ]
# method_scores <- method_scores[order(match(rownames(methos_scores), rownames(metadata))), ]
# metadata <- metadata[order(match(rownames(metadata), rownames(methos_scores))), ]
#
# # ------ genes - extra data -----
# genes = readxl::read_xlsx('/bigdata/almogangel/xCell2/dev_scripts/loai/genes_to_use_for_ml.xlsx') # Genes as features
# wanted_genes = genes$genes
# genes = t(TPMs[wanted_genes, rownames(totaldata_mini)]) %>% as.data.frame()
# genes <-as.data.frame(lapply(genes, normalize))
# colnames(genes) = paste0("gene_", colnames(genes))
#
#
# # ------ run model -----
# # >> Set run options
# models.out = list() # List to hold model results
# add_scores = c(TRUE) # FALSE mean just metadata
# seeds = seq(2245, by=888, length.out=number_of_iterations)
#
# case <- 'CR'# For RECIST: "PD" "SD" "PR" "CR" or "NoResponse"/"Response"
#
# model_features = create_lightgbm_data(scores = method_scores,
#                                       metadata = metadata,
#                                       genes = genes,
#                                       outcome = case,
#                                       algo = 'xcell2') # stoped here!!! how to define the predictor
# counter <- 1
# res <- list()
# for (s in seeds) {
#   set.seed(s)
#   print(paste("Running with seed:", s))
#
#
#
#   res[[counter]] <- predict_response_lightgbm(data = model_features,
#                                                 iteration_num = 1,
#                                                 case = case,
#                                                 remove.sd = FALSE,
#                                                 metadata = metadata,
#                                                 algo = 'xcell2',
#                                                 num_threads = 1)
#   counter = counter + 1
# }
#
# auc_scores <- c()
# auc_scores <- sapply(res, function(x) x)
# print(auc_scores)
# print(mean(auc_scores))
