library(tidyverse)
source("/home/sherlevy/shir/lightgbm_model_func.R")

number_of_iterations <- 10

# ------ load data and metadata -----
# >> Load tumors gene expression data
TPMs <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/TPMs.rds")
metadata <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/totaldata_mini.rds") # Load metadata
#metadata <- metadata[, !(names(metadata) %in% c("ResponseNumeric", "isTop"))]
# categorical columns: Gender, treatment, Cancer_Type


# ------ Run xCell2 -----
xcell2_kass = readRDS('/bigdata/almogangel/Loaikathon/newRuns_0708_TPM_xcell2_kass_tumor.rds') %>% t() %>% as.data.frame()
xcell2 = xcell2_kass
colnames(xcell2) <- clean_names(colnames(xcell2))
rownames(xcell2) <- clean_names(rownames(xcell2))

tide = read.csv('/bigdata/loainaom/newRuns_0708/dataFILES/dataTIDE.csv', header = TRUE) %>% as.data.frame()
rownames(tide) <- tide[,1]

easier <- readRDS('/bigdata/almogangel/xCell2_data/ICB_prediction/icb_easier.rds') %>% t() %>% as.data.frame()
easier <- t(easier)
rownames(easier) <- easier[,1]
impers <- readRDS('/bigdata/almogangel/xCell2_data/ICB_prediction/icb_impres.rds') %>% t() %>% as.data.frame()
impers <- t(impers)

# to use a different method replace xcellw with cibersort/immune_response_scores/quanTIseq_MINI/TIDE...
method_scores = xcell2[rownames(metadata),] %>% as.data.frame()
method_scores <- tide[rownames(metadata),] %>% as.data.frame()
method_scores <- easier[rownames(metadata),] %>% as.data.frame()
method_scores <- impers[rownames(metadata),] %>% as.data.frame()

# order rows
method_scores <- method_scores[rownames(method_scores) %in% rownames(metadata), ]
metadata <- metadata[rownames(metadata) %in% rownames(method_scores), ]
method_scores <- method_scores[order(match(rownames(method_scores), rownames(metadata))), ]
metadata <- metadata[order(match(rownames(metadata), rownames(method_scores))), ]

# ------ genes - extra data -----
genes = readxl::read_xlsx('/bigdata/almogangel/xCell2/dev_scripts/loai/genes_to_use_for_ml.xlsx') # Genes as features
wanted_genes = genes$genes
genes = t(TPMs[wanted_genes, rownames(metadata)]) %>% as.data.frame()
genes <-as.data.frame(lapply(genes, normalize))
colnames(genes) = paste0("gene_", colnames(genes))


# ------ run model -----
# >> Set run options
models.out = list() # List to hold model results
add_scores = c(TRUE) # FALSE mean just metadata
seeds = seq(2245, by=888, length.out=number_of_iterations)

case <- 'PD'# For RECIST: "PD" "SD" "PR" "CR" or "NoResponse"/"Response"

model_features = create_lightgbm_data(scores = method_scores,
                                      metadata = metadata,
                                      genes = genes,
                                      outcome = case,
                                      algo = 'xcell2') # stoped here!!! how to define the predictor
counter <- 1
res <- list()
for (s in seeds) {
  set.seed(s)
  print(paste("Running with seed:", s))



  res[[counter]] <- predict_response_lightgbm(data = data,#model_features,
                                                iteration_num = 1,
                                                case = case,
                                                remove.sd = FALSE,
                                                metadata = metadata,
                                                algo = 'xcell2',
                                                num_threads = 1)
  counter = counter + 1
}

auc_scores <- c()
auc_scores <- sapply(res, function(x) x)
print(auc_scores)
print(mean(auc_scores))
