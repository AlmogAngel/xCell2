library(tidyverse)
library(xCell2)

# -------------- Single-cell validation --------------
data("ts_labels_with_ontology")

# xCell2 -----
allData <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/")
allData <- allData[allData != "Liver.rds"] # Because contain only two cell types

xcell2.out.list <- lapply(allData, function(data){

  data.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", data))
  ref <- data.in$singleCellExpr
  labels <- tibble(ont = "NA",
                   label = data.in$singleCellLabels,
                   sample = colnames(ref),
                   dataset = data.in$singleCellSubjects)

  # Get ontology
  labels <- labels %>%
    mutate(ont = plyr::mapvalues(sample, ts_labels_with_ontology$sample, ts_labels_with_ontology$ont.fine, warn_missing = FALSE))
  labels <- as.data.frame(labels)

  # Run xCell2.0
  xcell2Sigs <- xCell2Train(ref = ref, labels = labels, data_type = "sc")
  xcell2.out.mat <- xCell2Analysis(bulk = data.in$bulk, xcell2sigs = xcell2Sigs)
  rownames(xcell2.out.mat) <- gsub("-", "_", rownames(xcell2.out.mat))


  # Calculate correlation
  truth <- data.in$bulkRatio
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(xcell2.out.mat[x,], truth[x,], method = "spearman")
  })
})
names(xcell2.out.list) <- allData

xCell2.allCors.tidy <- enframe(xcell2.out.list, value = "cor", name = "tissue") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0) %>%
  mutate(method = "xCell2") %>%
  select(method, tissue, cor, celltype)

saveRDS(xCell2.allCors.tidy, "/bigdata/almogangel/twelve_years_decon_paper/analysis/xcell2_correlations_080523.rds")

# Load this
xCell2.allCors.tidy <- readRDS("/bigdata/almogangel/twelve_years_decon_paper/analysis/xcell2_correlations_080523.rds")



# Load single-cell validation results for all other methods -----
allRes <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/results/accuracy/")
allRes <- allRes[allRes != "Liver.rds"] # Because contain only two cell types

allCors <- sapply(allRes, function(res){

  # Load results
  res.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/results/accuracy/", res))

  # Calculate correlation
  prop <- t(res.in$P)
  truth <- res.in$groundTruth
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(prop[x,], truth[x,], method = "spearman")
  })
})

allCors.tidy <- enframe(allCors, value = "cor") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  separate(name, into = c("method", "tissue"), sep="_", extra = "merge") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0)

# CIBERSORTx -----

# Make results - only run once
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

allData <- list.files("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/")
allData <- allData[allData != "Liver.rds"] # Because contain only two cell types

cbrx.out.list <- lapply(allData, function(data){

  data.in <- readRDS(paste0("/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/", data))
  mix <- data.in$bulk
  singleCellExpr <- data.in$singleCellExpr
  colnames(singleCellExpr) <- data.in$singleCellLabels
  refsample <- singleCellExpr

  cbrx.out.mat <- runCIBERSORTx(mix, refsample, dir = "/bigdata/almogangel/CIBERSORTx_docker")

  # Calculate correlation
  truth <- data.in$bulkRatio
  celltypes <- rownames(truth)
  sapply(celltypes, function(x){
    cor(cbrx.out.mat[x,], truth[x,], method = "spearman")
  })
})

cbrx.out.list <- all_validation_combined_results$cbrx.out.list
names(cbrx.out.list) <- allData

CBRx.allCors.tidy <- enframe(cbrx.out.list, value = "cor", name = "tissue") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(tissue = gsub(".rds", "", tissue)) %>%
  replace(is.na(.), 0) %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, tissue, cor, celltype)

saveRDS(CBRx.allCors.tidy, "/bigdata/almogangel/twelve_years_decon_paper/analysis/cibersortx_correlations.rds")

# Load this
CBRx.allCors.tidy <- readRDS("/bigdata/almogangel/twelve_years_decon_paper/analysis/cibersortx_correlations.rds")

CBRx.allCors.tidy <- CBRx.allCors.tidy %>%
  filter(tissue != "Liver")

# Combine and plot all correlation results -----------------

allCors.tidy <- rbind(allCors.tidy, CBRx.allCors.tidy, xCell2.allCors.tidy)

allCors.tidy <- allCors.tidy %>%
  filter(tissue != "Liver")

allCors.tidy.median <- allCors.tidy %>%
  group_by(method, tissue) %>%
  summarise(cor = median(cor))

method_sorted <- allCors.tidy.median %>%
  group_by(method) %>%
  summarise(cor = median(cor)) %>%
  arrange(cor) %>%
  pull(method)

allCors.tidy.median$method <- factor(allCors.tidy.median$method, levels = method_sorted)

allCors.tidy.median %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=tissue), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median of all cell types)", title = "Tabula Sapiens Validation", x = NULL, colour = NULL, fill = NULL)


xCell2.allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=tissue), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_wrap(~tissue, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Tabula Sapiens Validation - xCell2", x = NULL, colour = NULL, fill = NULL)


rbind(allCors.tidy, CBRx.allCors.tidy) %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=method), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                              "#424242", "#00F5FF", "#FF3E96", "#54FF9F")) +
  facet_wrap(~tissue, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Tabula Sapiens Validation", x = NULL, colour = NULL, fill = NULL)

# -------------- Sorted-cells validation --------------

celltype_conversion_long <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Choose reference
kass_blood_ref <- readRDS("/bigdata/almogangel/xCell2_data/dev_data/kass_blood_ref.rds")
ref <- kass_blood_ref$ref
labels <- kass_blood_ref$labels
lineage_file <- "/bigdata/almogangel/xCell2_data/dev_data/kass_blood_dependencies_checked.tsv"

# Subset to non-depdendent cell types
celltypes2use <- c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "natural killer cell", "monocyte",
                   "neutrophil", "B cell", "basophil", "eosinophil")
labels <- labels[labels$label %in% celltypes2use,]
ref <- ref[,colnames(ref) %in% labels$sample]

# Choose validation mixtures
val_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
blood_val <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813")
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"



# xCell2 -----
# xcell2sigs <- xCell2Train(ref, labels, data_type = "rnaseq", lineage_file = lineage_file)
# saveRDS(xcell2sigs, "/bigdata/almogangel/xCell2_data/data/kass_blood_sigs_110523.rds")
xcell2sigs <- readRDS("/bigdata/almogangel/xCell2_data/data/kass_blood_sigs_110523.rds")


xcell2.out.list <- lapply(blood_val, function(val){

  print(val)
  mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
  xcell2.out.mat <- xCell2Analysis(bulk = mix, xcell2sigs = xcell2sigs)
  #rownames(xcell2.out.mat) <- gsub("-", "_", rownames(xcell2.out.mat))

  # Calculate correlation
  truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
  truth <- truth[rownames(truth) != "GSM",]
  rows <- rownames(truth)
  truth <- apply(truth, 2, as.numeric)
  rownames(truth) <- rows
  truth[is.na(truth)] <- 0
  rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  celltypes <- intersect(rownames(truth), rownames(xcell2.out.mat))
  samples <- intersect(colnames(truth), colnames(xcell2.out.mat))
  sapply(celltypes, function(x){
    cor(xcell2.out.mat[x, samples], truth[x, samples], method = "spearman")
  })
})
names(xcell2.out.list) <- blood_val


xCell2.allCors.tidy <- enframe(xcell2.out.list, value = "cor", name = "dataset") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  replace(is.na(.), 0) %>%
  mutate(method = "xCell2") %>%
  select(method, dataset, cor, celltype)


# CIBERSORTx  -----
runCIBERSORTx <- function(ref, labels, mix, lineage_file, refName, use_existing_sigmat = TRUE, dir = "/bigdata/almogangel/CIBERSORTx_docker"){


if(!use_existing_sigmat){

  # Make pheno file
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
  dep_list <- getDependencies(lineage_file)
  all_celltypes <- unique(labels$label)
  pheno.mat <- t(sapply(all_celltypes, function(ct){
    ifelse(ct == labels$label, 1, 2)
  }))
  colnames(pheno.mat) <- labels$label
  # Insert zero to dependent cell types
  # for (ct in all_celltypes) {
  #   dep_cts <- unname(unlist(dep_list[[ct]]))
  #   dep_cts <- dep_cts[dep_cts != ct]
  #   pheno.mat[ct, dep_cts] <- 0
  # }

  pheno.df <- data.frame(cbind(rownames(pheno.mat), pheno.mat), check.names = FALSE)
  pheno_file <- paste0(dir, "/", refName, "_pheno.txt")
  write.table(pheno.df, file = pheno_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



  # Make reference file
  tmp <- ref
  colnames(tmp) <- labels$label
  tmp <- cbind("genes" = rownames(tmp), tmp)
  #all(colnames(tmp)[-1] == colnames(pheno.df)[-1])
  refsample_file <- paste0(dir, "/", refName, "_ref.txt")
  write.table(tmp, file = refsample_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

}


  # Make mixture file
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

  if(!use_existing_sigmat){
    cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                  token, " --refsample ", refsample_file, " --phenoclasses ", pheno_file,  " --mixture ", mix_file, " --single_cell FALSE --QN FALSE --rmbatchBmode TRUE --verbose TRUE 1> ",
                  results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

    # Run Docker via shell
    system(cmd, wait = TRUE)

    # Save signature matrix for future use:
      file.copy(from = paste0(results_dir, "/CIBERSORTx_", refName, "_pheno.CIBERSORTx_", refName, "_ref.bm.K999.txt"), to = paste0(dir, "/", refName, "_sigmat.txt"), overwrite = TRUE)
  }else{
    # Load existing signature matrix
    sigmat_path <- paste0(dir, "/", refName, "_sigmat.txt")
    cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                  token, " --sigmatrix ", sigmat_path,  " --mixture ", mix_file, " --single_cell FALSE --QN FALSE --rmbatchBmode TRUE --verbose TRUE 1> ",
                  results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")

    # Run Docker via shell
    system(cmd, wait = TRUE)
  }




  # Load results
  cibersortx_out <- t(read.table(paste0(results_dir,  "/CIBERSORTx_Adjusted.txt"), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
  cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


  return(cibersortx_out)
}


# Run CIBERSORTx
cbrx.out.list <- lapply(blood_val, function(val){

  mix <- as.matrix(read.table(paste0(val_dir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))

  cibersortx_out <- runCIBERSORTx(ref, labels, mix, lineage_file, refName = "kass_blood", use_existing_sigmat = TRUE, dir = "/bigdata/almogangel/CIBERSORTx_docker")

  # Calculate correlation
  truth <- as.matrix(read.table(paste0(truths_dir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
  truth <- truth[rownames(truth) != "GSM",]
  rows <- rownames(truth)
  truth <- apply(truth, 2, as.numeric)
  rownames(truth) <- rows
  truth[is.na(truth)] <- 0
  rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  celltypes <- intersect(rownames(truth), rownames(cibersortx_out))
  samples <- intersect(colnames(truth), colnames(cibersortx_out))
  sapply(celltypes, function(x){
    cor(cibersortx_out[x, samples], truth[x, samples], method = "spearman")
  })
})
names(cbrx.out.list) <- blood_val

CBRx.allCors.tidy <- enframe(cbrx.out.list, value = "cor", name = "dataset") %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  replace(is.na(.), 0) %>%
  mutate(method = "CIBERSORTx") %>%
  select(method, dataset, cor, celltype)


# Run other methods  -----



# Combine and plot all correlation results -----------------

allCors.tidy <- rbind(CBRx.allCors.tidy, xCell2.allCors.tidy)


allCors.tidy.median <- allCors.tidy %>%
  group_by(method, dataset) %>%
  summarise(cor = median(cor))

method_sorted <- allCors.tidy.median %>%
  group_by(method) %>%
  summarise(cor = median(cor)) %>%
  arrange(cor) %>%
  pull(method)

allCors.tidy.median$method <- factor(allCors.tidy.median$method, levels = method_sorted)

allCors.tidy.median %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median of all cell types)", title = "Sorted Cells Validation - Blood", x = NULL, colour = NULL, fill = NULL)


xCell2.allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_wrap(~dataset, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Sorted Cells Validation - Blood (xCell2)", x = NULL, colour = NULL, fill = NULL)


allCors.tidy %>%
  ggplot(., aes(x=celltype, y=cor)) +
  geom_jitter(aes(col=method), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                              "#424242", "#00F5FF", "#FF3E96", "#54FF9F")) +
  facet_wrap(~dataset, scales = "free_x") +
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "#1A1A1A", linetype = "dashed"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Sorted Cells Validation - Blood", x = NULL, colour = NULL, fill = NULL)


