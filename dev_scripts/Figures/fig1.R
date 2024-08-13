library(tidyverse)

set.seed(123)

source("/bigdata/almogangel/xCell2/dev_scripts/benchmarking_functions.R")

dir2use <- "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/mouse/"


params2use <- list(human2mouse = FALSE,
                   num_threads = 30,
                   min_pb_cells = 20,
                   min_pb_samples = 10,
                   min_sc_genes = 1e4,
                   use_ontology = TRUE,
                   return_signatures = FALSE,
                   return_analysis = FALSE,
                   use_sillover = TRUE,
                   spillover_alpha = 0.5,
                   top_spill_value = 0.5)

roundResults <- 3

xcell2_benchmark_res_file <- paste0("xcell2_benchmark_results_spillAlpha_", params2use$spillover_alpha, ".rds")

get_xcell2_benchmark_results(vals2remove = c(), save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 5, output_name = xcell2_benchmark_res_file)


# B - Ontology --------------------------------------

benchmark_correlations_ontology <- readRDS(paste0(dir2use, xcell2_benchmark_res_file))
benchmark_correlations_no_ontology <- readRDS(paste0(dir2use, "no_ontology/", xcell2_benchmark_res_file))

benchmark_correlations_ontology <- get_xcell2_correlations(vals2remove = c(), xCell2results = benchmark_correlations_ontology, round_results = roundResults, weight_rhos = FALSE)
benchmark_correlations_no_ontology <- get_xcell2_correlations(vals2remove = c(), xCell2results = benchmark_correlations_no_ontology, round_results = roundResults, weight_rhos = FALSE)


benchmark_correlations_ontology <- benchmark_correlations_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "yes")

benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "no")

# Remove cell type which the ontology have no effect
delta_cor_index <- (benchmark_correlations_ontology$cor - benchmark_correlations_no_ontology$cor) != 0
benchmark_correlations_ontology <- benchmark_correlations_ontology[delta_cor_index, ]
benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology[delta_cor_index, ]

# Plot
benchmark_correlations <- rbind(benchmark_correlations_no_ontology, benchmark_correlations_ontology)
benchmark_correlations$ontology <- factor(benchmark_correlations$ontology, levels = c("no", "yes"))
ggplot(benchmark_correlations, aes(x=ontology, y=cor, fill=ontology)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  scale_fill_manual(values = c("yes" = "#155289", "no" = "#B9DBF4")) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  ggpubr::stat_compare_means(comparisons = list(c("yes", "no")), paired = TRUE,
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = 1) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(x="", y="Spearman rho", fill="Ontology?") +
  scale_x_discrete(expand=c(1, 0))




# B - xCell (original) vs. xCell 2.0 --------------------------------------

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

getCors <- function(val_data, cts, sigs = NULL,
                    cyto.vals, is_xcell2,
                    xcell2objects_dir = dir2use){

  # Load truth table
  truth_mat <- cyto.vals$truth$blood[[val_data]]
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$tumor[[val_data]]
  }
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$other[[val_data]]
  }
  truth_cts <- rownames(truth_mat)

  # Load mixture
  mix <- cyto.vals$mixtures$blood[[val_data]]
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$tumor[[val_data]]
  }
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$other[[val_data]]
  }
  mix_ranked <- singscore::rankGenes(mix)

  if (is.null(sigs) & is_xcell2) {
    xcell2_obj <- readRDS(paste0(xcell2objects_dir, val_data, "_bp.xcell2object.rds"))
    sigs <- xcell2_obj@signatures
    split_names <- strsplit(names(sigs), "#")
    named_list <- setNames(sigs, sapply(split_names, `[`, 1))
    sigs <- split(sigs, names(named_list))
  }

  ct2use <- intersect(cts, names(sigs))
  cts_cors <- parallel::mclapply(ct2use, function(ct){

    # Score sigs
    xcell.scores <- sapply(sigs[[ct]], simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    rownames(xcell.scores) <- colnames(mix_ranked)
    xcell.scores <- round(rowMeans(xcell.scores), 3)

    # Get shared samples
    samples <- intersect(colnames(truth_mat), names(xcell.scores))
    xcell.scores <- xcell.scores[samples]
    truth <- truth_mat[ct, samples]

    # Calculate correlation
    cor(truth, xcell.scores, method = "spearman", use = "pairwise.complete.obs")

  }, mc.cores = 20)
  cts_cors <- unlist(cts_cors)
  names(cts_cors) <- ct2use

  return(cts_cors)
}

# Load xCell signatures
xcell.data <- xCell::xCell.data
xcell.sigs <- list()
for (i in 1:length(xcell.data$signatures@.Data)) {
  signame <- GSEABase::setName(xcell.data$signatures@.Data[[i]])
  ct <- gsub("%.*","", signame)
  ref <- gsub(".*%(.+)%.*", "\\1", signame)

  if (!ref %in% c("BLUEPRINT", "ENCODE")) {
    next
  }

  genes <- GSEABase::geneIds(xcell.data$signatures@.Data[[i]])

  if (!ct %in% names(xcell.sigs)) {
    xcell.sigs[[ct]] <- list()
  }
  xcell.sigs[[ct]][[paste0("sig-", length(xcell.sigs[[ct]]))]] <- genes

}
names(xcell.sigs) <- plyr::mapvalues(names(xcell.sigs), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
xcell.cts <- names(xcell.sigs)

# Get xCell correlations
xcell1_cors <- refval.tbl %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, is_xcell2 = FALSE, sigs = xcell.sigs, cyto.vals = cyto.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell") %>%
  select(method, val_dataset, celltype, cor)

# Get xCell 2.0 correlations
xcell2_cors <- refval.tbl %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, is_xcell2 = TRUE, cyto.vals = cyto.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell 2.0") %>%
  select(method, val_dataset, celltype, cor)


# Plot
res_combined <- rbind(xcell1_cors, xcell2_cors)
res_combined %>%
  ggplot(., aes(x=method, y=cor, fill=method)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  scale_fill_manual(values = c("xCell 2.0" = "#155289", "xCell" = "#B9DBF4")) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  ggpubr::stat_compare_means(comparisons = list(c("xCell", "xCell 2.0")),
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = 1,
                             paired = TRUE) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(x="", y="Spearman rho", fill="") +
  scale_x_discrete(expand=c(1, 0))



# B - xCell 2.0 spillover analysis --------------------------------------

makeGEPMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  gep_mat <- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- Rfast::rowMedians(as.matrix(ref[,type_samples]))
    }
  })
  rownames(gep_mat) <- rownames(ref)

  return(gep_mat)
}
getCellTypeCorrelation <- function(gep_mat, ref_type){

  celltypes <- colnames(gep_mat)

  if (ref_type != "sc") {

    # Use top 10% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.9, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }else{

    # Use top 1% most variable genes
    genes_var <- apply(gep_mat, 1, var)
    most_var_genes_cutoff <- quantile(genes_var, 0.99, na.rm=TRUE)
    gep_mat <- gep_mat[genes_var > most_var_genes_cutoff,]

  }


  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  # TODO: Change for loop to apply function to measure time
  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(gep_mat[,celltype_i], gep_mat[,celltype_j], method = "spearman")
  }

  return(cor_mat)
}

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
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
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())


# Generate results using different alphas
alphas <- seq(0, 1, by = 0.25)

if (FALSE) {
    for (alp in alphas) {
      params2use$spillover_alpha <- alp
      xcell2_benchmark_res_file <- paste0("xcell2_benchmark_results_spillAlpha_", params2use$spillover_alpha, ".rds")
      get_xcell2_benchmark_results(vals2remove = c(), save_object = TRUE, dir = dir2use,
                                   params = params2use, ncores = 1, output_name = xcell2_benchmark_res_file)
    }
}


# Load results files
files <- list.files(dir2use, pattern = "*spillAlpha*", full.names = TRUE)

# Get xCell2 correlations
if (FALSE) {
  spill_cor_res_table <- parallel::mclapply(files, function(file){

    alpha <- gsub(".rds", "", gsub("xcell2_benchmark_results_spillAlpha_", "", basename(file),))
    vals.refs.res$res <- readRDS(file)

    parallel::mclapply(1:nrow(vals.refs.res), function(i){


      valType <- vals.refs.res[i,]$val_type
      valDataset <- vals.refs.res[i,]$val_dataset
      refName <- vals.refs.res[i,]$ref_name
      ref.in <- readRDS(paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/", refName, "_ref.rds"))
      dep_list <- getDependencies(ref.in$lineage_file)
      gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
      ref_type <- ifelse(refName %in% c("ts_blood", "sc_pan_cancer"), "sc", "rnaseq")
      cor_mat <- getCellTypeCorrelation(gep_mat, ref_type)
      truth <- cyto.vals$truth[[valType]][[valDataset]]
      truth <- round(truth, roundResults)
      res <- vals.refs.res[i,]$res[[1]]
      res <- round(res, roundResults)

      celltypes <- intersect(rownames(res), rownames(truth))
      samples <- intersect(colnames(res), colnames(truth))

      ct2most_simillar <- sapply(celltypes, function(ct){
        celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
        names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
      })

      y <- vals.refs.res[i,] %>%
        select(method:shared_celltypes)

      spill_cor_res <- tibble(sig_ct = celltypes, most_sim_truth_ct = ct2most_simillar)
      spill_cor_res <- cbind(y, spill_cor_res)

      spill_cor_res %>%
        rowwise() %>%
        mutate(spill_cor = cor(res[sig_ct, samples], truth[most_sim_truth_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        mutate(direct_cor = cor(res[sig_ct, samples], truth[sig_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
        ungroup() %>%
        return(.)


    }, mc.cores = 20) %>%
      bind_rows() %>%
      mutate(spill_alpha = alpha)

  }, mc.cores = 1) %>%
    bind_rows()
  saveRDS(spill_cor_res_table, "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2_spillover_correlations.rds")
}

spill_cor_res_table <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2_spillover_correlations.rds")


spill_cor_res_table[is.na(spill_cor_res_table$spill_cor),]$spill_cor <- 0
spill_cor_res_table[is.na(spill_cor_res_table$direct_cor),]$direct_cor <- 0

spill_cor_res_table <- spill_cor_res_table %>%
  mutate(delta_cor = direct_cor - spill_cor)
spill_cor_res_table <- spill_cor_res_table %>%
  select(-c(shared_celltypes, val_data_type))

# Plot
spill_cor_res_table %>%
  pivot_longer(cols = c(direct_cor, spill_cor), names_to = "cor_type", values_to = "cor") %>%
  ggplot(., aes(x=spill_alpha, y=cor, fill=cor_type)) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("direct_cor"="#155289", "spill_cor"="#B9DBF4"),
                    labels = c("Direct", "Spill")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 12, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12)) +
  labs(x="Alpha", y="Pearson r", fill="Correlation type")


# Get other method spill correlation
if (FALSE) {
  other_methods_spillcor <- get_xcell2_correlations(cMethod = "pearson", spillcors = TRUE, round_results = 3)
  saveRDS(other_methods_spillcor, "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/other_methods_spillover_correlations.rds")
}

other_methods_spillcor <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/other_methods_spillover_correlations.rds")


other_methods_spillcor[is.na(other_methods_spillcor$spill_cor),]$spill_cor <- 0
other_methods_spillcor[is.na(other_methods_spillcor$direct_cor),]$direct_cor <- 0

other_methods_spillcor[other_methods_spillcor$spill_cor < 0,]$spill_cor <- 0
other_methods_spillcor[other_methods_spillcor$direct_cor < 0,]$direct_cor <- 0


other_methods_spillcor <- other_methods_spillcor %>%
  mutate(delta_cor = direct_cor - spill_cor)
other_methods_spillcor$spill_alpha <- "0"



# Combine
spill_res_final <- rbind(spill_cor_res_table, other_methods_spillcor[,colnames(spill_cor_res_table)])
spill_res_final$spill_alpha <- factor(spill_res_final$spill_alpha)



# Plot

spill_res_final %>%
  filter(spill_alpha %in% c(0, 0.5)) %>%
  mutate(method = ifelse(method == "xCell2", paste0(method, " (alpha = ", spill_alpha, ")"), method)) %>%
  # mutate(method = factor(method, levels = c("xCell2 (alpha = 0)", "xCell2 (alpha = 0.1)", "xCell2 (alpha = 0.2)", "xCell2 (alpha = 0.3)", "xCell2 (alpha = 0.4)",
  #                                            "xCell2 (alpha = 0.5)", "xCell2 (alpha = 0.6)", "xCell2 (alpha = 0.7)", "xCell2 (alpha = 0.8)", "xCell2 (alpha = 0.9)", "xCell2 (alpha = 1)"))) %>%
  ggplot(., aes(x=method, y=delta_cor, fill=method)) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("spill_cor"="tomato", "direct_cor"="#7AC5CD"),
                    labels = c("Direct", "Spill")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),
        legend.position = "bottom") +
  labs(x="Alpha", y="Pearson r", fill="Correlation type")

data2plot <- spill_res_final %>%
  filter(spill_alpha %in% c(0, 0.5)) %>%
  mutate(method = ifelse(method == "xCell2", paste0(method, " (alpha = ", spill_alpha, ")"), method))

method_levels <- data2plot %>% group_by(method) %>% summarise(delta_cor = median(delta_cor)) %>% arrange(delta_cor) %>% pull(method)
method_levels[1:3] <- rev(method_levels[1:3])
data2plot$method <- factor(data2plot$method, levels = method_levels)
data2plot$is_xcell <- ifelse(startsWith(as.character(data2plot$method), "xCell2"), "yes", "no")

ggplot(data2plot, aes(x=method, y=delta_cor, fill=is_xcell)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  coord_flip(ylim = c(0, 0.8)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Delta Pearson r") +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10))




# Archive  ---------------------
source("/bigdata/almogangel/xCell2/dev_scripts/benchmarking_functions.R")

benchmark_correlations_ontology <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2_analysis_benchmark_output.rds")
benchmark_correlations_no_ontology <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/no_ontology/xcell2_analysis_benchmark_output.rds")

train_ds <- c()
benchmark_correlations_ontology <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = benchmark_correlations_ontology, round_xcell2_results = 3, weight_rhos = FALSE)
benchmark_correlations_no_ontology <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = benchmark_correlations_no_ontology, round_xcell2_results = 3, weight_rhos = FALSE)


benchmark_correlations_ontology <- benchmark_correlations_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "yes")

benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "no")


delta_cor_index <- (benchmark_correlations_ontology$cor - benchmark_correlations_no_ontology$cor) != 0

benchmark_correlations_ontology <- benchmark_correlations_ontology[delta_cor_index, ]
benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology[delta_cor_index, ]


benchmark_correlations <- rbind(benchmark_correlations_no_ontology, benchmark_correlations_ontology)

ggplot(benchmark_correlations, aes(x=ontology, y=cor, fill=ontology)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  #geom_jitter(alpha = 0.2, size = 0.5) +
  scale_fill_manual(values = c("yes" = "tomato", "no" = "gray")) +
  coord_cartesian(ylim = c(0.2, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  ggpubr::stat_compare_means(comparisons = list(c("yes", "no")),
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = 1) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_blank()) +
  labs(x="", y="Spearman Rho", fill="Ontology?")



# xCell1 vs 2
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


# Load xCell1 signatures
xcell.data <- xCell::xCell.data
xcell.sigs <- list()
for (i in 1:length(xcell.data$signatures@.Data)) {
  signame <- GSEABase::setName(xcell.data$signatures@.Data[[i]])
  ct <- gsub("%.*","", signame)
  ref <- gsub(".*%(.+)%.*", "\\1", signame)

  if (!ref %in% c("BLUEPRINT", "ENCODE")) {
    next
  }

  genes <- GSEABase::geneIds(xcell.data$signatures@.Data[[i]])

  if (!ct %in% names(xcell.sigs)) {
    xcell.sigs[[ct]] <- list()
  }
  xcell.sigs[[ct]][[paste0("sig-", length(xcell.sigs[[ct]]))]] <- genes

}
names(xcell.sigs) <- plyr::mapvalues(names(xcell.sigs), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
xcell.cts <- names(xcell.sigs)


# Get correlations

getCors <- function(val_data, cts, sigs = NULL,
                    cyto.vals, is_xcell2 = FALSE, xcell2objects_dir = "/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/"){

  # Load truth table
  truth_mat <- cyto.vals$truth$blood[[val_data]]
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$tumor[[val_data]]
  }
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$other[[val_data]]
  }
  truth_cts <- rownames(truth_mat)

  # Load mixture
  mix <- cyto.vals$mixtures$blood[[val_data]]
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$tumor[[val_data]]
  }
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$other[[val_data]]
  }
  mix_ranked <- singscore::rankGenes(mix)

  if (is.null(sigs) & is_xcell2) {
    xcell2_obj <- readRDS(paste0(xcell2objects_dir, val_data, "_bp.xcell2object.rds"))
    sigs <- xcell2_obj@signatures
    split_names <- strsplit(names(sigs), "#")
    named_list <- setNames(sigs, sapply(split_names, `[`, 1))
    sigs <- split(sigs, names(named_list))
  }

  ct2use <- intersect(cts, names(sigs))
  cts_cors <- parallel::mclapply(ct2use, function(ct){

    # Score sigs
    xcell.scores <- sapply(sigs[[ct]], simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    rownames(xcell.scores) <- colnames(mix_ranked)
    xcell.scores <- round(rowMeans(xcell.scores), 3)

    # Get shared samples
    samples <- intersect(colnames(truth_mat), names(xcell.scores))
    xcell.scores <- xcell.scores[samples]
    truth <- truth_mat[ct, samples]

    # Calculate correlation
    cor(truth, xcell.scores, method = "spearman", use = "pairwise.complete.obs")

  }, mc.cores = 20)
  cts_cors <- unlist(cts_cors)
  names(cts_cors) <- ct2use

  return(cts_cors)
}

xcell1_cors <- refval.tbl %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, sigs = xcell.sigs, cyto.vals = cyto.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell") %>%
  select(method, val_dataset, celltype, cor)


xcell2_cors <- refval.tbl %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, is_xcell2 = TRUE, cyto.vals = cyto.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell 2.0") %>%
  select(method, val_dataset, celltype, cor)


res_combined <- rbind(xcell1_cors, xcell2_cors)

res_combined %>%
  ggplot(., aes(x=method, y=cor, fill=method)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  scale_fill_manual(values = c("xCell 2.0" = "tomato", "xCell" = "#7AC5CD")) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_minimal() +
  ggpubr::stat_compare_means(comparisons = list(c("xCell", "xCell 2.0")),
                             method = "wilcox.test",
                             label = "p.format",
                             label.y = 1,
                             paired = TRUE) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_blank()) +
  labs(x="", y="Spearman Rho", fill="")

res_combined %>%
  filter(method == "xCell 2.0") %>%
  filter(cor < -0.5)



# Archive

# A - Signatures generation and ontology integration

# BG_Blood_ts_blood
truth_mat <- cyto.vals$truth$blood$BG_blood


###
# Signatures generation heatmap
###

# Required from xCell2Train():
mix
cor_mat
gep_mat
signatures


ctoi <- "B cell"
# ctoi <- "CD4-positive, alpha-beta T cell"

signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

length(signatures_ctoi)
# [1] 274
length(unique(unlist(signatures_ctoi)))
# [1] 359


# Get signatures scores with reference
gep_mat_ranked <- singscore::rankGenes(gep_mat[,names(sort(cor_mat[ctoi,], decreasing = T))])
gep_mat_ranked_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(gep_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
})
rownames(gep_mat_ranked_scores) <- colnames(gep_mat_ranked)


# Heatmap
color_palette <- colorRampPalette(c("darkred" ,"red", "white", "lightblue1", "darkblue"))(50)

ct2use <- rownames(gep_mat_ranked_scores)

data <- gep_mat_ranked_scores[ct2use, names(sort(gep_mat_ranked_scores[1,]))]

rownames(data) <- gsub("CD4-positive, alpha-beta", "CD4+", rownames(data))
rownames(data) <- gsub("CD8-positive, alpha-beta", "CD8+", rownames(data))
rownames(data) <- gsub("thymus-derived", "", rownames(data))


ht_list <- ComplexHeatmap::Heatmap(data, show_row_dend = FALSE, col = color_palette,
                                   show_column_dend = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
                                   row_names_side = "left", row_names_gp = gpar(fontface = "bold", fontsize = 10),
                                   heatmap_legend_param = list(title = "Enrichment score", direction = "horizontal",
                                                               title_position = "lefttop", at = seq(0, 1, 0.2), color_bar_len = 4))
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")




###
# BP heatmap
###

signatures_cts <- signatures[startsWith(names(signatures), paste0(ct2use, "#"))]
genes_cts <- unique(unlist(signatures_cts))
data <- ref[genes_cts, labels$label %in% ct2use]
data <- data[rowSums(data) >0,]

color_palette <- c(colorRampPalette(c("darkred" ,"red"))(5),
                   colorRampPalette(c("white", "lightblue1", "darkblue"))(50))

column_annotation_data <- data.frame(Cell_Type = labels[labels$sample %in% colnames(data),]$label)
group_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                  "#A65628", "#6959CD", "#999999", "darkblue", "darkred", "#8DA0CB",
                  "#E78AC3", "#54FF9F")
group_colors <- setNames(group_colors, ct2use)


annotations <- HeatmapAnnotation(
  df = column_annotation_data,
  col = list(Cell_Type = group_colors),
  simple_anno_size = unit(0.5, "cm")
)


ht_list <- ComplexHeatmap::Heatmap(data, show_row_dend = FALSE, col = color_palette, show_row_names = FALSE, top_annotation = annotations,
                                   show_column_dend = FALSE, cluster_columns = TRUE, show_column_names = FALSE,
                                   heatmap_legend_param = list(title = "log2(TPM)", direction = "horizontal",
                                                               title_position = "lefttop", color_bar_len = 4))
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


###
# Ontology boxplot
###

# 1) Run createSignatures with dep_list
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c_onto <- sort(c, decreasing = TRUE)

# 2) Run createSignatures with dep_list = NULL
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c_noOnto <- sort(c, decreasing = TRUE)


# 3) Run createSignatures with dep_list = NULL - only bad sigs
signatures_ctoi <- signatures[colnames(data)[1:50]]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c_noOntoBad <- sort(c, decreasing = TRUE)


# Plot results
data <- tibble(cors = c(c_onto, c_noOnto, c_noOntoBad), label = c(rep("Ontology", length(c_onto)),
                                                                  rep("No Ontology", length(c_noOnto)),
                                                                  rep("No Ontology (sub)", length(c_noOntoBad))))

data <- tibble(cors = c(c_onto, c_noOnto), label = c(rep("Ontology", length(c_onto)),
                                                     rep("No Ontology", length(c_noOnto))))


ggpubr::compare_means(cors~label, data, method = "wilcox.test")

data %>%
  ggplot(., aes(x=label, y= cors, fill=label)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "", x = "", y = "Spearman Rho with ground truth", fill = "") +
  scale_fill_manual(values = c("lightblue1", "#E9C61D", "tomato")) +
  scale_y_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.2)) +
  theme(legend.text = element_text(face = "bold", size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ggpubr::stat_compare_means(comparisons = list(c("Ontology", "No Ontology")), method = "wilcox.test", label = "p.signif", paired = FALSE)

###
# Ontology waterfall
###

# (1) get_xcell2_res.R - Make signatures
# (2) fig3.R - Make correlations
data.no = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_gen_tuning_3/onto/grid_26_cors_noonto.rds")

data.no <- data.no %>%
  filter(method == "xCell2") %>%
  mutate(with_onto = "no") %>%
  select(ref, val, celltype, cor, with_onto)


data.yes = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/sig_gen_tuning_3/onto/cors.rds")

data.yes <- data.yes %>%
  filter(method == "xCell2") %>%
  mutate(with_onto = "yes") %>%
  select(ref, val, celltype, cor, with_onto)


data <- data.yes %>%
  full_join(data.no, by = c("ref", "val", "celltype")) %>%
  #filter(cor.x >= 0.3 | cor.y >= 0.3) %>%
  mutate(delta_cor = cor.x - cor.y) %>%
  select(ref, val, celltype, delta_cor)

data_sorted <- data %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  filter(delta_cor != 0) %>%
  mutate(x_axis = row_number())


red_color_palette <- colorRampPalette(c("darkred", "red"))(33)
blue_color_palette <- colorRampPalette(c("lightblue1", "darkblue"))(34)
color_palette <- c(red_color_palette, blue_color_palette)


# Create a waterfall plot
data_sorted %>%
  ggplot(., aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = color_palette, guide = FALSE) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  annotate("text", x = 100, y = mean(data_sorted$delta_cor), label = paste0("Mean delta: ", round(mean(data_sorted$delta_cor), 2)), vjust = -0.5) +
  theme(axis.title.x = element_blank()) +
  labs(x = "", y = "Delta Correlation", title = "") +
  scale_y_continuous(limits = c(round(min(data_sorted$delta_cor), 1), round(max(data_sorted$delta_cor), 1)+0.1),
                     breaks = seq(round(min(data_sorted$delta_cor), 1), round(max(data_sorted$delta_cor), 1), 0.1)) +
  theme_minimal()




###
# Signatures filtering boxplot
###

# Require within filterSignatures function:
best_sigs
ds_cors_list

# Score input mixture and calculate correlations (make sure to run createSignatures with dep_list)
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c <- sort(c, decreasing = TRUE)


data <- tibble(sig = names(c), rho = c, filtered = "All")
data2 <- tibble(sig = best_sigs, rho = c[best_sigs], filtered = "Filtered")


# Run the complete filterSignatures to have the filtered signature including the essential genes
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c <- sort(c, decreasing = TRUE)


data3 <- tibble(sig = names(signatures_ctoi), rho = c, filtered = "Filtered + Essential")


data_mix <- rbind(data, data2, data3)
data_mix$data_type <- "Input Dataset"
data_mix$dataset <- "GSE127813"

data_mix <- data_mix %>%
  select(dataset, rho, sig, filtered, data_type)



# Filtering correlations (from the function)

data <- enframe(ds_cors_list, name = "dataset") %>%
  unnest_longer(value, values_to = "rho", indices_to = "sig") %>%
  mutate(filtered = "All")

data2 <- data %>%
  filter(sig %in% best_sigs) %>%
  mutate(filtered = "Filtered")

data_val <- rbind(data, data2)
data_val$data_type <- "Filtering Datasets"

# Combine filtering and mixture results
data_final <- rbind(data_mix, data_val)


# Boxplot
color_palette <- c("lightblue1", "tomato")
# data_final %>%
#   ggplot(., aes(x=dataset, y=rho, fill=filtered)) +
#   geom_boxplot(width = 0.8) +
#   facet_grid(~data_type, scales = "free_x", space = "free") +
#   theme_minimal() +  # A minimal theme for a nicer look
#   labs(title = "", x = "", y = "Spearman Rho", fill = "Signatures") +
#   scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +  # Set y-axis limits
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
#         axis.text.y = element_text(size = 10),
#         legend.position = "bottom",
#         axis.title.y = element_text(face = "bold", size = 10),
#         title = element_text(face = "bold", size = 12),
#         legend.text = element_text(face = "bold", size = 10),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         strip.text = element_text(face="bold", size=10),
#         axis.line = element_line(colour = "black"))

data_final %>%
  filter(data_type != "Input Dataset") %>%
  ggplot(., aes(x=dataset, y=rho, fill=filtered)) +
  geom_boxplot(width = 0.8) +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Spearman Rho", fill = "Signatures") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 10),
        title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face="bold", size=10),
        axis.line = element_line(colour = "black"))


color_palette <- c("lightblue1", "#E9C61D", "tomato")
data_final %>%
  filter(data_type == "Input Dataset") %>%
  ggplot(., aes(x=dataset, y=rho, fill=filtered)) +
  geom_boxplot(width = 0.8) +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Spearman Rho", fill = "Signatures") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, 0.1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        axis.title.y = element_text(face = "bold", size = 10),
        title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face="bold", size=10),
        axis.line = element_line(colour = "black"))


# Filterting boxplot

no_filt <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun_nofilt/data_combined.rds")
filt <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun/data_combined.rds")

no_filt <- no_filt %>%
  filter(method == "xCell2") %>%
  select(ref, ref_rho) %>%
  mutate(is_filt = "no")

filt <- filt %>%
  filter(method == "xCell2") %>%
  select(ref, ref_rho) %>%
  mutate(is_filt = "yes")

data <- rbind(filt, no_filt)


data %>%
  ggplot(., aes(x=ref, y=ref_rho, fill=is_filt)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0.2, 0.9))







###
# filtering waterfall
###

# (1) get_xcell2_res.R - Make signatures
# (2) fig3.R - Make correlations

data.no = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun_nofilt/cors.rds")

data.no <- data.no %>%
  filter(method == "xCell2") %>%
  mutate(with_filt = "no") %>%
  select(ref, val, celltype, cor, with_filt)


data.yes = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_objects_23jun/cors.rds")

data.yes <- data.yes %>%
  filter(method == "xCell2") %>%
  mutate(with_filt = "yes") %>%
  select(ref, val, celltype, cor, with_filt)


data <- data.yes %>%
  full_join(data.no, by = c("ref", "val", "celltype")) %>%
  filter(cor.x >= 0.3 | cor.y >= 0.3) %>%
  mutate(delta_cor = cor.x - cor.y) %>%
  select(ref, val, celltype, delta_cor)

data_sorted <- data %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  filter(delta_cor != 0) %>%
  mutate(x_axis = row_number())


red_color_palette <- colorRampPalette(c("darkred", "red"))(37)
blue_color_palette <- colorRampPalette(c("lightblue1", "darkblue"))(34)
color_palette <- c(red_color_palette, blue_color_palette)


# Create a waterfall plot
data_sorted %>%
  ggplot(., aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = color_palette, guide = FALSE) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  annotate("text", x = 100, y = mean(data_sorted$delta_cor), label = paste0("Mean delta: ", round(mean(data_sorted$delta_cor), 2)), vjust = -0.5) +
  theme(axis.title.x = element_blank()) +
  labs(x = "", y = "Delta Correlation", title = "") +
  scale_y_continuous(limits = c(round(min(data_sorted$delta_cor), 1), round(max(data_sorted$delta_cor), 1)),
                     breaks = seq(round(min(data_sorted$delta_cor), 1), round(max(data_sorted$delta_cor), 1), 0.1)) +
  theme_minimal()




# B - Signatures filtering


refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

getCors <- function(val_data, ct, sigs, cyto.vals){

  # Load truth table
  truth_mat <- cyto.vals$truth$blood[[val_data]]
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$tumor[[val_data]]
  }
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$other[[val_data]]
  }
  truth_cts <- rownames(truth_mat)

  # Load mixture
  mix <- cyto.vals$mixtures$blood[[val_data]]
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$tumor[[val_data]]
  }
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$other[[val_data]]
  }
  mix_ranked <- singscore::rankGenes(mix)

  # Score sigs
  xcell.scores <- sapply(sigs[[1]], simplify = TRUE, function(sig){
    singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })
  rownames(xcell.scores) <- colnames(mix_ranked)

  # Get shared samples
  samples <- intersect(colnames(truth_mat), rownames(xcell.scores))
  xcell.scores <- xcell.scores[samples,]
  truth <- truth_mat[ct, samples]

  # Calculate correlations
  xcell.sigs.cors <- apply(xcell.scores, 2, function(x){
    cor(truth, x, method = "spearman", use = "pairwise.complete.obs")
  })

  return(xcell.sigs.cors)
}

combineRhos2 <- function(rhos, sample_sizes = NULL){

  if (length(rhos) == 1) {
    return(rhos)
  }

  if (is.null(sample_sizes)) {
    sample_sizes <- rep(1, length(rhos))
  }


  if (length(sample_sizes) != length(rhos)) {
    # sample_sizes <- rep(sample_sizes, length(rhos))
    stop("values and weights must have the same length")
  }


  rhos[rhos == 1] <- 0.999999999
  rhos[rhos == -1] <- -0.999999999

  # Fisher's Z Transformation
  z_values <- 0.5 * log((1 + rhos) / (1 - rhos))


  weights <- sample_sizes
  z_mean <- sum(weights * z_values) / sum(weights)

  # Back Transformation
  rho_weighted_mean <- (exp(2 * z_mean) - 1) / (exp(2 * z_mean) + 1)
  return(rho_weighted_mean)

}


# Load xCell's BLUEPRINT-ENCODE signatures
celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

xcell.data <- xCell::xCell.data

xcell.sigs <- list()
for (i in 1:length(xcell.data$signatures@.Data)) {
  signame <- GSEABase::setName(xcell.data$signatures@.Data[[i]])
  ct <- gsub("%.*","", signame)
  ref <- gsub(".*%(.+)%.*", "\\1", signame)

  if (!ref %in% c("BLUEPRINT", "ENCODE")) {
    next
  }

  genes <- GSEABase::geneIds(xcell.data$signatures@.Data[[i]])

  if (!ct %in% names(xcell.sigs)) {
    xcell.sigs[[ct]] <- list()
  }
  xcell.sigs[[ct]][[paste0("sig-", length(xcell.sigs[[ct]]))]] <- genes

}

names(xcell.sigs) <- plyr::mapvalues(names(xcell.sigs), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
xcell.cts <- names(xcell.sigs)


# Get xCell2's signatures data (before/after filtering

# xcell2.sigs <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.filteredSigsCors.newEssentialv2valtype.3aprNewVal.rds")

# xcell2.sigs <- xcell2.sigs %>%
#   bind_rows() %>%
#   filter(ref == "bp") %>%
#   mutate(label = ifelse(is_filtered == "yes", "xCell 2.0 (filtered)", "xCell 2.0 (all)")) %>%
#   select(-c(ref, is_filtered))


all_cors # from figure 3

xcell2.sigs <- all_cors %>%
  bind_rows() %>%
  filter(ref == "bp" & method == "xCell2") %>%
  mutate(sig_name = "xCell2",
         label = "xCell 2.0") %>%
  select(val, celltype, sig_name, cor, label)


# Get correlations for first xCell version signatures

xcell.sigs.data <- xcell2.sigs %>%
  select(val, celltype) %>%
  unique() %>%
  rowwise() %>%
  mutate(xcell_sigs = list(xcell.sigs[[celltype]]))

xcell.sigs.data$cors <- NA
for (i in 1:nrow(xcell.sigs.data)) {

  val=xcell.sigs.data[i,]$val
  ct=xcell.sigs.data[i,]$celltype
  sigs=xcell.sigs.data[i,]$xcell_sigs
  cors <- getCors(val, ct, sigs, cyto.vals)
  xcell.sigs.data[i,]$cors <- list(cors)
}

xcell.sigs.data <- xcell.sigs.data %>%
  unnest_longer(cors, values_to = "cor", indices_to = "sig_name") %>%
  mutate(label = "xCell") %>%
  select(-xcell_sigs) %>%
  select(val, celltype, sig_name, cor, label)


# Plot results

results <- rbind(xcell.sigs.data, xcell2.sigs)

# Identify common val and celltype combinations
common_combinations <- reduce(list(xcell.sigs.data, xcell2.sigs),
                              ~inner_join(.x, .y, by = c("val", "celltype"))) %>%
  select(val, celltype) %>%
  distinct()


# Filter each tibble for common combinations and combine them
results <- list(xcell.sigs.data, xcell2.sigs) %>%
  map(~inner_join(.x, common_combinations, by = c("val", "celltype"))) %>%
  bind_rows()

color_palette <- c("lightblue1", "tomato")


# data <- results %>%
#   group_by(val, celltype, label) %>%
#   summarise(cor = mean(cor))
data <- results

pval <- ggpubr::compare_means(cor~label, data, method = "wilcox.test")


data %>%
  ungroup() %>%
  ggplot(., aes(x=label, y=cor, fill=label)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(data = data, comparisons = list(c("xCell", "xCell 2.0")), method = "wilcox.test", label = "p.format", paired = FALSE) +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Mean Spearman Rho", fill = "") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face="bold", size=10),
        axis.line = element_line(colour = "black"))


results %>%
  group_by(val, celltype, label) %>%
  summarise(cor = mean(cor)) %>%
  group_by(val, label) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=val, y=cor, fill=label)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Mean Spearman Rho", fill = "") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(0, 1.001), breaks = seq(0, 1, 0.1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),  # Make x-axis labels bold and bigger, adjust as needed
        legend.position = "right",  # Adjust legend position if needed
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        axis.line = element_line(colour = "black"))


### Heatmaps ###

celltypes <- results %>% select(val, celltype) %>%  unique() %>% pull(celltype) %>% table() %>% sort(.,decreasing = T) %>% names()
top_celltype <- celltypes[1:10]

xcell.hm <- results %>%
  filter(label == "xCell" & celltype %in% top_celltype) %>%
  group_by(val, celltype) %>%
  summarise(cor = mean(cor)) %>%
  pivot_wider(names_from = "val", values_from = cor)

ds <- colnames(xcell.hm)[-1]

cts <- xcell.hm$celltype
xcell.hm <- as.matrix(xcell.hm[,-1])
rownames(xcell.hm) <- cts

xcell.hm <- xcell.hm[top_celltype,ds]

color_palette <- c(colorRampPalette(c("darkred" ,"red", "#f07167", "white"))(40),
                   colorRampPalette(c("lightblue1", "blue", "darkblue"))(50))


ht_list <- ComplexHeatmap::Heatmap(xcell.hm, show_row_dend = FALSE, col = color_palette,
                                   show_column_dend = TRUE, cluster_columns = FALSE, show_column_names = TRUE,
                                   row_names_side = "left", row_names_gp = gpar(fontface = "bold", fontsize = 10),
                                   heatmap_legend_param = list(title = "Correlation", direction = "horizontal",
                                                               title_position = "lefttop", at = seq(-1, 1, 0.2), color_bar_len = 6))
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")






xcell2.hm <- results %>%
  filter(label == "xCell 2.0" & celltype %in% top_celltype) %>%
  group_by(val, celltype) %>%
  summarise(cor = mean(cor)) %>%
  pivot_wider(names_from = "val", values_from = cor)

cts <- xcell2.hm$celltype
xcell2.hm <- as.matrix(xcell2.hm[,-1])
rownames(xcell2.hm) <- cts

xcell2.hm <- xcell2.hm[top_celltype,ds]

ComplexHeatmap::Heatmap(xcell2.hm, show_row_dend = FALSE, col = color_palette,
                        show_column_dend = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
                        row_names_side = "left", row_names_gp = gpar(fontface = "bold", fontsize = 10), show_heatmap_legend = FALSE)


# Archive

###
# Regularization tuning
###

signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
scores_tmp <- (scores_tmp^(1/b)) / a
scores_tmp <- scale(scores_tmp)
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]

# Remove "i"
getCorss <- function(scores = scores_tmp, f = fracs, i, c){

  if (nrow(c) < 3) {
    return(NA)
  }
  s <- scores[,rownames(c)]
  p <- as.vector((s %*% c))
  cor <- cor(p, f, method = "spearman", use = "pairwise.complete.obs")
  return(cor)
}


reg_res <- lapply(seq(0, 1, 0.1), function(alp){

  cv_fit <- glmnet::cv.glmnet(X_scaled, Y, nfolds = nfold, grouped = grouped, alpha = alp, family = "gaussian")

  lambdas <- lapply(1:length(cv_fit$lambda), function(i){

    if (cv_fit$lambda[i] == cv_fit$lambda.1se) {
      is_1se <- "yes"
    }else{
      is_1se <- "no"
    }

    coefs <- coef(cv_fit, s = cv_fit$lambda[i])
    intercept <- coefs[1,]
    coefs <- as.matrix(coefs[which(coefs[-1, ] != 0) + 1,])
    tibble(lambda = cv_fit$lambda[i], intercept = intercept, reg_coef = list(coefs), is_1se = is_1se)
  }) %>%
    bind_rows()

  lambdas <- lambdas %>%
    rowwise() %>%
    mutate(cor = getCorss(i = intercept, c = reg_coef)) %>%
    drop_na() %>%
    mutate(alpha = alp)


  return(lambdas)

})
reg_res <- bind_rows(reg_res)


reg_res %>%
  mutate(alpha = as.character(alpha)) %>%
  ggplot(., aes(x=alpha, y=cor, fill=alpha)) +
  geom_boxplot()

reg_res %>%
  filter(is_1se == "yes")


c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
mean(c)















#
data.in = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.filteredSigsCors.newEssentialv2valtype.18mar.rds")
data.in = bind_rows(data.in)

data <- data.in %>%
  group_by(ref, val, celltype, is_filtered) %>%
  summarise(cor = mean(cor)) %>%
  # group_by(ref, val, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  # group_by(ref, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  spread(key = is_filtered, value = cor) %>%
  mutate(delta_cor = `yes` - `no`) %>%
  #select(ref, val, celltype, delta_cor) %>%
  #select(ref, val, delta_cor) %>%
  select(ref, delta_cor) %>%
  drop_na()


data_sorted <- data %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  mutate(x_axis = row_number())

# Create a waterfall plot
ggplot(data_sorted, aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "yellow", midpoint = 0) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Delta Cor", title = "Waterfall Plot of Delta Cor")










