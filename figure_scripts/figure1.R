library(tidyverse)
library(ggtext)

source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")


params2use <- list(
                   num_threads = 30,
                   min_pb_cells = 20,
                   min_pb_samples = 10,
                   min_sc_genes = 1e4,
                   use_ontology = TRUE,
                   return_signatures = FALSE,
                   return_analysis = FALSE,
                   use_sillover = TRUE,
                   spillover_alpha = 0.5
                   )


# B - Ontology vs No Ontology --------------------------------------


# Run xCell2 with ontology
dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1B/ontology"
xcell2_res_ontology <- paste0(dir2use, "/fig1b_xcell2_ontology_res.rds")
get_xcell2_benchmark_results(save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 4, output_name = xcell2_res_ontology)


# Run xCell2 without ontology
dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1B/no_ontology"
xcell2_res_no_ontology <- paste0(dir2use, "/fig1b_xcell2_no_ontology_res.rds")
params2use$use_ontology <- FALSE
get_xcell2_benchmark_results(save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 4, output_name = xcell2_res_no_ontology)


benchmark_correlations_ontology <- readRDS(xcell2_res_ontology)
benchmark_correlations_no_ontology <- readRDS(xcell2_res_no_ontology)

# Calculate correlations
benchmark_correlations_ontology <- get_xcell2_correlations(xCell2results = benchmark_correlations_ontology, round_results = 3)
benchmark_correlations_no_ontology <- get_xcell2_correlations(xCell2results = benchmark_correlations_no_ontology, round_results = 3)

# Prepare data
benchmark_correlations_ontology <- benchmark_correlations_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "yes")

benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology %>%
  filter(method == "xCell2") %>%
  mutate(ontology = "no")

# Remove cell type which the ontology have no effect
delta_cor_index <- (benchmark_correlations_ontology$ref_cor - benchmark_correlations_no_ontology$ref_cor) != 0
benchmark_correlations_ontology <- benchmark_correlations_ontology[delta_cor_index, ]
benchmark_correlations_no_ontology <- benchmark_correlations_no_ontology[delta_cor_index, ]

# Plot
benchmark_correlations <- rbind(benchmark_correlations_no_ontology, benchmark_correlations_ontology)
benchmark_correlations$ontology <- factor(benchmark_correlations$ontology, levels = c("no", "yes"))
fig1b <- ggplot(benchmark_correlations, aes(x=ontology, y=ref_cor, fill=ontology)) +
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


# fig1b: 575w x 600h
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1b.png", plot = fig1b, device = "png", width = 5.75, height = 6.00, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1b.pdf", plot = fig1b, device = "pdf", width = 5.75, height = 6.00)


# C - xCell (original) vs. xCell 2.0 --------------------------------------


getCors <- function(val_data, cts, sigs = NULL,
                    cyto.vals, is_xcell2,
                    xcell2objects_dir = "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1B/ontology"){
  
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
  
  ref_obj_files <- list.files(xcell2objects_dir, full.names = TRUE)
  ref_obj_files <- ref_obj_files[grepl("_bp.xcell2object.rds", ref_obj_files)]
  if (is.null(sigs) & is_xcell2) {
    xcell2_obj <- readRDS(ref_obj_files[grepl(val_data, ref_obj_files)])
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

celltype_conversion <- read_tsv("/bigdata/almogangel/xCell2_data/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

refs.vals.matched.human <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/refs.vals.matched.human.rds")
human.vals <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/human_validation.rds")


# Load xCell's Blueprint-ENCODE signatures
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

# Calculate correlations (xCell)
xcell1_cors <- refs.vals.matched.human %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, is_xcell2 = FALSE, sigs = xcell.sigs, cyto.vals = human.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell") %>%
  select(method, val_dataset, celltype, cor)

# Calculate correlations (xCell 2.0)
xcell2_cors <- refs.vals.matched.human %>%
  filter(ref_name == "bp") %>%
  rowwise() %>%
  mutate(cor = list(getCors(val_data = val_dataset, cts = shared_celltypes, is_xcell2 = TRUE, cyto.vals = human.vals))) %>%
  unnest_longer(cor, indices_to = "celltype") %>%
  mutate(method = "xCell 2.0") %>%
  select(method, val_dataset, celltype, cor)


# Plot
res_combined <- rbind(xcell1_cors, xcell2_cors)
fig1c <- ggplot(res_combined, aes(x=method, y=cor, fill=method)) +
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

# fig1c: 575w x 600h
print(fig1c)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1c.png", plot = fig1c, device = "png", width = 5.75, height = 6.00, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1c.pdf", plot = fig1c, device = "pdf", width = 5.75, height = 6.00)


# D - xCell 2.0 spillover analysis --------------------------------------


# Generate results using different alphas
alphas <- seq(0, 1, by = 0.25)
dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D"
params2use$use_ontology <- TRUE

for (alp in alphas) {
  print(paste0("Calculating results using alpha = ", alp))
  params2use$spillover_alpha <- alp
  xcell2_benchmark_res_file <- paste0(dir2use, "/fig1d_xcell2_res_alpha_", alp, ".rds")
  get_xcell2_benchmark_results(save_object = TRUE, dir = paste0(dir2use, "/alpha_", alp),
                               params = params2use, ncores = 4, output_name = xcell2_benchmark_res_file)
}


# Load results files
files <- list.files(dir2use, pattern = "fig1d_xcell2_res_alpha_*", full.names = TRUE)

# Calculate direct and spill correlation for all alpha values
spill_cor_res_table <- parallel::mclapply(files, function(file){
  
  alpha <- gsub(".rds", "", gsub("fig1d_xcell2_res_alpha_", "", basename(file),))
  refs.vals.matched.human.xcell2$res <- readRDS(file)
  
  parallel::mclapply(1:nrow(refs.vals.matched.human.xcell2), function(i){
    
    # Load ref, val, and ground truth
    valType <- refs.vals.matched.human.xcell2[i,]$val_type
    valDataset <- refs.vals.matched.human.xcell2[i,]$val_dataset
    refName <- refs.vals.matched.human.xcell2[i,]$ref_name
    refType <- refs.vals.matched.human.xcell2[i,]$ref_type
    ref.in <- human.refs[[refType]][[refName]]
    dep_list <- getDependencies(ref.in$lineage_file)
    gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
    cor_mat <- getCellTypeCorrelation(gep_mat, refType)
    truth <- human.vals$truth[[valType]][[valDataset]]
    truth <- round(truth, 3)
    res <- refs.vals.matched.human.xcell2[i,]$res[[1]]
    res <- round(res, 3)
    
    celltypes <- intersect(rownames(res), rownames(truth))
    samples <- intersect(colnames(res), colnames(truth))
    
    # Using correlation find the most similar cell type
    ct2most_simillar <- sapply(celltypes, function(ct){
      celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
      names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
    })
    
    y <- refs.vals.matched.human.xcell2[i,] %>%
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

# saveRDS(spill_cor_res_table, "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D/xcell2_spillover_correlations.rds")

# Plot
spill_cor_res_table <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D/xcell2_spillover_correlations.rds")
spill_cor_res_table[is.na(spill_cor_res_table$spill_cor),]$spill_cor <- 0
spill_cor_res_table[is.na(spill_cor_res_table$direct_cor),]$direct_cor <- 0

spill_cor_res_table <- spill_cor_res_table %>%
  mutate(delta_cor = direct_cor - spill_cor)
spill_cor_res_table <- spill_cor_res_table %>%
  select(-c(shared_celltypes, val_data_type))

spill_cor_res_table <- spill_cor_res_table %>%
  pivot_longer(cols = c(direct_cor, spill_cor), names_to = "cor_type", values_to = "cor")

fig1d <- ggplot(spill_cor_res_table, aes(x=spill_alpha, y=cor, fill=cor_type)) +
  coord_cartesian(ylim = c(-0.2, 1)) +
  geom_boxplot(width = .5, show.legend = TRUE, position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("direct_cor"="#155289", "spill_cor"="#B9DBF4"),
                    labels = c("Direct", "Spill")) +
  scale_y_continuous(breaks = seq(-0.2, 1, by = 0.1),
                     labels = function(y) {
                       ifelse(y == 0, "<b>0</b>", as.character(y))
                     }) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 14, face = "bold"),
        axis.text.y = element_markdown(size = 10), # Use element_markdown() to render bold text
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(color = ifelse(seq(-0.2, 1, by = 0.1) == 0, "black", "lightgray"),
                                          size = ifelse(seq(-0.2, 1, by = 0.1) == 0, 1, 0.3))) +
  labs(x="Alpha", y="Pearson r", fill="Correlation type")

# fig1d: 1200w x 600h
print(fig1d)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1d.png", plot = fig1d, device = "png", width = 12, height = 6.00, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig1d.pdf", plot = fig1d, device = "pdf", width = 12, height = 6.00)
