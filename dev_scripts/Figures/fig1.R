############################################################
#
############################################################

library(tidyverse)
library(xCell2)
library(parallel)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)



setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))


# Load references matrices
refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
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

vals.refs.res.sc <- sc.refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(sc.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())

# Load other methods results
cyto.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".cyto.res.rds", full.names = TRUE), function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)


# xCell 2.0 settings
set.seed(123)
thisseed <- 123
cores2use <- 10
# gene settings
useTopVar <- TRUE
nTopVar <- 5000
ProteinCoding <- FALSE
ProteinCodingSC <- TRUE
genesGroups <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
genesGroupsSC <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY")
# signature settings
load_sigs <- FALSE
#sigs_suffix <- "x"
scores_results <- FALSE
minpbcells <- 30
minpbgroups <- 10
weight_genes <- TRUE
# simulations settings
simNoise <- NULL
fLvel <- "high"
# sim_method <- c("ref_multi", "ref_thin", "ref_mix_thin")
sim_method <- "ref_mix_thin"
simFracs <- c(0, seq(0.01, 0.25, 0.002), seq(0.3, 1, 0.05))
# xCell2Analysis
tranform <- TRUE
spillover <- FALSE
nSims <- 20


xCell2results <- lapply(1:nrow(vals.refs.res), function(i){

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


  # xCell2CleanGenes
  if (refType == "sc") {
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }else{
    shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
    ref <- shared_cleaned_genes$ref
    mix <- shared_cleaned_genes$mix
  }

  # Load signatures?
  if (load_sigs) {
    sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
  }else{
    sigsFile <- NULL
  }

  # xCell2Train
  if (sim_method == "ref_mix_thin") {
    mix2use <- mix
  }else{
    mix2use <- NULL
  }
  sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                               sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                               nCores = cores2use, simMethod = sim_method, mix = mix2use, sim_fracs = simFracs, filtLevel = fLvel, ct_sims = nSims)

  # xCell2Analysis
  res_mat <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
  res_mat <- res_mat[vals.refs.res[i,]$shared_celltypes[[1]], ]
  res_mat
})

saveRDS(xCell2results, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.fig1.cyto.res.rds")

f1 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.fig1.cyto.res.rds")

# allLevel1plots <- list()
# allLevel2plots <- list()
# for (i in 1:nrow(vals.refs.res)) {
#
#   # Load data
#   val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
#   print(val_ref)
#   mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
#   ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
#   labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
#   lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
#   refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
#   valType <- vals.refs.res[i,]$val_type
#   valDataset <- vals.refs.res[i,]$val_dataset
#   refName <- vals.refs.res[i,]$ref_name
#
#
#   # xCell2CleanGenes
#   if (refType == "sc") {
#     shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = useTopVar, use_protein_coding = ProteinCodingSC, n_var_genes = nTopVar, gene_groups = genesGroupsSC)
#     ref <- shared_cleaned_genes$ref
#     mix <- shared_cleaned_genes$mix
#   }else{
#     shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
#     ref <- shared_cleaned_genes$ref
#     mix <- shared_cleaned_genes$mix
#   }
#
#   # Load signatures?
#   if (load_sigs) {
#     sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
#   }else{
#     sigsFile <- NULL
#   }
#
#   # xCell2Train
#   if (sim_method == "ref_mix_thin") {
#     mix2use <- mix
#   }else{
#     mix2use <- NULL
#   }
#   sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
#                                sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
#                                nCores = cores2use, simMethod = sim_method, mix = mix2use, sim_fracs = simFracs, filtLevel = fLvel, ct_sims = nSims)
#
#   # xCell2Analysis
#   res_mat <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0)
#   res_mat <- res_mat[vals.refs.res[i,]$shared_celltypes[[1]], ]
#
#
#   # Load other methods results
#   cyto.Res.tmp <- cyto.Res %>%
#     filter(ref_name == refName, val_dataset == valDataset)
#   cyto.Res.tmp[cyto.Res.tmp$method == "xCell2",]$res[[1]] <- res_mat
#
#
#   # Calculate correlation with ground truth
#   truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
#
#
#   level1plots <- lapply(1:nrow(cyto.Res.tmp), function(j){
#
#     res <- cyto.Res.tmp[j,]$res[[1]]
#     celltypes <- intersect(rownames(res), rownames(truth_mat))
#     method <- cyto.Res.tmp[j,]$method[[1]]
#
#     df <- lapply(celltypes, function(ct){
#
#       truth <- truth_mat[ct,]
#       res <- res[ct,]
#
#       samples <- intersect(names(res), names(truth))
#
#       tibble(celltype = ct, truth = truth[samples], prediction = res[samples])
#
#     }) %>%
#       bind_rows()
#
#     cor_results <- df %>%
#       group_by(celltype) %>%
#       summarize(
#         cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
#         p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
#         n = sum(!is.na(truth) & !is.na(prediction))
#       ) %>%
#       mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
#       mutate(label = paste(" rho: ", round(cor, 3), "\n p: ", sprintf("%.2e", p_value), "\n n: ", n, sep = ""))
#
#     df <- df %>%
#       left_join(cor_results, by = "celltype")
#
#     # Plot correlations - level 1
#     p <- df %>%
#       drop_na() %>%
#       ggplot(., aes(x = truth, y = prediction, color = celltype)) +
#       geom_point(alpha = 0.7, size = 2) +
#       geom_smooth(method = "lm", se = TRUE, color = "black") +
#       facet_wrap(~celltype, scales = "free", strip.position = "top") +
#       theme_minimal(base_size = 14) +
#       theme(
#         plot.title = element_text(size = 16, face = "bold"),
#         plot.subtitle = element_text(size = 14),
#         strip.text = element_text(size = 14, face = "bold"),
#         strip.background = element_rect(fill = NA, colour = "black", size = 1),
#         panel.border = element_rect(color = "black", fill = NA),
#         panel.background = element_rect(fill = "gray95"),
#         legend.position = "none"
#       ) +
#       labs(title = method, subtitle = paste0("Ref: ", refName, ", Val: ", valDataset),
#            x = "Truth", y = "prediction") +
#       geom_text(data = cor_results, aes(label = label, x = -Inf, y = Inf),
#                 hjust = 0, vjust = 1.2, size = 4, color = "black")
#
#     p
#   })
#   names(level1plots) <- cyto.Res.tmp$method
#
#   if (!refName[[1]] %in% names(allLevel1plots)) {
#     allLevel1plots[[refName]] <- list()
#   }
#
#   if (!valDataset %in% names(allLevel1plots[[refName]])) {
#     allLevel1plots[[refName]][[valDataset]] <- list()
#   }
#
#   allLevel1plots[[refName]][[valDataset]] <- level1plots
#
#
#   getCors <- function(res, truth = truth_mat){
#
#     celltypes <- intersect(rownames(res), rownames(truth))
#
#     df <- lapply(celltypes, function(ct){
#
#       truth <- truth[ct,]
#       res <- res[ct,]
#
#       samples <- intersect(names(res), names(truth))
#
#       tibble(celltype = ct, truth = truth[samples], prediction = res[samples])
#
#     }) %>%
#       bind_rows()
#
#     df %>%
#       group_by(celltype) %>%
#       summarize(
#         cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
#         p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
#         n = sum(!is.na(truth) & !is.na(prediction))
#       ) %>%
#       mutate(cor = ifelse(is.na(cor), 0, cor))
#
#   }
#
#
#   cyto.Res.tmp <- cyto.Res.tmp %>%
#     rowwise() %>%
#     mutate(cors = list(getCors(res)))
#
#
#   df <- cyto.Res.tmp %>%
#     select(method, cors) %>%
#     unnest(cols = c(cors))
#
#   # Reordering the factor levels based on median of correlation for each method
#   df <- df %>%
#     mutate(method = factor(method, levels = names(sort(tapply(cor, method, median, na.rm = TRUE), decreasing = TRUE))))
#
#   # Adding a column to specify color based on method
#   df <- df %>%
#     mutate(isxcell = ifelse(method == "xCell2", "tomato", "gray"))
#
#   mincor <- min(df$cor)-0.1
#
#   # Defining distinct colors for cell types
#   num_colors <- length(unique(df$celltype))
#   color_palette <- brewer.pal(num_colors, "Set2")
#
#   # Plotting
#   p2 <- ggplot(df, aes(x = method, y = cor)) +
#     geom_boxplot(aes(fill = isxcell), outlier.shape = NA) +
#     geom_point(aes(color = celltype), position = position_dodge(width = 0.75), size = 3) +
#     scale_fill_manual(values = c("gray", "tomato")) +
#     scale_color_manual(values = color_palette) +
#     theme_minimal() +
#     labs(x = "", y = "Spearman Rho", color = "") +
#     scale_y_continuous(limits = c(mincor, 1), breaks = seq(-1, 1, 0.2)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, face = ifelse(levels(df$method) == "xCell2", "bold", "plain"))) +
#     theme(legend.position = "right") +
#     theme(legend.key.width = unit(1, "cm")) +
#     theme(legend.text = element_text(size = 10)) +
#     theme(legend.title = element_text(size = 10, face = "bold")) +
#     guides(fill = FALSE)
#
#   # Reordering the factor levels based on median of correlation for each celltype
#   df <- df %>%
#     mutate(celltype = factor(celltype, levels = names(sort(tapply(cor, celltype, median, na.rm = TRUE), decreasing = FALSE))))
#
#   mincor <- min(df$cor)-0.1
#   color_palette[levels(df$method) == "xCell2"] <- "red"
#
#   p3 <- ggplot(df, aes(x = celltype, y = cor)) +
#     geom_boxplot(outlier.shape = NA, fill = "lightgray") +
#     geom_point(aes(color = method, shape = isxcell), position = position_dodge(width = 0.75), size = 3) +
#     scale_color_manual(values = color_palette) +
#     scale_shape_manual(values = c("gray" = 16, "tomato" = 8)) +  # Assigns shapes based on isxcell, using backticks for TRUE and FALSE
#     theme_minimal() +
#     labs(x = "", y = "Spearman Rho", color = "") +
#     scale_y_continuous(limits = c(mincor, 1), breaks = seq(-1, 1, 0.2)) +
#     theme(axis.text.x = element_text(hjust = 1)) +
#     theme(legend.position = "right") +
#     theme(legend.key.width = unit(1, "cm")) +
#     theme(legend.text = element_text(size = 10)) +
#     theme(legend.title = element_text(size = 10, face = "bold")) +
#     guides(shape = FALSE) +
#     coord_flip()
#
#   refvalname <- paste0("Ref: ", refName, ", Val: ", valDataset)
#
#   p23 <- gridExtra::grid.arrange(p2, p3, ncol = 1,
#                                  top = textGrob(refvalname,
#                                                 gp = gpar(fontface = "bold", fontsize = 14)))
#
#
#   if (!refvalname %in% names(allLevel2plots)) {
#     allLevel2plots[[refvalname]] <- list
#   }
#   allLevel2plots[[refvalname]] <- p23
# }

