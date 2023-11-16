#####################################################
# xCell2 present improved signatures generation
#####################################################

# Goal: to generate boxplot that show that xCell2's signatures are better than xCell given the same reference

library(tidyverse)

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")


# Load xCell's signatures --------
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


# Get xCell2's signatures from various validations ---------------



refList <- list(rna_seq = c(mixed = "bp"))

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
  mutate(method = "xCell2", .before = everything()) %>%
  filter(ref_name == "bp")


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
sigs_suffix <- "xcell1v2"
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
spillover <- TRUE
nSims <- 20


xcell2.sigs <- lapply(1:nrow(vals.refs.res), function(i){

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


  shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
  ref <- shared_cleaned_genes$ref
  mix <- shared_cleaned_genes$mix

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
  sigs <-  xCell2::xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                               sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed,
                               nCores = cores2use, simMethod = sim_method, mix = mix2use, sim_fracs = simFracs, filtLevel = fLvel, ct_sims = nSims)

  saveRDS(sigs, file = paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds"))
  sigs
})

xcell2.sigs.filt <- lapply(1:nrow(vals.refs.res), function(i){

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


  shared_cleaned_genes <-  xCell2::xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = ProteinCoding, gene_groups = genesGroups)
  ref <- shared_cleaned_genes$ref
  mix <- shared_cleaned_genes$mix

  # Load signatures?
  if (TRUE) {
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

  saveRDS(sigs@signatures, file = paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigsFiltered_", sigs_suffix, ".rds"))
  sigs
})





# TODO ----------------------------------------

xcell2.sigs.files <- list.files("/bigdata/almogangel/xCell2_data/dev_data/sigs", pattern = "*_sigs_xcell1v2", full.names = TRUE)

cor.res <- lapply(xcell2.sigs.files, function(val.file){

  val <- gsub("_bp_.*", "", basename(val.file))
  print(val)

  # Load truth table
  truth_mat <- cyto.vals$truth$blood[[val]]
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$tumor[[val]]
  }
  if (is.null(truth_mat)) {
    truth_mat <- cyto.vals$truth$other[[val]]
  }
  truth_cts <- rownames(truth_mat)

  # Load xCell2 signatures
  xcell2.sigs <- readRDS(val.file)
  xcell2.cts <- unique(gsub("#.*", "", names(xcell2.sigs)))
  all(xcell.cts %in% xcell2.cts & xcell2.cts %in% xcell.cts)

  mix <- cyto.vals$mixtures$blood[[val]]
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$tumor[[val]]
  }
  if (is.null(mix)) {
    mix <- cyto.vals$mixtures$other[[val]]
  }
  mix_ranked <- singscore::rankGenes(mix)

  ct2use <- intersect(intersect(xcell.cts, xcell2.cts), truth_cts)
  cor.res <- lapply(ct2use, function(ct){

    truth <- truth_mat[ct,]

    ct_sigs <- xcell.sigs[[ct]]
    xcell.scores <- sapply(ct_sigs, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })
    rownames(xcell.scores) <- colnames(mix_ranked)
    samples <- intersect(colnames(truth_mat), rownames(xcell.scores))
    xcell.scores <- xcell.scores[samples,]
    xcell.sigs.cors <- apply(xcell.scores, 2, function(x){
      cor(truth[samples], x, method = "spearman", use = "pairwise.complete.obs")
    })

    ct_sigs <- xcell2.sigs[startsWith(names(xcell2.sigs), ct)]
    xcell2.scores <- sapply(ct_sigs, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    })
    rownames(xcell2.scores) <- colnames(mix_ranked)
    samples <- intersect(colnames(truth_mat), rownames(xcell2.scores))
    xcell2.scores <- xcell2.scores[samples,]
    xcell2.sigs.cors <- apply(xcell2.scores, 2, function(x){
      cor(truth[samples], x, method = "spearman", use = "pairwise.complete.obs")
    })

    tibble(validation = val, celltype = ct, sig = c(names(xcell.sigs.cors), names(xcell2.sigs.cors)), cor = c(xcell.sigs.cors, xcell2.sigs.cors),
           method = c(rep("xCell", length(xcell.sigs.cors)), rep("xCell2", length(xcell2.sigs.cors))))

  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# saveRDS(cor.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell1vs2.sig.cors.rds")
cor.res <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell1vs2.sig.cors.rds")


# Add filtered signatures
xcell2.sigs.filt.files <- list.files("/bigdata/almogangel/xCell2_data/dev_data/sigs", pattern = "*_sigsFiltered_*", full.names = TRUE)

sigs.filt <- lapply(xcell2.sigs.filt.files, function(file){
  sig.in <- readRDS(file)
  names(sig.in)
})
sigs.filt <- unique(unlist(sigs.filt))

cor.res <- cor.res %>%
  filter(sig %in% sigs.filt) %>%
  mutate(method = "xCell2 - filtered") %>%
  rbind(., cor.res)

library(ggsignif)

cor.res %>%
  ggplot(., aes(x=validation, y=cor, fill=method)) +
  geom_boxplot() +
  coord_flip() +  # Flips the coordinates
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "",x = "", y = "",fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for x-axis labels if needed
    legend.position = "bottom")

library(RColorBrewer)
color_palette <- Polychrome::createPalette(22,  c("#ff0000", "#00ff00", "#0000ff"))

cor.res %>%
  group_by(validation, celltype, method) %>%
  summarise(cor = mean(cor)) %>%
  group_by(validation, method) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=method, y=cor, fill=method)) +
  geom_boxplot() +
  geom_jitter(position = position_dodge(width = 0.8), aes(color = validation), size = 3) +  # Adjusted geom_jitter
  scale_fill_brewer(palette = "Set1") +  # Using a color palette from RColorBrewer
  scale_color_manual(values = as.character(color_palette)) +
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "",x = "", y = "",fill = "", color = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for x-axis labels if needed
        legend.position = "bottom") +
  guides(fill = FALSE)

library(ggpubr)
library(rstatix)

cor.res$method <- factor(cor.res$method, levels= c("xCell", "xCell2", "xCell2 - filtered"))

stat.test <- cor.res %>%
  group_by(method, validation, celltype) %>%
  summarise(cor = mean(cor)) %>%
  ungroup() %>%
  t_test(cor ~ method)
stat.test <- stat.test %>% add_xy_position(x = "method")

ggboxplot(cor.res, x = "method", y = "cor", fill = "method",
                   palette = c("#00AFBB", "#E7B800", "#FC4E07"), xlab = "", ylab = "Spearman rho") +
  theme(legend.position = "none") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01,
    y.position = c(1.2, 1.4, 1.2), bracket.shorten = 0.05)



top_ct <- cor.res %>% select(validation,celltype) %>% unique() %>% group_by(celltype) %>% summarise(n=n()) %>% arrange(-n)
top_ct <- top_ct[1:6,]$celltype


cor.res %>%
  #filter(celltype %in% top_ct) %>%
  ggplot(., aes(x=method, y=cor, fill=method)) +
  geom_boxplot() +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "",x = "", y = "",fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for x-axis labels if needed
        legend.position = "bottom") +
  facet_wrap(~celltype)  # Create separate plots for each celltype





df <- lapply(celltypes, function(ct){

  truth <- truth_mat[ct,]
  res <- res_mat[ct,]

  samples <- intersect(names(res), names(truth))

  tibble(celltype = ct, truth = truth[samples], prediction = res[samples])

}) %>%
  bind_rows()

cor_results <- df %>%
  group_by(celltype) %>%
  summarize(
    cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
    p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
    n = sum(!is.na(truth) & !is.na(prediction))
  ) %>%
  mutate(cor = ifelse(is.na(cor), 0, cor)) %>%
  mutate(label = paste(" rho: ", round(cor, 3), "\n p: ", sprintf("%.2e", p_value), "\n n: ", n, sep = ""))
cor_results





fracs <- truth_mat["monocyte",colnames(mix_ranked)]
c <- apply(scores, 2, function(x){
  cor(x, fracs, method = "spearman")
})
