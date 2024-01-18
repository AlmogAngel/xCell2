# xCell2Analysis
res_mat <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = TRUE, spillover = FALSE, spillover_alpha = 0.5)
#res_mat <- res_mat[vals.refs.res[1,]$shared_celltypes[[1]], ]
#res_mat


truth_mat <- cyto.vals$truth$blood$GSE65133


celltypes <- intersect(rownames(res_mat), rownames(truth_mat))
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

markers <- read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/lm22_markers.txt", sep = "\t", header = T)
sig <- markers[markers$label == "monocyte",]$marker

sim_ranked <- singscore::rankGenes(mix)
sig <- c("FOXP3", "IL2RA", "MS4A6A", "PTPRCAP", "TIMD4")
# sig <- names(gene_passed)
r <- singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
names(r) <- colnames(sim_ranked)
t <- truth_mat["regulatory T cell", names(r)]
cor(t, r, method = "spearman", use = "pairwise.complete.obs")


# signatures <- sigs@signatures
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#ctoi <- rownames(truth_mat)[9]
# signatures_ctoi <- type_sigs
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
# mix_ranked <- singscore::rankGenes(samples_sim_mat)
# rownames(mix_ranked) <- rownames(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})

scores_tmp <- apply(scores_tmp, 2, function(sample){
  normalize(sample)
})

fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]
# fracs=as.numeric(gsub(".*mix%%", "", colnames(mix_ranked)))
c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
sort(c)


sort(table(unlist(signatures[names(sort(c, decreasing = T)[1:10])])), decreasing = T)





dim(sim_list[[ctoi]])




truth_mat <- cyto.vals$truth$blood$BG_blood
celltypes <- unique(gsub("#.*", "", names(signatures)))
all_ctoi_res <- lapply(celltypes, function(ctoi){
  if (!ctoi %in% rownames(truth_mat)) {
    return(NA)
  }
  print(ctoi)
  signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
  mix_ranked <- singscore::rankGenes(mix)
  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
  })
  fracs <- truth_mat[ctoi,colnames(mix_ranked)]
  c <- apply(scores, 2, function(x){
    cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
  })
  return(c)
})
names(all_ctoi_res) <- celltypes


# Check signatures correlations with all other validations
mixs <- cyto.vals$mixtures$blood
truths <- cyto.vals$truth$blood
vals <- names(truths)
vals <- vals[vals %in% c("GSE107011", "GSE107572", "GSE127813", "GSE53655")]
celltypes <- unique(gsub("#.*", "", names(signatures)))

all_res <- lapply(vals, function(val){

  mix <- mixs[[val]]
  truth_mat <- truths[[val]]

  ct2use <- intersect(celltypes, rownames(truth_mat))
  samples2use <- intersect(colnames(mix), colnames(truth_mat))
  mix <- mix[,samples2use]
  truth_mat <- truth_mat[,samples2use]
  mix_ranked <- singscore::rankGenes(mix)

  all_ctoi_res <- lapply(ct2use, function(ctoi){

    signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
    signatures_ctoi[unlist(lapply(signatures_ctoi, function(sig){all(sig %in% genes2use)}))]
    signatures_ctoi <- signatures_ctoi[which(lengths(signatures_ctoi) > 3)]
    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    fracs <- truth_mat[ctoi,]
    c <- apply(scores, 2, function(x){
      cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
    })
    return(c)
  })
  names(all_ctoi_res) <- ct2use

  return(  list(val = all_ctoi_res))

})
names(all_res) <- vals

top_sigs <- enframe(all_res) %>%
  unnest(value) %>%
  unnest(value)  %>%
  unnest_longer(value) %>%
  separate(value_id, into = c("celltype"), sep = "#", remove = FALSE) %>%
  group_by(name, celltype) %>%
  top_frac(n = 0.2, wt = value) %>%
  summarise(top_sigs = list(value_id))

top_sigs_cts <- unique(top_sigs$celltype)


xx=enframe(all_ctoi_res, name = "celltype") %>%
  filter(celltype %in% top_sigs_cts) %>%
  unnest_longer(value, values_to = "cor", indices_to = "signature") %>%
  mutate(sig_status = "before")

yy=xx %>%
  filter(signature %in% unique(unlist(top_sigs$top_sigs))) %>%
  mutate(sig_status = "after")

z=rbind(xx,yy)
z %>%
  group_by(celltype, sig_status) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=sig_status, y= cor, fill=sig_status)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
# res_mat=readRDS("/bigdata/almogangel/xCell2_data/dev_data/sigs/GSE20300_bp_resMat_8dec.rds")
# res_mat=readRDS("/bigdata/almogangel/xCell2_data/dev_data/sigs/BG_blood_lm22_resMatRefMulti_8dec.rds")
getCors <- function(res, truth = truth_mat){

  celltypes <- intersect(rownames(res), rownames(truth))

  df <- lapply(celltypes, function(ct){

    truth <- truth[ct,]
    res <- res[ct,]

    samples <- intersect(names(res), names(truth))

    tibble(celltype = ct, truth = truth[samples], prediction = res[samples])

  }) %>%
    bind_rows()

  df %>%
    group_by(celltype) %>%
    summarize(
      cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
      p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
      n = sum(!is.na(truth) & !is.na(prediction))
    ) %>%
    mutate(cor = ifelse(is.na(cor), 0, cor))

}

truth_mat <- cyto.vals$truth$blood$SDY420
getCors(res = res_mat, truth = truth_mat)







# Scoring mixture
param <- BiocParallel::MulticoreParam(workers = 30)
celltypes <- unique(labels$label)
all_ctoi_res <- BiocParallel::bplapply(celltypes, function(ctoi){

  signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
  mix_ranked <- singscore::rankGenes(mix)

  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
  })
  rownames(scores) <- colnames(mix)

  return(scores)
}, BPPARAM = param)
names(all_ctoi_res) <- celltypes

truth_mat <- cyto.vals$truth$blood$BG_blood

for(ctoi in celltypes){

  if (!ctoi %in% rownames(truth_mat)) {
    next
  }

  ctoi_scores <- all_ctoi_res[[ctoi]]
  means_sorted <- sort(colMeans(ctoi_scores))

  fracs <- truth_mat[ctoi,rownames(ctoi_scores)]
  cors <- apply(ctoi_scores, 2, function(x){
    cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
  })
  cors <- cors[names(means_sorted)]

  print(paste0(ctoi, ":", cor(means_sorted, cors)))

}
