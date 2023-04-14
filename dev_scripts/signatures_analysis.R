library(tidyverse)
library(gridExtra)

# Choose signatures to analyze
xcell2_sigs <- readRDS("/home/almogangel/xCell2_git/Data/benchmarking_data/xcell2ref_ts_main.rds")


# Load benchmarking truths and mixtures
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"
mix_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"
ds <- gsub(".tsv", "", list.files(truths_dir))


scoreMixtures <- function(ctoi, mixture_ranked, signatures_collection){
  signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), paste0(ctoi, "#"))]

  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })

  # In case some signatures contain genes that are all not in the mixtures
  if (is.list(scores)) {
    signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
  }
  colnames(scores) <- names(signatures_ctoi)
  rownames(scores) <- colnames(mixture_ranked)
  return(t(scores))
}

getSigsCor <- function(datasets, signatures_collection){

  all_celltypes <- unique(gsub(pattern = "#.*", "", names(signatures_collection)))
  all_ds <- sapply(datasets, function(x) NULL)

  for (file in datasets) {

    # Load mixture
    mix <- read.table(paste0(mix_dir, file, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

    # Load truth
    truth <- read.table(paste0(truths_dir, file, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

    samples2use <- intersect(colnames(truth), colnames(mix))
    mix <- mix[,samples2use]
    truth <- truth[,samples2use]

    mix_ranked <- singscore::rankGenes(mix)

    if (!all(colnames(truth) == colnames(mix_ranked))) {
      errorCondition(paste0("Error with dataset: ", file))
    }


    all_celltypes_cor <- sapply(all_celltypes, function(ctoi){

      scores_ctoi <- scoreMixtures(ctoi, mix_ranked, signatures_collection)
      truth_ctoi <- truth[ctoi, colnames(scores_ctoi)]
      truth_ctoi <- truth_ctoi[1, !is.na(truth_ctoi)]
      truth_ctoi <- truth_ctoi[1, truth_ctoi != ""]
      scores_ctoi <- scores_ctoi[,names(truth_ctoi)]

      if (!all(names(truth_ctoi) == colnames(scores_ctoi))) {
        errorCondition(paste0("Error with dataset: ", file))
      }

      truth_ctoi <- as.numeric(truth_ctoi)

      if (all(truth_ctoi == 0)) {
        NULL
      }else{
        apply(scores_ctoi, 1, function(x){cor(x, truth_ctoi, method = "spearman")})
      }


    })

    all_ds[[file]] <- all_celltypes_cor
  }

  all_ds_cors <- unlist(all_ds)

  cor_final <- tibble(id = names(all_ds_cors), cor = all_ds_cors) %>%
    separate(id, into = c("dataset", "celltype", "signature"), sep = "\\.", extra = "merge")


  return(cor_final)

}

getTopSigs <- function(cor_final, top = 0.1){
  top_sigs <- blood_ref_sigs_cors %>%
    group_by(celltype, signature) %>%
    summarise(median_cor = median(cor)) %>%
    filter(median_cor >= quantile(median_cor, 1-top)) %>%
    pull(signature)

  return(top_sigs)
}



# Filter by Grubb's test - all
grubbs.filtered <- xcell2_blood_ref@score_mat %>%
  filter(signature %in% filtered_sigs) %>%
  group_by(signature_ct, signature) %>%
  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  filter(grubbs_statistic >= quantile(grubbs_statistic, 0.9)) %>%
  pull(signature)

# Filter by Grubb's test - similarity
ct_similarity <- xcell2_blood_ref@score_mat %>%
  rowwise() %>%
  mutate(similarity = xcell2_blood_ref@correlationMatrix[signature_ct, sample_ct]) %>%
  select(signature_ct, sample_ct, similarity) %>%
  unique() %>%
  group_by(signature_ct) %>%
  mutate(similarity_level = ifelse(similarity >= quantile(similarity, 0.8), "high", "low")) %>%
  mutate(similarity_level = ifelse(signature_ct == sample_ct, "same", similarity_level))

grubbs.filtered.sim <- xcell2_blood_ref@score_mat %>%
  left_join(ct_similarity, by = c("signature_ct", "sample_ct")) %>%
  filter(similarity_level != "low") %>%
  group_by(signature_ct, signature) %>%
  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  filter(grubbs_statistic >= quantile(grubbs_statistic, 0.9)) %>%
  pull(signature)


# Get signatures correlations with the validation datasets
blood_ref_sigs_cors <- getSigsCor(datasets = validation_ds_blood, signatures_collection = xcell2_blood_ref@all_signatures)

# saveRDS(blood_ref_sigs_cors, "../xCell2.0/blood_ref_sigs_cors.rds")
cytometry_validation <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813")

blood_ref_sigs_cors %>%
  ungroup() %>%
  mutate(passed_filter = factor(ifelse(signature %in% grubbs.filtered, "yes", "no"), levels = c("yes", "no"))) %>%
  group_by(dataset, celltype, passed_filter) %>%
  summarise(median_cor = median(cor)) %>%
  group_by(dataset, celltype) %>%
  arrange(passed_filter, .by_group = T) %>%
  mutate(delta_median =  median_cor - lag(median_cor)) %>%
  drop_na() %>%
  select(dataset, celltype, delta_median) %>%
  group_by(celltype) %>%
  summarise(delta_median = median(delta_median)) %>%
  arrange(-delta_median)


blood_ref_sigs_cors %>%
  ungroup() %>%
  #left_join(grubbs, by = "signature") %>%
  filter(dataset %in% cytometry_validation) %>%
  mutate(passed_filter = factor(ifelse(signature %in% names(xcell2_blood_ref@filtered_signatures) , "yes", "no"), levels = c("yes", "no"))) %>%
  #filter(dataset == "BG_blood") %>%
  ggplot(., aes(x=grubbs_statistic, y=cor, fill= dataset)) +
  geom_point(aes(col=passed_filter)) +
  stat_smooth(method = "lm") +
  scale_color_manual(values=c("#66CD00", "#CD3333")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.2)) +
  facet_wrap(~celltype, scales = "free_x")


no_filt <- blood_ref_sigs_cors %>%
  ungroup() %>%
  filter(dataset %in% cytometry_validation) %>%
  mutate(filter_type = "All")

grubbs_filt <- blood_ref_sigs_cors %>%
  ungroup() %>%
  filter(dataset %in% cytometry_validation & signature %in% grubbs.filtered) %>%
  mutate(filter_type = "Simulations+Grubbs")

simu_filt <- blood_ref_sigs_cors %>%
  ungroup() %>%
  filter(dataset %in% cytometry_validation & signature %in% filtered_sigs) %>%
  mutate(filter_type = "Simulations")

blood_ref_sigs_cors %>%
  mutate(passed_filter = factor(ifelse(signature %in% names(xcell2_blood_ref@filtered_signatures) , "yes", "no"), levels = c("yes", "no"))) %>%
  ggplot(., aes(x=dataset, y=cor)) +
  geom_boxplot(aes(fill=passed_filter)) +
  scale_fill_manual(values=c("#66CD00", "#CD3333", "blue")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~celltype, scales = "free_x") +
  labs(y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


rbind(no_filt, rbind(grubbs_filt, simu_filt)) %>%
  mutate(filter_type = factor(filter_type, levels = c("All", "Simulations", "Simulations+Grubbs"))) %>%
  filter(celltype %in% unique(simu_filt$celltype)[37:40]) %>%
  ggplot(., aes(x=dataset, y=cor)) +
  geom_boxplot(aes(fill=filter_type)) +
  scale_fill_manual(values=c("#66CD00", "#CD3333", "blue")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~celltype, scales = "free_x") +
  labs(y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


rbind(no_filt, rbind(grubbs_filt, simu_filt)) %>%
  mutate(filter_type = factor(filter_type, levels = c("All", "Simulations", "Simulations+Grubbs"))) %>%
  filter(celltype %in% unique(simu_filt$celltype)) %>%
  ggplot(., aes(x=filter_type, y=cor)) +
  geom_boxplot(aes(fill=filter_type)) +
  stat_summary(fun.y="mean")+
  scale_fill_manual(values=c("#66CD00", "#CD3333", "blue")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~dataset, scales = "free_x") +
  labs(y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


# sig_score.df <- xcell2_blood_ref@score_mat  %>%
#   filter(signature_ct == "Transitional memory T-helpers") %>%
#   dplyr::select(signature, sample_ct, score) %>%
#   pivot_wider(names_from = sample_ct, values_from = score) %>%
#   column_to_rownames(var = "signature")
# sigs_to_use <- rownames(sig_score.df)
# celltype_order <- names(sort(xcell2_blood_ref@correlationMatrix["Transitional memory T-helpers",], decreasing = TRUE))
# celltype_order <- celltype_order[celltype_order %in% names(sig_score.df)]
# sigs_to_use <- sigs_to_use[sigs_to_use %in% names(xcell2_blood_ref@filtered_signatures)]
# sig_score.df <- sig_score.df[sigs_to_use, celltype_order]
# sig_score.df <- sig_score.df[ ,colSums(is.na(sig_score.df)) == 0]
# pheatmap::pheatmap(sig_score.df, cluster_rows=F, cluster_cols=F, scale = "row", col= RColorBrewer::brewer.pal(11, "RdBu"))

