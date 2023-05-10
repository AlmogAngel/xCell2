library(tidyverse)
library(gridExtra)

# Score xCell2.0 signatures

scoreMixture <- function(bulk, signatures_collection){

  celltypes <- unique(unlist(lapply(names(signatures_collection), function(x) {strsplit(x, "#")[[1]][1]})))
  bulk_ranked <- singscore::rankGenes(bulk)

  scores.list <- lapply(celltypes, function(ctoi){
    signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), paste0(ctoi, "#"))]

    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    if (is.list(scores)) {
      signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(bulk_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
    }

    colnames(scores) <- names(signatures_ctoi)
    rownames(scores) <- colnames(bulk)
    t(scores)
  })
  names(scores.list) <- celltypes

  return(scores.list)

}

getCors <- function(sigs.out.list, truth, file){

  all_celltypes <- intersect(names(sigs.out.list), rownames(truth))

  all_celltypes_cor <- lapply(all_celltypes, function(ctoi){


    truth_ctoi <- truth[ctoi, ]
    truth_ctoi <- truth_ctoi[1, !is.na(truth_ctoi)]
    truth_ctoi <- truth_ctoi[1, truth_ctoi != ""]

    results <- sigs.out.list[[ctoi]]
    shared_samples <- intersect(colnames(truth_ctoi), colnames(results))


    cors <- apply(results, 1, function(x){
      cor(x[shared_samples], as.numeric(truth_ctoi[,shared_samples]), method = "spearman")
    })

    cors <- sort(cors, decreasing = T)
    cors

  })
  all_celltypes_cor <- unlist(all_celltypes_cor)
  names(all_celltypes_cor) <- paste0(file, "#", names(all_celltypes_cor))


  return(all_celltypes_cor)

}


# Load benchmarking truths and mixtures
truths_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/"
mix_dir <- "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/"

# Load celltype conversion
celltype_conversion_long <- read_tsv("/bigdata/almogangel/xCell2/dev_data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))


blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813")

all_ds_cors <- c()
for (file in blood_ds) {
  bulk <- as.matrix(read.table(paste0(mix_dir, file, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))

  # Load truth
  truth <- read.table(paste0(truths_dir, file, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)


  # Run xCell2.0
  sigs.out.list <- scoreMixture(bulk, signatures_collection)
  # Add mean score of signatures
  for (i in 1:length(sigs.out.list)) {
    mean_score <- colMeans(sigs.out.list[[i]])
    sigs.out.list[[i]] <- rbind(mean_score, sigs.out.list[[i]])
    rownames(sigs.out.list[[i]])[1] <- paste0(names(sigs.out.list[i]), "#", "mean_score")
  }
  cors.out <- getCors(sigs.out.list, truth, file)

  all_ds_cors <- c(all_ds_cors, cors.out)
}



# Analyze worst signatures
all_sigs_results <- tibble(id = names(all_ds_cors), cor = all_ds_cors, performance = "All sigs") %>%
  separate(id, into = c("dataset", "signature"), sep = "#", extra = "merge") %>%
  separate(signature, into = c("celltype", "sig_info"), sep = "#", remove = FALSE) %>%
  separate(sig_info, into = c("remove", "probs", "diff", "n_genes"), sep = "_", extra = "drop", remove = FALSE) %>%
  select(-c(sig_info, remove))

worst_sigs <- all_sigs_results %>%
  group_by(dataset, celltype) %>%
  slice_min(prop = 0.05, order_by = cor) %>%
  ungroup() %>%
  mutate(performance = paste0("Bottom 5% - ", dataset))

ds2use <- "BG_blood"

rbind(all_sigs_results, worst_sigs) %>%
  filter(dataset == ds2use & performance == "All sigs" | signature %in% worst_sigs[worst_sigs$dataset == ds2use,]$signature & performance != "All sigs") %>%
  ggplot(., aes(x=celltype,y=cor, fill=performance)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  scale_fill_manual(values=c("black", "red", "green", "purple", "blue")) +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))

ctoi <- "T cell"
rbind(all_sigs_results, worst_sigs) %>%
  filter(celltype == ctoi) %>%
  filter(dataset == ds2use & performance == "All sigs" | signature %in% worst_sigs[worst_sigs$dataset == ds2use,]$signature & performance != "All sigs") %>%
  rowwise() %>%
  mutate(dataset = if(performance == "All sigs") paste0(ds2use, " - All sigs") else dataset) %>%
  ggplot(., aes(x=dataset, y=cor, fill=performance)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  scale_fill_manual(values=c("black", "red", "green", "purple", "blue")) +
  labs(title = ctoi, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))

all_sigs_results %>%
  drop_na() %>%
  mutate(n_genes = factor(n_genes, levels = as.character(sort(as.numeric(unique(all_sigs_results$n_genes)), decreasing = T)))) %>%
  ggplot(., aes(x=celltype,y=cor, fill=n_genes)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  facet_wrap(~dataset, scales = "free_x") +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))

all_sigs_results %>%
  drop_na() %>%
  mutate(n_genes = factor(n_genes, levels = as.character(sort(as.numeric(unique(all_sigs_results$n_genes)), decreasing = T)))) %>%
  ggplot(., aes(x=dataset,y=cor, fill=n_genes)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


all_sigs_results %>%
  drop_na() %>%
  mutate(probs = factor(probs, levels = as.character(sort(as.numeric(unique(all_sigs_results$probs)), decreasing = T)))) %>%
  ggplot(., aes(x=celltype,y=cor, fill=probs)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  facet_wrap(~dataset, scales = "free_x") +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))

all_sigs_results %>%
  drop_na() %>%
  mutate(probs = factor(probs, levels = as.character(sort(as.numeric(unique(all_sigs_results$probs)), decreasing = T)))) %>%
  ggplot(., aes(x=dataset,y=cor, fill=probs)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


all_sigs_results %>%
  drop_na() %>%
  mutate(diff = factor(diff, levels = as.character(sort(as.numeric(unique(all_sigs_results$diff)), decreasing = T)))) %>%
  ggplot(., aes(x=celltype,y=cor, fill=diff)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  facet_wrap(~dataset, scales = "free_x") +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))

all_sigs_results %>%
  drop_na() %>%
  mutate(diff = factor(diff, levels = as.character(sort(as.numeric(unique(all_sigs_results$diff)), decreasing = T)))) %>%
  ggplot(., aes(x=dataset,y=cor, fill=diff)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


all_sigs_results %>%
  drop_na() %>%
  mutate(probs_diff = paste0(probs, "-", diff)) %>%
  ggplot(., aes(x=probs_diff, y=cor, fill=probs_diff)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  labs(title = ds2use, y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))


# Analyze signatures

cor_final <- tibble(id = names(all_ds_cors), cor = all_ds_cors, method = "xCell2.0 - before filtering") %>%
  separate(id, into = c("dataset", "signature"), sep = "#", extra = "merge")

# For top signature with simbu simulations
top_sigs_by_sim <- cor_final %>%
  filter(signature %in% unlist(top_sigs)) %>%
  mutate(method = "xCell2.0 - top 10% simulations")
cor_final <- rbind(top_sigs_by_sim, cor_final)

# For top signature with RF regulaztion
top_sigs_by_rf <- cor_final %>%
  filter(signature %in% unlist(top_sigs_rf)) %>%
  mutate(method = "xCell2.0 - top 10% RF")
cor_final <- rbind(top_sigs_by_rf, cor_final)

cor_final <- cor_final %>%
  separate(signature, into = c("celltype", "is_mean"), sep = "#", extra = "drop", remove = FALSE) %>%
  mutate(method = ifelse(is_mean == "mean_score", "xCell2.0 - mean", method))

# Merge cor_final with benchmark_bulk_validation.R results to compare signatures to the alternative method
xcell2sigs_cors <- cor_final %>%
  mutate(celltype = gsub("-", "_", celltype)) %>%
  mutate(celltype = gsub(" ", "_", celltype)) %>%
  mutate(celltype = gsub(",", "_", celltype)) %>%
  select(dataset, method, celltype, cor)

all_methods <- alternative_methods_results.tbl %>%
  select(-results) %>%
  rbind(xcell2sigs_cors)

celltypes2use <- intersect(unique(xcell2sigs_cors$celltype), unique(alternative_methods_results.tbl$celltype))


all_methods %>%
  filter(celltype %in% celltypes2use) %>%
  #mutate(passed_filter = factor(ifelse(signature %in% names(signatures_collection_filtered) , "yes", "no"), levels = c("yes", "no"))) %>%
  ggplot(., aes(x=dataset, y=cor, col=method,
                alpha = method %in% c("xCell2.0 - length 8", "xCell2.0 - length 100"),
                size = method %in% c("xCell2.0 - length 8", "xCell2.0 - length 100"))) +
  geom_jitter(position=position_jitter(width=0.2)) +
  scale_color_manual(values=c("blue", "red", "green", "orange", "black", "purple", "gray20", "#98F5FF", "#FF00FF", "#FFA500")) +
  scale_size_ordinal(range = c(3, 1), guide = FALSE) +
  scale_alpha_manual(values = c(1, 0.1), guide = FALSE) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~celltype, scales = "free_x") +
  labs(y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))














blood_ref_sigs_cors <- getSigsCor(datasets = cytometry_validation, signatures_collection = xcell2_sigs)


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

