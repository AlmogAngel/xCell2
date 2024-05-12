library(tidyverse)


# A - Signatures generation and ontology integration --------------------------------------------------

# GSE107572_bp
truth_mat <- cyto.vals$truth$blood$GSE107572


###
# Signatures generation heatmap
###

# Required from xCell2Train():
mix
cor_mat
gep_mat
signatures


ctoi <- "CD8-positive, alpha-beta T cell"
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]

length(signatures_ctoi)
# [1] 730
length(unique(unlist(signatures_ctoi)))
# [1] 390


# Get signatures scores with reference
gep_mat_ranked <- singscore::rankGenes(gep_mat[,names(sort(cor_mat[ctoi,], decreasing = T))])
gep_mat_ranked_scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(gep_mat_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
})
rownames(gep_mat_ranked_scores) <- colnames(gep_mat_ranked)


# Heatmap
color_palette <- colorRampPalette(c("darkred" ,"red", "white", "lightblue1", "darkblue"))(50)

ct2use <- c(ctoi, "central memory CD8-positive, alpha-beta T cell",
            "effector memory CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell",
            "central memory CD4-positive, alpha-beta T cell", "effector memory CD4-positive, alpha-beta T cell",
            "regulatory T cell", "eosinophil", "naive B cell", "plasma cell", "neutrophil", "monocyte",  "epithelial cell", "fibroblast")

data <- gep_mat_ranked_scores[ct2use, names(sort(gep_mat_ranked_scores[1,]))]
rownames(data) <- c("CD8+ T", "CM CD8+ T", "EM CD8+ T", "CD4+ T", "CM CD4+ T", "EM CD4+ T",
                    "T-reg", "Eosinophil", "Naive B", "Plasma cell", "Neutrophil", "Monocyte",  "Epithelial", "Fibroblast")

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

# Plot results
data <- tibble(cors = c(c_onto, c_noOnto), label = c(rep("Ontology", length(c_onto)), rep("No Ontology", length(c_noOnto))))

ggpubr::compare_means(cors~label, data, method = "wilcox.test")

data %>%
  ggplot(., aes(x=label, y= cors, fill=label)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "", x = "", y = "Spearman Rho with ground truth", fill = "") +
  scale_fill_manual(values = c("lightblue1", "darkblue")) +
  scale_y_continuous(limits = c(-0.8, 1.001), breaks = seq(-0.8, 1, 0.2)) +
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

data_mix <- rbind(data, data2)
data_mix$data_type <- "Input Dataset"
data_mix$dataset <- "GSE107572"

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
color_palette <- c("gray", "darkblue")
data_final %>%
  ggplot(., aes(x=dataset, y=rho, fill=filtered)) +
  geom_boxplot(width = 0.8) +
  facet_grid(~data_type, scales = "free_x", space = "free") +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Spearman Rho", fill = "Signatures") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(-0.52, 1), breaks = seq(-0.5, 1, 0.25)) +  # Set y-axis limits
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



# B --------------------------------------------------


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

xcell2.sigs <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.filteredSigsCors.newEssentialv2valtype.3aprNewVal.rds")

xcell2.sigs <- xcell2.sigs %>%
  bind_rows() %>%
  filter(ref == "bp") %>%
  mutate(label = ifelse(is_filtered == "yes", "xCell 2.0 (filtered)", "xCell 2.0 (all)")) %>%
  select(-c(ref, is_filtered))


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

color_palette <- c("#EE3A8C", "gray", "darkblue")


results %>%
  group_by(val, celltype, label) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=val, y=cor, fill=label)) +
  geom_boxplot() +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Mean Spearman Rho", fill = "") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(-0, 1.001), breaks = seq(0, 1, 0.1)) +  # Set y-axis limits
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




# Archive code  --------------------------------------------------

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
