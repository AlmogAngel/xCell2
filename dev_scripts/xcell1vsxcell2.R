#####################################################
# xCell2 vs. xCell1 signatures comparison
#####################################################


library(tidyverse)
library(xCell2)

refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")

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


# Get xCell2's signatures data (before/after filtering) ---------------

xcell2.sigs <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.filteredSigsCors.newEssentialv2.13mar.rds")

xcell2.sigs <- xcell2.sigs %>%
  bind_rows() %>%
  filter(ref == "bp") %>%
  mutate(label = ifelse(is_filtered == "yes", "xCell2 (filtered)", "xCell2")) %>%
  select(-c(ref, is_filtered))



# Get correlations for first xCell version signatures ----------------------------------------

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

# Get xCell2's signatures data (model predictions) ---------------
all_cors # from fig1.R

xcell2.sigs.model <- all_cors %>%
  filter(method == "xCell2" & ref == "bp" & celltype %in%  unique(xcell.sigs.data$celltype)) %>%
  mutate(sig_name = "model",
         label = "xCell2 (model)") %>%
  select(val, celltype, sig_name, cor, label)


# Plot results --------------------

#results <- rbind(xcell.sigs.data, xcell2.sigs, xcell2.sigs.model)
#results <- rbind(xcell.sigs.data, xcell2.sigs)

# Identify common val and celltype combinations
# common_combinations <- reduce(list(xcell.sigs.data, xcell2.sigs, xcell2.sigs.model),
#                               ~inner_join(.x, .y, by = c("val", "celltype"))) %>%
#   select(val, celltype) %>%
#   distinct()

common_combinations <- reduce(list(xcell.sigs.data, xcell2.sigs),
                              ~inner_join(.x, .y, by = c("val", "celltype"))) %>%
  select(val, celltype) %>%
  distinct()

# Filter each tibble for common combinations and combine them
# results <- list(xcell.sigs.data, xcell2.sigs, xcell2.sigs.model) %>%
#   map(~inner_join(.x, common_combinations, by = c("val", "celltype"))) %>%
#   bind_rows()

results <- list(xcell.sigs.data, xcell2.sigs) %>%
  map(~inner_join(.x, common_combinations, by = c("val", "celltype"))) %>%
  bind_rows()


color_palette <- c("#00AFBB", "#E7B800", "#FC4E07", "#828282")


results %>%
  group_by(val, celltype, label) %>%
  summarise(cor = combineRhos2(cor)) %>%
  ggplot(., aes(x=label, y=cor, fill=label)) +
  geom_boxplot() +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "",x = "", y = "",fill = "") +
  scale_fill_manual(values = as.character(color_palette)) +
  scale_y_continuous(limits = c(-1, 1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for x-axis labels if needed
        legend.position = "bottom")



results %>%
  group_by(val, celltype, label) %>%
  summarise(cor = combineRhos2(cor)) %>%
  group_by(val, label) %>%
  summarise(cor = combineRhos2(cor)) %>%
  ggplot(., aes(x=label, y=cor, fill=label)) +
  geom_boxplot() +
  geom_jitter(alpha=0.3) +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Mean Spearman Rho", fill = "") +
  scale_fill_manual(values = as.character(color_palette), guide = FALSE) +  # Remove fill legend
  scale_y_continuous(limits = c(0.2, 1.001), breaks = seq(0.2, 1, 0.1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),  # Make x-axis labels bold and bigger, adjust as needed
        legend.position = "bottom",  # Adjust legend position if needed
        legend.title = element_blank(),  # Optionally, hide the legend title if it's still showing
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        axis.line = element_line(colour = "black"))


# Micro (filtering) -------------------------
data <- results %>%
  group_by(val, celltype, label) %>%
  summarise(cor = combineRhos2(cor)) %>%
  group_by(val, label)


wide_data <- data %>%
  pivot_wider(names_from = label, values_from = cor)

matching_cases <- wide_data %>%
  filter(`xCell2 (model)` > `xCell2 (filtered)`,
         `xCell2 (filtered)` > xCell2,
         xCell2 > xCell)

# GSE120444 - CD8-positive, alpha-beta T cell

ctoi <- "CD8-positive, alpha-beta T cell"
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
mix_ranked_tmp <- singscore::rankGenes(mix)
scores_tmp <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked_tmp, upSet = sig, centerScore = FALSE)$TotalScore)
})
truth_mat <- cyto.vals$truth$other$GSE120444
fracs <- truth_mat[ctoi,colnames(mix_ranked_tmp)]

c <- apply(scores_tmp, 2, function(x){
  cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
})
c <- sort(c, decreasing = TRUE)


data <- tibble(sig = names(c), rho = c, filtered = "All")
filt_sigs <- names(out$filt_sigs[startsWith(names(out$filt_sigs), paste0(ctoi, "#"))])
data2 <- tibble(sig = filt_sigs, rho = c[filt_sigs], filtered = "Filtered")

data_mix <- rbind(data, data2)
data_mix$data_type <- "Mixture Dataset"
data_mix$dataset <- "GSE120444"

data_mix<- data_mix %>%
  select(dataset, rho, sig, filtered, data_type)

# filtering (from the function)

data <- enframe(ds_cors_list, name = "dataset") %>%
  unnest_longer(value, values_to = "rho", indices_to = "sig") %>%
  mutate(filtered = "All")

data2 <- data %>%
  filter(sig %in% names(out$filt_sigs)) %>%
  mutate(filtered = "Filtered")

data_val <- rbind(data, data2)
data_val$data_type <- "Filtering Datasets"

data_final <- rbind(data_mix, data_val)

color_palette <- c("#00AFBB", "#E7B800")


data_final %>%
  ggplot(., aes(x=dataset, y=rho, fill=filtered)) +
  geom_boxplot(width = 0.8) +
  facet_grid(~data_type, scales = "free_x", space = "free") +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "", x = "", y = "Spearman Rho", fill = "") +
  scale_fill_manual(values = as.character(color_palette)) +  # Remove fill legend
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.25)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),  # Make x-axis labels bold and bigger, adjust as needed
        title = element_text(face = "bold", size = 14),  # Optionally, hide the legend title if it's still showing
        legend.position = "right",  # Adjust legend position if needed
        legend.text = element_text(face = "bold", size = 14),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_rect(fill = '#F0F8FF', color = NA),  # Remove panel background
        strip.text = element_text(face="bold", size=14),
        axis.line = element_line(colour = "black"))

library(patchwork)

# Combine the plots side by side
combined_plot <- p1 + p2

# Display the combined plot
print(combined_plot)

# Micro (model) -----------------

set.seed(123)
noise_level <- 0.5
ctoi_sim_list <- lapply(1:10, function(i){

  # Sample 3-10 control cell types (if possible)
  if (length(controls) > 3) {
    controls2use <- sample(controls, sample(3:min(10, length(controls)), 1), replace = FALSE)
  }else{
    controls2use <- controls
  }

  # Generate controls matrix
  if (length(controls2use) > 1) {
    controls_mat <- gep_mat_linear[,controls2use]
    numbers <- runif(ncol(controls_mat)) # Get random numbers for fractions
    controls_mat_frac <- sapply(1-sim_fracs, function(s){
      fracs <- numbers / sum(numbers) * s
      rowSums(controls_mat %*% diag(fracs))
    })
  }else{
    controls_mat <- matrix(rep(gep_mat_linear[,controls2use], length(sim_fracs)), byrow = FALSE, ncol = length(sim_fracs), dimnames = list(rownames(gep_mat_linear), sim_fracs))
    controls_mat_frac <- controls_mat %*% diag(sim_fracs)
  }


  # Combine fractions
  simulation <- ctoi_mat_frac + controls_mat_frac
  colnames(simulation) <- paste0("mix", "%%", sim_fracs)


  # Add noise
  if (!is.null(noise_level)) {
    eps_sigma <- sd(as.vector(ref)) * noise_level
    epsilon_gaussian <- array(rnorm(prod(dim(simulation)), 0, eps_sigma), dim(simulation))
    simulation_noised <- 2^(log2(simulation+1) + epsilon_gaussian)
    return(simulation_noised)
  }else{
    simulation
  }

})

ctoi <- "CD8-positive, alpha-beta T cell"
signatures_filtered <- out$filt_sigs
signatures_ctoi <- signatures_filtered[gsub("#.*", "", names(signatures_filtered)) %in% ctoi]


xx <- lapply(ctoi_sim_list, function(sim){
  sim_ranked <- singscore::rankGenes(sim)
  colnames(sim_ranked) <- make.unique(colnames(sim_ranked), sep = "%")
  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(sim_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })

  apply(scores, 2, function(x){
    cor(x, sim_fracs, method = "spearman")
  })

})
xx










# Waterfall?

cors_data <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.filteredSigsCors.22jan.rds")

data <- bind_rows(cors_data)


data2 <- data %>%
  group_by(ref, val, celltype, is_filtered) %>%
  summarise(cor = mean(cor)) %>%
  spread(key = is_filtered, value = cor) %>%
  mutate(delta_cor = `yes` - `no`) %>%
  select(ref, val, celltype, delta_cor)


data_sorted <- data2 %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  mutate(x_axis = row_number())

# Create a waterfall plot
ggplot(data_sorted, aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "white", midpoint = 0) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Delta Cor", title = "Waterfall Plot of Delta Cor")


data2 %>%
  filter(delta_cor < 0) %>%
  group_by(val) %>%
  summarise(n = n())



library(RColorBrewer)
color_palette <- Polychrome::createPalette(22,  c("#ff0000", "#00ff00", "#0000ff"))

results %>%
  ggplot(., aes(x=label, y=cor, fill=label)) +
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

results$label <- factor(results$label, levels= c("xCell", "xCell2", "xCell2 (filtered)"))

stat.test <- results %>%
  ungroup() %>%
  t_test(cor ~ label)
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
