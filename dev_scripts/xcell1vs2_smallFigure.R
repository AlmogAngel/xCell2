# Signature generation small figure -----

# Make signatures with BP and BG_blood
signatures

# Use another reference to check signatures
ref_tmp <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/references/kass_blood_ref.rds")

# Take the shared cell types
cts2use <- intersect(unique(labels$label), unique(ref_tmp$labels$label))

scores1 <- parallel::mclapply(cts2use, function(ct1){

  print(paste0("ct1-", ct1))

  ref_tmp_ct1 <- ref_tmp$ref[,ref_tmp$labels$label == ct1]
  if (class(ref_tmp_ct1)[1] == "numeric") {
    ref_tmp_ct1 <- cbind(ref_tmp_ct1, ref_tmp_ct1)
    colnames(ref_tmp_ct1) <- c("1", "2")
  }


  ref_tmp_ct1_ranked <- singscore::rankGenes(ref_tmp_ct1)


  scores2 <- sapply(cts2use, function(ct2){
    print(paste0("ct2-", ct2))

    signatures_ct2 <- signatures[startsWith(names(signatures), paste0(ct2, "#"))]

    scores <- colMedians(sapply(signatures_ct2, simplify = TRUE, function(sig){
      suppressWarnings(singscore::simpleScore(ref_tmp_ct1_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
    }))

    scores

  })

  scores2


}, mc.cores = 20)
names(scores1) <- cts2use


scores1.tbl <- enframe(scores1, name = "ref_celltype") %>%
  unnest(value) %>%
  unnest_longer(value, indices_to = "sig", values_to = "score") %>%
  separate(sig, into = "sig_celltype", extra = "drop", sep = "#", remove = FALSE)


# Assuming your data is in a variable called scores1.tbl
# 1. Aggregate data
aggregated_data <- scores1.tbl %>%
  group_by(ref_celltype, sig_celltype) %>%
  summarize(mean_score = median(score, na.rm = TRUE))

# 2. Reshape data to wide format
wide_data <- aggregated_data %>%
  pivot_wider(names_from = sig_celltype, values_from = mean_score, values_fill = list(mean_score = NA))
ref_celltypes <- wide_data$ref_celltype

# Remove the ref_celltype column to create a matrix
heatmap_matrix <- as.matrix(wide_data[-1])
rownames(heatmap_matrix) <- ref_celltypes

# 3. Create heatmap
colMain <- colorRampPalette(brewer.pal(9, "Blues"))(30)
colMain2 <- c(colMain[1], colMain[5], colMain[10], colMain[15], colMain[17], colMain[30])


tcells <- c("CD8-positive, alpha-beta memory T cell", "CD4-positive, alpha-beta T cell", "regulatory T cell")

colnames(heatmap_matrix) <- gsub("class switched ", "", colnames(heatmap_matrix))
rownames(heatmap_matrix) <- gsub("class switched ", "", rownames(heatmap_matrix))

bcells <- c("plasma cell", "memory B cell", "naive B cell")
pheatmap::pheatmap(heatmap_matrix[bcells, bcells], cluster_rows = F, cluster_cols = F,
         show_rownames = TRUE, show_colnames = TRUE, scale = "column", color = colMain2, fontsize = 16)





# Signature filterin small figure ------

# Make signatures with BP and GSE127813

ds_cors_list

# CD8-positive, alpha-beta T cell
best_sigs

x <- enframe(ds_cors_list, "val") %>%
  unnest_longer(value, indices_to = "sig", values_to = "score") %>%
  separate(sig, into = "sig_celltype", extra = "drop", sep = "#", remove = FALSE)
x$is_filt <- "No"

y <- x[x$sig %in% best_sigs,]
y$is_filt <- "Yes"

data <- rbind(x, y)

data$is_filt <- factor(data$is_filt, levels = c("No", "Yes"))

color_palette <- c("#FC4E07", "#00AFBB")

data %>%
  ggplot(., aes(x=val, y=score, fill=is_filt)) +
  geom_boxplot() +
  theme_minimal() +  # A minimal theme for a nicer look
  labs(title = "CD8+ T cells",x = "", y = "Score", fill = "Filtered?") +
  scale_fill_manual(values = as.character(color_palette)) +
  scale_y_continuous(limits = c(-0.5, 1)) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20), # Center, bold, and increase size of title
        axis.title.y = element_text(size = 15)) # Increase size of y-axis title
