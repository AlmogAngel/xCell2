library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

# From fig1
all_cors



data <- all_cors

data <- data %>%
  filter(!method %in% c("quanTIseq"))

data <- data %>%
  mutate(refMethod = paste0(ref, "-", method)) %>%
  group_by(refMethod, val) %>%
  summarise(cors_list = list(cor),
            n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
  rowwise() %>%
  mutate(cor = combineRhos(rhos = cors_list, sample_sizes = log(n_ct_samples), use_median = FALSE, summarize = TRUE)) %>%
  dplyr::select(refMethod, val, cor)


# Reshape the data to a wide format
data_wide <- data %>%
  pivot_wider(names_from = refMethod, values_from = cor)

# Ensure the data is numeric
vals <- pull(data_wide[,1])
data_wide <- data_wide %>% dplyr::select(-val)  # Remove non-numeric column 'val'

# Convert to matrix
data_matrix <- as.matrix(data_wide)
rownames(data_matrix) <- vals

column_annotation <- data.frame(ref = gsub("-.*", "", colnames(data_matrix)))
rownames(column_annotation) <- colnames(data_matrix)

gap_indices <- which(diff(as.numeric(factor(gsub("-.*", "", colnames(data_matrix))))) != 0)

# Define the color gradients for both halves
colors_pos <- colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(50)  # From dark red to green with yellow in the middle
colors_neg <- colorRampPalette(c("black", "darkred"))(50)  # From black to dark red

# Combine the two color sets ensuring that dark red is at the center
colors <- c(colors_neg, colors_pos[length(colors_pos):1])


#data_matrix[data_matrix<0] <- 0


# Generate the heatmap
# pheatmap::pheatmap(t(data_matrix),
#          annotation_row = column_annotation,
#          color = colors_pos,
#          breaks = seq(0, 1, 0.02),
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          labels_row = gsub(".*-", "", colnames(data_matrix)),
#          gaps_row = gap_indices,
#          legend_breaks = seq(-1, 1, by = 0.2))
#


#

prefixes <- sapply(strsplit(colnames(data_matrix), "-"), `[`, 1)
group_factor <- factor(prefixes, levels = unique(prefixes))
data_matrix <- t(data_matrix)

# Order heatmap
bp_vals <- all_cors %>% filter(method == "xCell2" & ref == "bp") %>% pull(val) %>% unique()
lm22_vals <- all_cors %>% filter(method == "xCell2" & ref == "lm22") %>% pull(val) %>% unique()
kassblood_vals <- all_cors %>% filter(method == "xCell2" & ref == "kass_blood") %>% pull(val) %>% unique()
tsblood_vals <- all_cors %>% filter(method == "xCell2" & ref == "ts_blood") %>% pull(val) %>% unique()
kasstumor_vals <- all_cors %>% filter(method == "xCell2" & ref == "kass_tumor") %>% pull(val) %>% unique()
pancancer_vals <- all_cors %>% filter(method == "xCell2" & ref == "sc_pan_cancer") %>% pull(val) %>% unique()
only_general <- c(bp_vals, lm22_vals)
only_general <- unique(only_general[!only_general %in% c(kassblood_vals, tsblood_vals, kasstumor_vals,pancancer_vals)])
only_blood <- c(kassblood_vals, tsblood_vals)
only_blood <- unique(only_blood[!only_blood %in% c(only_general, kasstumor_vals,pancancer_vals)])
only_tumor <- c(kasstumor_vals, pancancer_vals)
only_tumor <- unique(only_tumor[!only_tumor %in% c(only_general, only_blood)])
vals_ordered <- c(only_general, only_blood, only_tumor)


methods_orders <- data_combined %>% filter(ref == "bp") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "bp")]
suffixes <- sub("bp-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_bp <- rnames[order(order_index)]

methods_orders <- data_combined %>% filter(ref == "lm22") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "lm22")]
suffixes <- sub("lm22-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_lm22 <- rnames[order(order_index)]

methods_orders <- data_combined %>% filter(ref == "kass_blood") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "kass_blood")]
suffixes <- sub("kass_blood-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_kass_blood <- rnames[order(order_index)]

methods_orders <- data_combined %>% filter(ref == "ts_blood") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "ts_blood")]
suffixes <- sub("ts_blood-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_ts_blood <- rnames[order(order_index)]

methods_orders <- data_combined %>% filter(ref == "kass_tumor") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "kass_tumor")]
suffixes <- sub("kass_tumor-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_kass_tumor <- rnames[order(order_index)]

methods_orders <- data_combined %>% filter(ref == "sc_pan_cancer") %>% group_by(ref, method) %>%  summarise(cor = median(ref_rho)) %>% arrange(-cor) %>% pull(method)
rnames <- rownames(data_matrix)[startsWith(rownames(data_matrix), "sc_pan_cancer")]
suffixes <- sub("sc_pan_cancer-", "", rnames)
order_index <- match(suffixes, methods_orders)
sorted_rnames_sc_pan_cancer<- rnames[order(order_index)]


rows_orders <- c(sorted_rnames_bp, sorted_rnames_lm22, sorted_rnames_kass_blood, sorted_rnames_ts_blood, sorted_rnames_kass_tumor,sorted_rnames_sc_pan_cancer)

ht_list <- ComplexHeatmap::Heatmap(data_matrix[rows_orders,vals_ordered],
        col = colorRampPalette(c("darkred" ,"red", "#FF6A6A", "white", "lightblue1","lightblue2", "blue", "darkblue"))(50),
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_side = "left",
        column_names_side = "top",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        na_col = "lightgray",
        split = group_factor,
        row_gap = unit(5, "mm"),
        heatmap_legend_param = list(title = "Weighted Spearman Correlation", direction = "horizontal",
                                    title_position = "topcenter", at = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), color_bar_len = 4)
        )
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")








# Order heatmap
bp_vals <- all_cors %>% filter(method == "xCell2" & ref == "bp") %>% pull(val) %>% unique()
lm22_vals <- all_cors %>% filter(method == "xCell2" & ref == "lm22") %>% pull(val) %>% unique()
kassblood_vals <- all_cors %>% filter(method == "xCell2" & ref == "kass_blood") %>% pull(val) %>% unique()
tsblood_vals <- all_cors %>% filter(method == "xCell2" & ref == "ts_blood") %>% pull(val) %>% unique()
kasstumor_vals <- all_cors %>% filter(method == "xCell2" & ref == "kass_tumor") %>% pull(val) %>% unique()
pancancer_vals <- all_cors %>% filter(method == "xCell2" & ref == "sc_pan_cancer") %>% pull(val) %>% unique()
only_general <- unique(c(bp_vals, lm22_vals))

only_blood <- c(kassblood_vals, tsblood_vals)
only_blood <- unique(only_blood[!only_blood %in% c(only_general, kasstumor_vals,pancancer_vals)])
only_tumor <- c(kasstumor_vals, pancancer_vals)
only_tumor <- unique(only_tumor[!only_tumor %in% c(only_general, only_blood)])
vals_ordered <- c(only_general, only_blood, only_tumor)



rows_orders <- c(sorted_rnames_bp, sorted_rnames_lm22)

prefixes <- sapply(strsplit(rows_orders, "-"), `[`, 1)
group_factor <- factor(prefixes, levels = unique(prefixes))

ht_list <- ComplexHeatmap::Heatmap(data_matrix[rows_orders, only_general],
                                   col = colorRampPalette(c("darkred" ,"red", "#FF6A6A", "white", "lightblue1","lightblue2", "blue", "darkblue"))(50),
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   row_names_side = "left",
                                   column_names_side = "top",
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   na_col = "lightgray",
                                   split = group_factor,
                                   row_gap = unit(5, "mm"),
                                   heatmap_legend_param = list(title = "Weighted Spearman Correlation", direction = "horizontal",
                                                               title_position = "topcenter", at = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), color_bar_len = 4)
)
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")



only_blood <- unique(c(kassblood_vals, tsblood_vals))

rows_orders <- c(sorted_rnames_kass_blood, sorted_rnames_ts_blood)

prefixes <- sapply(strsplit(rows_orders, "-"), `[`, 1)
group_factor <- factor(prefixes, levels = unique(prefixes))

ht_list <- ComplexHeatmap::Heatmap(data_matrix[rows_orders, only_blood],
                                   col = colorRampPalette(c("darkred" ,"red", "#FF6A6A", "white", "lightblue1","lightblue2", "blue", "darkblue"))(50),
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   row_names_side = "left",
                                   column_names_side = "top",
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   na_col = "lightgray",
                                   split = group_factor,
                                   row_gap = unit(5, "mm"),
                                   heatmap_legend_param = list(title = "Weighted Spearman Correlation", direction = "horizontal",
                                                               title_position = "topcenter", at = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), color_bar_len = 4)
)
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")




only_tumor <- unique(c(kasstumor_vals, pancancer_vals))

rows_orders <- c(sorted_rnames_kass_tumor, sorted_rnames_sc_pan_cancer)

prefixes <- sapply(strsplit(rows_orders, "-"), `[`, 1)
group_factor <- factor(prefixes, levels = unique(prefixes))

ht_list <- ComplexHeatmap::Heatmap(data_matrix[rows_orders, only_tumor],
                                   col = colorRampPalette(c("darkred" ,"red", "#FF6A6A", "white", "lightblue1","lightblue2", "blue", "darkblue"))(50),
                                   show_row_names = TRUE,
                                   show_column_names = TRUE,
                                   row_names_side = "left",
                                   column_names_side = "top",
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   na_col = "lightgray",
                                   split = group_factor,
                                   row_gap = unit(5, "mm"),
                                   heatmap_legend_param = list(title = "Weighted Spearman Correlation", direction = "horizontal",
                                                               title_position = "topcenter", at = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), color_bar_len = 4)
)
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
