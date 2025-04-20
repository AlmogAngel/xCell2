library(tidyverse)


source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")

makeBoxPlot <- function(data.in, ref2use, is_mouse = FALSE){
  
  
  bold_xCell2_labels <- function(labels) {
    lapply(labels, function(label) {
      modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
      if (startsWith(modified_label, "xCell2")) {
        bquote(bold(.(modified_label)))
      } else {
        modified_label
      }
    })
  }
  
  data.in <- data.in %>%
    filter(ref %in% ref2use)
  
  if (is_mouse) {
    new_labels <- c("igd" = "ImmGenData", "mouse_rnaseq_data" = "MouseRNAseqData",
                    "tm_blood"= "Tabula Muris Blood")
  }else{
    new_labels <- c("bp" = "Blueprint-Encode", "kass_blood" = "Blood Immune Compendium",
                    "kass_tumor" = "TME Compendium", "lm22" = "LM22", "sc_pan_cancer" = "Pan Cancer",
                    "ts_blood"= "Tabula Sapiens Blood")
  }
  
  data.in$ref <- factor(data.in$ref, levels = ref2use)
  
  # Number of columns for the plot
  ncol2use <- ifelse(length(ref2use) == 2, 1, length(ref2use))
  
  # Round correlations
  data.in$ref_cor <- round(data.in$ref_cor, 3)
  
  # Rank the boxplots by the "adjusted median" = median - alpha × IQR
  adjusted_score_iqr <- function(x) {
    m <- median(x, na.rm = TRUE)
    iqr_val <- IQR(x, na.rm = TRUE)
    alpha <- 0.1
    m - alpha * iqr_val
  }
  
  
  p <- ggplot(data.in, aes(x = tidytext::reorder_within(x = method, by = ref_cor, within = ref, sep = "#", fun = adjusted_score_iqr), y = ref_cor, fill = is_xcell2)) +
    geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
    scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
    theme_minimal() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_flip(ylim = c(0, 1)) +
    facet_wrap(~ ref, scales = "free", ncol = ncol2use, labeller = labeller(ref = new_labels)) +
    tidytext::scale_x_reordered() +
    scale_x_discrete(labels=bold_xCell2_labels) +
    labs(x="", y="Weighted Spearman rho") +
    guides(fill=FALSE) +
    theme(strip.text = element_text(size = 20, face = "bold"),
          axis.title = element_text(face = "bold", size = 20),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 14))
  
  
  return(p)
}
makeHeatmap <- function(data.in, boxplot.data, ref2use, all_types_ref, is_mouse = FALSE, remove_legend = TRUE){
  
  bold_xCell2_labels <- function(labels) {
    lapply(labels, function(label) {
      modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
      if (startsWith(modified_label, "xCell2")) {
        bquote(bold(.(modified_label)))
      } else {
        modified_label
      }
    })
  }
  color_palette <- rev(c(colorRampPalette(c("tomato4", "tomato", "orange", "gold", "yellow", "lightyellow", "white"))(50)))
  
  data.in <- data.in %>%
    filter(ref %in% all_types_ref)
  
  data.in <- data.in %>%
    mutate(val = factor(val, levels = unique(val)))
  
  # Use boxplot data to order rows in heatmap
  boxplot.data <- boxplot.data %>%
    filter(ref %in% ref2use)
  
  complete_data <- data.in %>%
    complete(method, ref, val) %>%
    filter(ref %in% ref2use)
  

  
  complete_data$method <- factor(complete_data$method, levels = as.character(boxplot.data %>%
                                                                               group_by(method) %>%
                                                                               summarise(median_cor = median(ref_cor)) %>%
                                                                               arrange(median_cor) %>%
                                                                               pull(method)))
  
  
  complete_data <- complete_data %>%
    mutate(ref_val_cor_median = ifelse(ref_val_cor_median < 0, 0, ref_val_cor_median))
  
  
  p <- ggplot(complete_data, aes(x = val, y = method, fill = ref_val_cor_median)) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "#424242",
                         limit = c(0, 1), name="Median Weighted Spearman rho") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels=bold_xCell2_labels) +
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", size = 20),
          axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 8, face = "bold"),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 16, face = "bold"),
          legend.key.width = unit(1, "cm"),
          legend.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position title above and center it
    labs(x = "", y = "", title = "") +
    coord_equal()
  
  if (remove_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
  
}

# A - Validation datasets treemap --------------------

# /bigdata/almogangel/xCell2/dev_scripts/Figures/treemap.R


# B - Benchmarking results (Human) --------------------


# Load benchmarking result file from fig1D
xCell2results <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D/fig1d_xcell2_res_alpha_0.5.rds")


# Get correlations
benchmark_correlations_spr <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "spearman")
benchmark_correlations_val_spr <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3,
                                                      weight_cors = TRUE, by_val = TRUE, ref2use = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"), cMethod = "spearman")

benchmark_correlations_prs <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "pearson")
benchmark_correlations_val_prs <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3,
                                                          weight_cors = TRUE, by_val = TRUE, ref2use = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"), cMethod = "pearson")

# Boxplot
# fig3b_bp_lm22_spr: 800w x 1000h
fig3b_bp_lm22_spr <- makeBoxPlot(data.in = benchmark_correlations_spr, ref2use = c("bp", "lm22"))
print(fig3b_bp_lm22_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_lm22_spr.png", plot = fig3b_bp_lm22_spr, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_lm22_spr.pdf", plot = fig3b_bp_lm22_spr, device = "pdf", width = 8, height = 10)

fig3b_bp_lm22_prs <- makeBoxPlot(data.in = benchmark_correlations_prs, ref2use = c("bp", "lm22"))
print(fig3b_bp_lm22_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_lm22_prs.png", plot = fig3b_bp_lm22_prs, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_lm22_prs.pdf", plot = fig3b_bp_lm22_prs, device = "pdf", width = 8, height = 10)

# fig3b_tmec_panc_spr: 800w x 1000h
fig3b_tmec_panc_spr <- makeBoxPlot(data.in = benchmark_correlations_spr, ref2use = c("kass_tumor", "sc_pan_cancer"))
print(fig3b_tmec_panc_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_panc_spr.png", plot = fig3b_tmec_panc_spr, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_panc_spr.pdf", plot = fig3b_tmec_panc_spr, device = "pdf", width = 8, height = 10)

fig3b_tmec_panc_prs <- makeBoxPlot(data.in = benchmark_correlations_prs, ref2use = c("kass_tumor", "sc_pan_cancer"))
print(fig3b_tmec_panc_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_panc_prs.png", plot = fig3b_tmec_panc_prs, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_panc_prs.pdf", plot = fig3b_tmec_panc_prs, device = "pdf", width = 8, height = 10)

# fig3b_tmec_panc_spr: 800w x 1000h
fig3b_bloodc_panc_spr <- makeBoxPlot(data.in = benchmark_correlations_spr, ref2use = c("kass_blood"))
print(fig3b_bloodc_panc_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_panc_spr.png", plot = fig3b_bloodc_panc_spr, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_panc_spr.pdf", plot = fig3b_bloodc_panc_spr, device = "pdf", width = 8, height = 10)

fig3b_bloodc_panc_prs <- makeBoxPlot(data.in = benchmark_correlations_prs, ref2use = c("kass_blood"))
print(fig3b_bloodc_panc_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_panc_prs.png", plot = fig3b_bloodc_panc_prs, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_panc_prs.pdf", plot = fig3b_bloodc_panc_prs, device = "pdf", width = 8, height = 10)

# fig3b_tmec_panc_spr: 800w x 1000h
fig3b_tsblood_panc_spr <- makeBoxPlot(data.in = benchmark_correlations_spr, ref2use = c("ts_blood"))
print(fig3b_tsblood_panc_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_panc_spr.png", plot = fig3b_tsblood_panc_spr, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_panc_spr.pdf", plot = fig3b_tsblood_panc_spr, device = "pdf", width = 8, height = 10)

fig3b_tsblood_panc_prs <- makeBoxPlot(data.in = benchmark_correlations_prs, ref2use = c("ts_blood"))
print(fig3b_tsblood_panc_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_panc_prs.png", plot = fig3b_tsblood_panc_prs, device = "png", width = 8, height = 10, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_panc_prs.pdf", plot = fig3b_tsblood_panc_prs, device = "pdf", width = 8, height = 10)


# Heatmaps

# fig3b_bp_val_spr: 1400w x 1100h
fig3b_bp_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("bp"), all_types_ref = c("bp", "lm22")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_bp_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_val_spr.png", plot = fig3b_bp_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_val_spr.pdf", plot = fig3b_bp_val_spr, device = "pdf", width = 14, height = 11)

fig3b_bp_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("bp"), all_types_ref = c("bp", "lm22")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_bp_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_val_prs.png", plot = fig3b_bp_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bp_val_prs.pdf", plot = fig3b_bp_val_prs, device = "pdf", width = 14, height = 11)

# fig3b_lm22_val_spr: 1400w x 1100h
fig3b_lm22_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("lm22"), all_types_ref = c("bp", "lm22")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_lm22_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_lm22_val_spr.png", plot = fig3b_lm22_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_lm22_val_spr.pdf", plot = fig3b_lm22_val_spr, device = "pdf", width = 14, height = 11)

fig3b_lm22_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("lm22"), all_types_ref = c("bp", "lm22")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_lm22_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_lm22_val_prs.png", plot = fig3b_lm22_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_lm22_val_prs.pdf", plot = fig3b_lm22_val_prs, device = "pdf", width = 14, height = 11)

# fig3b_tmec_val_spr: 1400w x 1100h
fig3b_tmec_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("kass_tumor"), all_types_ref = c("kass_tumor", "sc_pan_cancer")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_tmec_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_val_spr.png", plot = fig3b_tmec_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_val_spr.pdf", plot = fig3b_tmec_val_spr, device = "pdf", width = 14, height = 11)

fig3b_tmec_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("kass_tumor"), all_types_ref = c("kass_tumor", "sc_pan_cancer")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_tmec_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_val_prs.png", plot = fig3b_tmec_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tmec_val_prs.pdf", plot = fig3b_tmec_val_prs, device = "pdf", width = 14, height = 11)

# fig3b_panc_val_spr: 1400w x 1100h
fig3b_panc_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("sc_pan_cancer"), all_types_ref = c("kass_tumor", "sc_pan_cancer")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_panc_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_panc_val_spr.png", plot = fig3b_panc_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_panc_val_spr.pdf", plot = fig3b_panc_val_spr, device = "pdf", width = 14, height = 11)

fig3b_panc_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("sc_pan_cancer"), all_types_ref = c("kass_tumor", "sc_pan_cancer")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_panc_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_panc_val_prs.png", plot = fig3b_panc_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_panc_val_prs.pdf", plot = fig3b_panc_val_prs, device = "pdf", width = 14, height = 11)

# fig3b_bloodc_val_spr: 1400w x 1100h
fig3b_bloodc_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("kass_blood"), all_types_ref = c("kass_blood")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_bloodc_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_val_spr.png", plot = fig3b_bloodc_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_val_spr.pdf", plot = fig3b_bloodc_val_spr, device = "pdf", width = 14, height = 11)

fig3b_bloodc_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("kass_blood"), all_types_ref = c("kass_blood")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_bloodc_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_val_prs.png", plot = fig3b_bloodc_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_bloodc_val_prs.pdf", plot = fig3b_bloodc_val_prs, device = "pdf", width = 14, height = 11)

# fig3b_tsblood_val_spr: 1400w x 1100h
fig3b_tsblood_val_spr <- makeHeatmap(data.in = benchmark_correlations_val_spr, boxplot.data = benchmark_correlations_spr, ref2use = c("ts_blood"), all_types_ref = c("ts_blood")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_tsblood_val_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_val_spr.png", plot = fig3b_tsblood_val_spr, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_val_spr.pdf", plot = fig3b_tsblood_val_spr, device = "pdf", width = 14, height = 11)

fig3b_tsblood_val_prs <- makeHeatmap(data.in = benchmark_correlations_val_prs, boxplot.data = benchmark_correlations_prs, ref2use = c("ts_blood"), all_types_ref = c("ts_blood")) + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
print(fig3b_tsblood_val_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_val_prs.png", plot = fig3b_tsblood_val_prs, device = "png", width = 14, height = 11, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3b_tsblood_val_prs.pdf", plot = fig3b_tsblood_val_prs, device = "pdf", width = 14, height = 11)


# C - All benchmarking results (human/mouse) --------------------

# Human
# Boxplot - all refs
method_order <- benchmark_correlations_spr %>%
  group_by(method) %>%
  summarise(ref_cor = median(ref_cor)) %>%
  arrange(-ref_cor) %>%
  pull(method)

benchmark_correlations_spr$method <- factor(benchmark_correlations_spr$method, levels = rev( method_order))

pastel_colors <- c("#CD5555", "#3A5FCD", "#43CD80", "#FFB90F", "#BF3EFF", "#FF6EB4")

new_labels <- c("bp" = "Blueprint-Encode", "kass_blood" = "Blood Immune Compendium",
                "kass_tumor" = "TME Compendium", "lm22" = "LM22", "sc_pan_cancer" = "Pan Cancer",
                "ts_blood"= "Tabula Sapiens Blood")

bold_xCell2_labels <- function(labels) {
  lapply(labels, function(label) {
    modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
    if (startsWith(modified_label, "xCell2")) {
      bquote(bold(.(modified_label)))
    } else {
      modified_label
    }
  })
}

fig3c_human_spr <- benchmark_correlations_spr %>%
  ggplot(., aes(x=method, y=ref_cor, col=ref)) +
  geom_jitter(width = 0.2, size=1.5, alpha=0.5) +
  geom_boxplot(aes(group = method), color = "black", fill = "gray", alpha = 0.3, outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = pastel_colors, labels = new_labels) +
  #stat_summary(fun = median, geom = "point", aes(group = method), size = 5) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_flip(ylim = c(0, 1)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho", col="Reference:") +
  guides(fill=FALSE) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    legend.position = "right")

# fig3c_human_spr: 1000w x 700h
print(fig3c_human_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_human_spr.png", plot = fig3c_human_spr, device = "png", width = 10, height = 7, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_human_spr.pdf", plot = fig3c_human_spr, device = "pdf", width = 10, height = 7)


method_order <- benchmark_correlations_prs %>%
  group_by(method) %>%
  summarise(ref_cor = median(ref_cor)) %>%
  arrange(-ref_cor) %>%
  pull(method)

benchmark_correlations_prs$method <- factor(benchmark_correlations_prs$method, levels = rev( method_order))


fig3c_human_prs <- benchmark_correlations_prs %>%
  ggplot(., aes(x=method, y=ref_cor, col=ref)) +
  geom_jitter(width = 0.2, size=1.5, alpha=0.5) +
  geom_boxplot(aes(group = method), color = "black", fill = "gray", alpha = 0.3, outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = pastel_colors, labels = new_labels) +
  #stat_summary(fun = median, geom = "point", aes(group = method), size = 5) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_flip(ylim = c(0, 1)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho", col="Reference:") +
  guides(fill=FALSE) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    legend.position = "right")

# fig3c_human_prs: 1000w x 700h
print(fig3c_human_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_human_prs.png", plot = fig3c_human_prs, device = "png", width = 10, height = 7, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_human_prs.pdf", plot = fig3c_human_prs, device = "pdf", width = 10, height = 7)


# Set mouse = TRUE !!!!
source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")

# Run xCell2 with mouse validations
params2use <- list(human2mouse = FALSE,
                   num_threads = 30,
                   min_pb_cells = 20,
                   min_pb_samples = 10,
                   min_sc_genes = 1e4,
                   use_ontology = TRUE,
                   return_signatures = FALSE,
                   return_analysis = FALSE,
                   use_sillover = TRUE,
                   spillover_alpha = 0.5,
                   top_spill_value = 0.5)

dir2use <- "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig2C"
xcell2_res_mouse <- paste0(dir2use, "/fig2c_xcell2_mouse_res.rds")

get_xcell2_benchmark_results(save_object = TRUE, dir = dir2use,
                             params = params2use, ncores = 4, output_name = xcell2_res_mouse)


xCell2results <- readRDS(xcell2_res_mouse)
benchmark_correlations_mouse_spr <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "spearman")
benchmark_correlations_mouse_prs <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "pearson")


# Plot
bold_xCell2_labels <- function(labels) {
  lapply(labels, function(label) {
    modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
    if (startsWith(modified_label, "xCell2")) {
      bquote(bold(.(modified_label)))
    } else {
      modified_label
    }
  })
}

method_order <- benchmark_correlations_mouse_spr %>%
  group_by(method) %>%
  summarise(ref_cor = median(ref_cor)) %>%
  arrange(-ref_cor) %>%
  pull(method)

benchmark_correlations_mouse_spr$method <- factor(benchmark_correlations_mouse_spr$method, levels = rev( method_order))

pastel_colors <- c("#CD5555", "#3A5FCD", "#43CD80")

new_labels <- c("igd" = "ImmGenData", "mouse_rnaseq_data" = "MouseRNAseqData",
                "tm_blood"= "Tabula Muris Blood")

fig3c_mouse_spr <- benchmark_correlations_mouse_spr %>%
  ggplot(., aes(x=method, y=ref_cor, col=ref)) +
  geom_jitter(width = 0.2, size=3, alpha=0.8) +
  geom_boxplot(aes(group = method), color = "black", fill = "gray", alpha = 0.3, outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = pastel_colors, labels = new_labels) +
  #stat_summary(fun = median, geom = "point", aes(group = method), size = 5) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_flip(ylim = c(0, 1)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho", col="Reference:") +
  guides(fill=FALSE) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    legend.position = "right")

# fig3c_mouse_spr: 1000w x 700h
print(fig3c_mouse_spr)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_mouse_spr.png", plot = fig3c_mouse_spr, device = "png", width = 10, height = 7, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_mouse_spr.pdf", plot = fig3c_mouse_spr, device = "pdf", width = 10, height = 7)



method_order <- benchmark_correlations_mouse_prs %>%
  group_by(method) %>%
  summarise(ref_cor = median(ref_cor)) %>%
  arrange(-ref_cor) %>%
  pull(method)

benchmark_correlations_mouse_prs$method <- factor(benchmark_correlations_mouse_prs$method, levels = rev( method_order))


fig3c_mouse_prs <- benchmark_correlations_mouse_prs %>%
  ggplot(., aes(x=method, y=ref_cor, col=ref)) +
  geom_jitter(width = 0.2, size=3, alpha=0.8) +
  geom_boxplot(aes(group = method), color = "black", fill = "gray", alpha = 0.3, outlier.shape = NA, width = 0.5) +
  scale_color_manual(values = pastel_colors, labels = new_labels) +
  #stat_summary(fun = median, geom = "point", aes(group = method), size = 5) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  coord_flip(ylim = c(0, 1)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho", col="Reference:") +
  guides(fill=FALSE) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    legend.position = "right")

# fig3c_mouse_prs: 1000w x 700h
print(fig3c_mouse_prs)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_mouse_prs.png", plot = fig3c_mouse_prs, device = "png", width = 10, height = 7, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3c_mouse_prs.pdf", plot = fig3c_mouse_prs, device = "pdf", width = 10, height = 7)


# D - Spillover ------------------------

# Set back mouse = FALSE !!!!
source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")

# Get other method spill correlations

other_methods_spillcor <- parallel::mclapply(1:nrow(human.predicted), function(i){
  
  # Load ref, val, and ground truth
  valType <- human.predicted[i,]$val_type
  valDataset <- human.predicted[i,]$val_dataset
  refName <- human.predicted[i,]$ref_name
  refType <- human.predicted[i,]$ref_type
  ref.in <- human.refs[[refType]][[refName]]
  dep_list <- getDependencies(ref.in$lineage_file)
  gep_mat <- makeGEPMat(ref.in$ref, ref.in$labels)
  cor_mat <- getCellTypeCorrelation(gep_mat, refType)
  truth <- human.vals$truth[[valType]][[valDataset]]
  truth <- round(truth, 3)
  res <- human.predicted[i,]$res[[1]]
  res <- round(res, 3)
  
  celltypes <- intersect(rownames(res), rownames(truth))
  samples <- intersect(colnames(res), colnames(truth))
  
  # Using correlation find the most similar cell type
  ct2most_simillar <- sapply(celltypes, function(ct){
    celltypes2use <-  celltypes[!celltypes %in% c(ct, unlist(dep_list[ct]))]
    names(sort(cor_mat[ct, celltypes2use], decreasing = TRUE))[1]
  })
  
  y <- human.predicted[i,] %>%
    select(method:val_dataset)
  
  spill_cor_res <- tibble(sig_ct = celltypes, most_sim_truth_ct = ct2most_simillar)
  spill_cor_res <- cbind(y, spill_cor_res)
  
  spill_cor_res %>%
    rowwise() %>%
    mutate(spill_cor = cor(res[sig_ct, samples], truth[most_sim_truth_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
    mutate(direct_cor = cor(res[sig_ct, samples], truth[sig_ct, samples], method = "pearson", use = "pairwise.complete.obs")) %>%
    ungroup() %>%
    return(.)
  
  
}, mc.cores = 20) 
other_methods_spillcor <- bind_rows(other_methods_spillcor)
saveRDS(other_methods_spillcor, "/bigdata/almogangel/xCell2_dev/paper_figures/data/fig2D/other_methods_spillcor.rds")


other_methods_spillcor <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig2D/other_methods_spillcor.rds")

other_methods_spillcor[is.na(other_methods_spillcor$spill_cor),]$spill_cor <- 0
other_methods_spillcor[is.na(other_methods_spillcor$direct_cor),]$direct_cor <- 0
other_methods_spillcor[other_methods_spillcor$spill_cor < 0,]$spill_cor <- 0
other_methods_spillcor[other_methods_spillcor$direct_cor < 0,]$direct_cor <- 0


other_methods_spillcor <- other_methods_spillcor %>%
  mutate(delta_cor = direct_cor - spill_cor)
other_methods_spillcor$spill_alpha <- "0"

# Load xCell2 spillover results
spill_cor_res_table <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D/xcell2_spillover_correlations.rds")
spill_cor_res_table[is.na(spill_cor_res_table$spill_cor),]$spill_cor <- 0
spill_cor_res_table[is.na(spill_cor_res_table$direct_cor),]$direct_cor <- 0
spill_cor_res_table <- spill_cor_res_table %>%
  mutate(delta_cor = direct_cor - spill_cor)
spill_cor_res_table <- spill_cor_res_table[,colnames(other_methods_spillcor)]

# Combine
spill_res_final <- rbind(spill_cor_res_table, other_methods_spillcor[,colnames(spill_cor_res_table)])
spill_res_final$spill_alpha <- factor(spill_res_final$spill_alpha)

# Plot
data2plot <- spill_res_final %>%
  filter(spill_alpha %in% c(0, 0.5)) %>%
  mutate(method = ifelse(method == "xCell2", paste0(method, " (alpha = ", spill_alpha, ")"), method))

method_levels <- data2plot %>% group_by(method) %>% summarise(delta_cor = median(delta_cor)) %>% arrange(delta_cor) %>% pull(method)
method_levels[1:3] <- rev(method_levels[1:3])
data2plot$method <- factor(data2plot$method, levels = method_levels)
data2plot$is_xcell <- ifelse(startsWith(as.character(data2plot$method), "xCell2"), "yes", "no")

bold_xCell2_labels <- function(labels) {
  lapply(labels, function(label) {
    modified_label <- str_remove(label, "#.*")  # Remove everything after "#"
    if (startsWith(modified_label, "xCell2")) {
      bquote(bold(.(modified_label)))
    } else {
      modified_label
    }
  })
}


fig3d <- ggplot(data2plot, aes(x=method, y=delta_cor, fill=is_xcell)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  coord_flip(ylim = c(0, 0.8)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Delta Pearson r") +
  guides(fill=FALSE) +
  theme(axis.title = element_text(face = "bold", size = 20),
        strip.text = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 20, face = "bold")
  )

# fig3d: 650w x 800h
print(fig3d)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3d.png", plot = fig3d, device = "png", width = 6.5, height = 8, dpi = 300)
ggsave(filename = "/bigdata/almogangel/xCell2_dev/paper_figures/fig3d.pdf", plot = fig3d, device = "pdf", width = 6.5, height = 8)

