library(tidyverse)

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

  p <- ggplot(data.in, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
    geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
    scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
    #geom_jitter(alpha = 0.2, size = 0.5) +
    theme_minimal() +
    #scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    #coord_flip() +
    coord_flip(ylim = c(0, 1)) +
    facet_wrap(~ ref, scales = "free", ncol = 1, labeller = labeller(ref = new_labels)) +
    tidytext::scale_x_reordered() +
    scale_x_discrete(labels=bold_xCell2_labels) +
    labs(x="", y="Weighted Spearman rho") +
    guides(fill=FALSE) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10))

  return(p)
}
makeHeatmap <- function(data.in, boxplot.data = benchmark_correlations, ref2use, is_mouse = FALSE){

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
    filter(ref %in% ref2use)

  data.in$val <- factor(data.in$val, levels = unique(data.in$val))

  boxplot.data <- boxplot.data %>%
    filter(ref %in% ref2use)

  complete_data <- data.in %>%
    complete(method, ref, val) %>%
    group_by(method, ref, val) %>%
    summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')


  complete_data$method <- factor(complete_data$method, levels = boxplot.data %>%
                                   group_by(method) %>%
                                   summarise(median_rho = median(ref_rho)) %>%
                                   arrange(median_rho) %>%
                                   pull(method))

  complete_data.splitted <- split(complete_data, complete_data$ref)


  plot_list <- lapply(ref2use, function(r){

    tmp_data <- complete_data.splitted[[r]]
    tmp_data$method <- factor(tmp_data$method, levels = benchmark_correlations %>%
                                filter(ref == r) %>%
                                group_by(method) %>%
                                summarise(median_rho = median(ref_rho)) %>%
                                arrange(median_rho) %>%
                                pull(method))

    ggplot(tmp_data, aes(x = val, y = method, fill = median_ref_rho)) +
      geom_tile(color = "black") +
      scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "#424242",
                           limit = c(0, 1), name="Median Weighted Spearman rho") +
      scale_x_discrete(position = "top") +
      scale_y_discrete(labels=bold_xCell2_labels) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 8, face = "bold"),
            axis.text.y = element_text(size = 12, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            legend.position = "bottom",
            legend.title = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position title above and center it
      labs(x = "", y = "", title = "") +
      coord_equal()


  })

  plot_list[2:length(plot_list)] <- plot_list[2:length(plot_list)] %>%
    map(~ .x + theme(axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.x=element_blank()))

  plot_list[1] <- plot_list[1] %>%
    map(~ .x +  guides(fill = FALSE))

  if (is_mouse) {
    plot_list[2] <- plot_list[2] %>%
      map(~ .x +  guides(fill = FALSE))
  }


  p <- gridExtra::arrangeGrob(egg::ggarrange(plots=plot_list, ncol=1),
                              heights=c(20,1))
  return(p)

}

# A - Validation datasets treemap --------------------

# /bigdata/almogangel/xCell2/dev_scripts/Figures/treemap.R

# B - Benchmarking results (Human) --------------------

source("/bigdata/almogangel/xCell2/dev_scripts/benchmarking_functions.R")

# Load benchmarking result file from fig1.R
xCell2results <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/xcell2_benchmark_results_spillAlpha_0.5.rds")

# Get correlations
benchmark_correlations <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_rhos = TRUE, cMethod = "spearman")


# Boxplot
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

new_labels <- c("bp" = "Blueprint-Encode", "kass_blood" = "Blood Immune Compendium",
                "kass_tumor" = "TME Compendium", "lm22" = "LM22", "sc_pan_cancer" = "Pan Cancer",
                "ts_blood"= "Tabula Sapiens Blood")

benchmark_correlations$ref <- factor(benchmark_correlations$ref, levels = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"))

p1 <- ggplot(benchmark_correlations, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  #geom_jitter(alpha = 0.2, size = 0.5) +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  #coord_flip() +
  coord_flip(ylim = c(0, 1)) +
  facet_wrap(~ ref, scales = "free", ncol = 1, labeller = labeller(ref = new_labels)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho") +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10))

p1


# Heatmap

benchmark_correlations_val <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3,
                                                      weight_rhos = TRUE, by_val = TRUE, ref2use = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"))


complete_data <- benchmark_correlations_val %>%
  complete(method, ref, val) %>%
  group_by(method, ref, val) %>%
  summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')


column_order <- c("BG_blood", "GSE106898", "GSE107011", "GSE107572",
                  "GSE107990", "GSE127813", "GSE130824",
                  "GSE59654","SDY311", "SDY420", "GSE65133", "GSE64385",
                  "GSE77343", "GSE77344", "SDY67",
                  "GSE65135", "GSE20300", "GSE115823",
                  "GSE64655", "GSE93722", "GSE120444",
                  "WU_ccRCC_RCCTC", "ccRCC_cytof_CD45+", "NSCLC_cytof")
complete_data$val <- factor(complete_data$val, levels = column_order)



complete_data.splitted <- split(complete_data, complete_data$ref)

color_palette <- rev(c(colorRampPalette(c("tomato4", "tomato", "orange", "gold", "yellow", "lightyellow", "white"))(50)))


plot_list <- lapply(c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"), function(r){
  tmp_data <- complete_data.splitted[[r]]
  tmp_data$method <- factor(tmp_data$method, levels = benchmark_correlations %>%
                                                 filter(ref == r) %>%
                                                 group_by(method) %>%
                                                 summarise(median_rho = median(ref_rho)) %>%
                                                 arrange(median_rho) %>%
                                                 pull(method))

  ggplot(tmp_data, aes(x = val, y = method, fill = median_ref_rho)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "#424242",
                         limit = c(0, 1), name="Median Weighted\nSpearman rho") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels=bold_xCell2_labels) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 8, face = "bold"),
          axis.text.y = element_text(size = 12, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold")) +
    #guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position title above and center it
    labs(x = "", y = "", title = "")

})

plot_list[2:length(plot_list)] <- plot_list[2:length(plot_list)] %>%
  map(~ .x + theme(axis.text.x=element_blank(),
                   axis.title.x=element_blank(),
                   axis.ticks.x=element_blank()))

plot_list[c(1,2,4,5,6)] <- plot_list[c(1,2,4,5,6)] %>%
  map(~ .x +  guides(fill = FALSE))



p2 <- gridExtra::arrangeGrob(egg::ggarrange(plots=plot_list, ncol=1),
                 heights=c(20,1))
p2

# C - Benchmarking results (mouse) --------------------

source("/bigdata/almogangel/xCell2/dev_scripts/benchmarking_functions.R")

# Load benchmarking result file from fig1.R
xCell2results <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/benchmark_xcell1_params/mouse/xcell2_benchmark_results_spillAlpha_0.5.rds")

# Get correlations
benchmark_correlations <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_rhos = TRUE, cMethod = "spearman")


# Boxplot
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

new_labels <- c("igd" = "ImmGenData", "mca_blood" = "Mouse Cell Atlas - Immune/Blood")

#benchmark_correlations$ref <- factor(benchmark_correlations$ref, levels = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"))

p1 <- ggplot(benchmark_correlations, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  #geom_jitter(alpha = 0.2, size = 0.5) +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  #coord_flip() +
  coord_flip(ylim = c(0, 1)) +
  facet_wrap(~ ref, scales = "free", ncol = 1, labeller = labeller(ref = new_labels)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman rho") +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10))

p1


# Heatmap

benchmark_correlations_val <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3,
                                                      weight_rhos = TRUE, by_val = TRUE, ref2use = c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"))


complete_data <- benchmark_correlations_val %>%
  complete(method, ref, val) %>%
  group_by(method, ref, val) %>%
  summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')


column_order <- c("BG_blood", "GSE106898", "GSE107011", "GSE107572",
                  "GSE107990", "GSE127813", "GSE130824",
                  "GSE59654","SDY311", "SDY420", "GSE65133", "GSE64385",
                  "GSE77343", "GSE77344", "SDY67",
                  "GSE65135", "GSE20300", "GSE115823",
                  "GSE64655", "GSE93722", "GSE120444",
                  "WU_ccRCC_RCCTC", "ccRCC_cytof_CD45+", "NSCLC_cytof")
complete_data$val <- factor(complete_data$val, levels = column_order)



complete_data.splitted <- split(complete_data, complete_data$ref)

color_palette <- rev(c(colorRampPalette(c("tomato4", "tomato", "orange", "gold", "yellow", "lightyellow", "white"))(50)))


plot_list <- lapply(c("bp", "lm22", "kass_blood", "ts_blood", "kass_tumor", "sc_pan_cancer"), function(r){
  tmp_data <- complete_data.splitted[[r]]
  tmp_data$method <- factor(tmp_data$method, levels = benchmark_correlations %>%
                              filter(ref == r) %>%
                              group_by(method) %>%
                              summarise(median_rho = median(ref_rho)) %>%
                              arrange(median_rho) %>%
                              pull(method))

  ggplot(tmp_data, aes(x = val, y = method, fill = median_ref_rho)) +
    geom_tile() +
    scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "#424242",
                         limit = c(0, 1), name="Median Weighted\nSpearman rho") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(labels=bold_xCell2_labels) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 8, face = "bold"),
          axis.text.y = element_text(size = 12, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold")) +
    #guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +  # Position title above and center it
    labs(x = "", y = "", title = "")

})

plot_list[2:length(plot_list)] <- plot_list[2:length(plot_list)] %>%
  map(~ .x + theme(axis.text.x=element_blank(),
                   axis.title.x=element_blank(),
                   axis.ticks.x=element_blank()))

plot_list[c(1,2,4,5,6)] <- plot_list[c(1,2,4,5,6)] %>%
  map(~ .x +  guides(fill = FALSE))



p2 <- gridExtra::arrangeGrob(egg::ggarrange(plots=plot_list, ncol=1),
                             heights=c(20,1))
p2


# D - Spillover ------------------------

# In fig1.R "spillover analysis" section



# Archive ------------------------------------



# Spillover
benchmark_correlations <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_results = 3, weight_rhos = TRUE, cMethod = "pearson", spillcors = TRUE)



new_labels <- c("bp" = "Blueprint-Encode", "kass_blood" = "Blood Immune Compendium",
                "kass_tumor" = "TME Compendium", "lm22" = "LM22", "sc_pan_cancer" = "Pan Cancer (scRNA-Seq)",
                "ts_blood"= "Tabula Sapiens Blood (scRNA-Seq)")

p2 <- ggplot(benchmark_correlations, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  #geom_jitter(alpha = 0.2, size = 0.5) +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  #coord_flip() +
  coord_flip(ylim = c(0, 1)) +
  facet_wrap(~ ref, scales = "free", ncol = 2, labeller = labeller(ref = new_labels)) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  labs(x="", y="Weighted Spearman Rho") +
  guides(fill=FALSE) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10))

p2









# Archive

benchmark_correlations_ref_bp <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "bp")
benchmark_correlations_ref_lm22 <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "lm22")
benchmark_correlations_ref_bp_lm22 <- rbind(benchmark_correlations_ref_bp, benchmark_correlations_ref_lm22)

ggplot(benchmark_correlations_ref_bp_lm22, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  labs(x="", y="", title = "") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_flip() +
  facet_wrap(~ ref, scales = "free", ncol = 1, ) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  guides(fill=FALSE) +
  theme(axis.text = element_text(size = 14, hjust = 1))


benchmark_correlations_val_bp_lm22 <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = TRUE, ref2use = c("bp", "lm22"))

color_palette <- rev(c(colorRampPalette(c("tomato", "orange", "yellow", "lightyellow", "white", "gray", "darkgray"))(50)))

method_order_upper <- rev(c("xCell2", "MCPcounter", "dtangle", "BayesPrism", "CIBERSORTx", "EPIC", "DeconRNASeq"))
method_order_lower <- rev(c("xCell2", "CIBERSORTx", "dtangle", "MCPcounter", "DeconRNASeq", "BayesPrism", "EPIC"))

complete_data <- benchmark_correlations_val_bp_lm22 %>%
  complete(method, ref, val) %>%
  group_by(method, ref, val) %>%
  summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')

complete_data <- complete_data %>%
  mutate(method = case_when(
    ref == "bp" ~ factor(method, levels = method_order_upper),
    ref == "lm22" ~ factor(method, levels = method_order_lower),
    TRUE ~ method  # In case other refs exist
  ))

complete_data %>%
  filter(ref == "bp") %>%
  mutate(method = factor(method, levels = method_order_upper)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")

complete_data %>%
  filter(ref == "lm22") %>%
  mutate(method = factor(method, levels = method_order_lower)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")


benchmark_correlations_ref_kb <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "kass_blood")
benchmark_correlations_ref_tsb <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "ts_blood")
benchmark_correlations_ref_kb_tsb <- rbind(benchmark_correlations_ref_kb, benchmark_correlations_ref_tsb)

ggplot(benchmark_correlations_ref_kb_tsb, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  labs(x="", y="", title = "") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_flip() +
  facet_wrap(~ ref, scales = "free", ncol = 1, ) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  guides(fill=FALSE) +
  theme(axis.text = element_text(size = 14, hjust = 1))

benchmark_correlations_val_kb_tsb <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = TRUE, ref2use = c("kass_blood", "ts_blood"))

color_palette <- rev(c(colorRampPalette(c("tomato", "orange", "yellow", "lightyellow", "white", "gray", "darkgray"))(50)))

method_order_upper <- rev(c("xCell2", "dtangle", "CIBERSORTx", "BayesPrism", "MCPcounter", "DeconRNASeq", "EPIC"))
method_order_lower <- rev(c("xCell2", "dtangle", "BayesPrism", "MCPcounter", "EPIC", "CIBERSORTx", "DeconRNASeq"))

complete_data <- benchmark_correlations_val_kb_tsb %>%
  complete(method, ref, val) %>%
  group_by(method, ref, val) %>%
  summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')

complete_data %>%
  filter(ref == "kass_blood") %>%
  mutate(method = factor(method, levels = method_order_upper)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")

complete_data %>%
  filter(ref == "ts_blood") %>%
  mutate(method = factor(method, levels = method_order_lower)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")


benchmark_correlations_ref_kt <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "kass_tumor")
benchmark_correlations_ref_spc <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = FALSE, ref2use = "sc_pan_cancer")
benchmark_correlations_ref_kt_spc <- rbind(benchmark_correlations_ref_kt, benchmark_correlations_ref_spc)

ggplot(benchmark_correlations_ref_kt_spc, aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#", fun = median), y = ref_rho, fill = is_xcell2)) +
  geom_boxplot(width=.5, outlier.colour=NA, coef = 0, show.legend = F, position = "dodge") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  labs(x="", y="", title = "") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  coord_flip() +
  facet_wrap(~ ref, scales = "free", ncol = 1, ) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels=bold_xCell2_labels) +
  guides(fill=FALSE) +
  theme(axis.text = element_text(size = 14, hjust = 1))


benchmark_correlations_val_kt_spc <- get_xcell2_correlations(vals2remove = train_ds, xCell2results = xCell2results, round_xcell2_results = 3, weight_rhos = TRUE, by_val = TRUE, ref2use = c("kass_tumor", "sc_pan_cancer"))

color_palette <- rev(c(colorRampPalette(c("tomato", "orange", "yellow", "lightyellow", "white", "gray", "darkgray"))(50)))

method_order_upper <- rev(c("xCell2", "dtangle", "MCPcounter", "CIBERSORTx", "BayesPrism", "EPIC", "DeconRNASeq"))
method_order_lower <- rev(c("dtangle", "xCell2", "MCPcounter", "CIBERSORTx", "BayesPrism", "EPIC", "DeconRNASeq"))

complete_data <- benchmark_correlations_val_kt_spc %>%
  complete(method, ref, val) %>%
  group_by(method, ref, val) %>%
  summarise(median_ref_rho = ifelse(is.na(median(ref_rho)), NA, median(ref_rho)), .groups = 'drop')

complete_data %>%
  filter(ref == "kass_tumor") %>%
  mutate(method = factor(method, levels = method_order_upper)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")

complete_data %>%
  filter(ref == "sc_pan_cancer") %>%
  mutate(method = factor(method, levels = method_order_lower)) %>%
  ggplot(., aes(x = val, y = method, fill = median_ref_rho)) +
  geom_tile() +
  facet_wrap(~ ref, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = color_palette, space = "Lab", na.value = "black",
                       limit = c(0, 1), name="Median Weighted Rho") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 12, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "", title = "")



