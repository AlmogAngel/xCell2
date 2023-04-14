
library(tidyverse)
library(Rtsne)

# New -------
all_cors.out <- bind_rows(all_cors.out)

ct2use <- names(sort(table(all_cors.out$celltype), decreasing = TRUE))[1:18]
ct2use <- ct2use[!ct2use %in% c("Other", "Non plasma B-cells")]

median_cors <- all_cors.out %>%
  filter(celltype %in% ct2use) %>%
  # filter(ds_type == "Blood datasets" & !method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  filter(cor_method == "Spearman") %>%
  dplyr::select(method, celltype, cor_score) %>%
  group_by(method, celltype) %>%
  summarise(median_cor = median(cor_score, na.rm = T))

df <- pivot_wider(median_cors, names_from = method, values_from = median_cor)
df <- data.frame(df[,-1], row.names = df$celltype, check.names = FALSE)
df[df < 0] <- NA
methods_sorted <- names(sort(apply(df, 2, function(x){median(x, na.rm = T)}), decreasing = F))

median_cors$method <- factor(median_cors$method, levels = methods_sorted)

median_cors %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=median_cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=celltype), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1), labels = as.character(seq(0,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(1,3,5,6,4,9,2,8,9,11)],
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median)", title = "All Validation Datasets", x = NULL, colour = NULL, fill = NULL)




#


ct2use <- names(sort(table(cors.out$celltype), decreasing = TRUE))[1:10]
ct2use <- ct2use[!ct2use %in% c("Other", "Non plasma B-cells")]


median_cors <- cors.out %>%
  filter(celltype %in% ct2use) %>%
  # filter(ds_type == "Blood datasets" & !method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  filter(cor_method == "Spearman") %>%
  dplyr::select(method, celltype, cor_score) %>%
  group_by(method, celltype) %>%
  summarise(median_cor = median(cor_score, na.rm = T))

df <- pivot_wider(median_cors, names_from = method, values_from = median_cor)
df <- data.frame(df[,-1], row.names = df$celltype, check.names = FALSE)
df[df < 0] <- NA
methods_sorted <- names(sort(apply(df, 2, function(x){median(x, na.rm = T)}), decreasing = F))

median_cors$method <- factor(median_cors$method, levels = methods_sorted)

median_cors %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=median_cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=celltype), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  #scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(1,3,5,6,4,9,2,8,9,11)],
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median)", title = "BG Blood dataset", x = NULL, colour = NULL, fill = NULL)



# ILANIT 2023 ----
# tSNE for reference datasets ----

# BP ref
set.seed(123)
tSNE_fit <- t(bp@assays@data$logcounts) %>%
  Rtsne(num_threads = 0, check_duplicates = FALSE)

tSNE_fit$Y %>%
  as.data.frame() %>%
  dplyr::rename(tSNE1=1,
         tSNE2=2) %>%
  cbind(label = bp@colData$label.main) %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = label))+
  geom_point(size=3) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", RColorBrewer::brewer.pal(8, "Set3"))) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

# Tumor ref
set.seed(123)
tSNE_fit_tumor <- t(tumor_ref) %>%
  Rtsne(num_threads = 0, check_duplicates = FALSE)
saveRDS(tSNE_fit_tumor, "../xCell2.0/tsne_tumor_fit.rds")

tSNE_fit_tumor$Y %>%
  as.data.frame() %>%
  dplyr::rename(tSNE1=1,
                tSNE2=2) %>%
  cbind(label = tumor_labels$label) %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = label))+
  geom_point(size=2) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", RColorBrewer::brewer.pal(8, "Set3"))) +
                                theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))

# Blood ref
set.seed(123)
tSNE_fit_blood <- t(blood_ref) %>%
  Rtsne(num_threads = 0, check_duplicates = FALSE)
saveRDS(tSNE_fit_blood, "../xCell2.0/tsne_blood_fit.rds")

tSNE_fit_blood$Y %>%
  as.data.frame() %>%
  dplyr::rename(tSNE1=1,
                tSNE2=2) %>%
  cbind(label = blood_labels$label) %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = label))+
  geom_point(size=3) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", RColorBrewer::brewer.pal(8, "Set3"))) +
                                theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"))


# Heatmap of validation datasets -----
mix_dir <- "../xCell2.0/Kassandara_data/expressions/"
ds <- gsub(".tsv", "", list.files(mix_dir))

mix <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", ds[1], ".tsv"), check.names = FALSE, row.names = 1, header = TRUE)
ds <- ds[-1]

for(d in ds){
  mix_tmp <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", d, ".tsv"), check.names = FALSE, row.names = 1, header = TRUE)
  if (nrow(mix_tmp) != nrow(mix)) {
    genes2use <- interaction(rownames(mix_tmp), rownames(mix))
    mix <- mix[genes2use,]
    mix_tmp <- mix_tmp[genes2use,]
  }
  mix <- cbind(mix, mix_tmp)
}

pheatmap::pheatmap(mix, cluster_rows=F, cluster_cols=F, scale = "none", col= RColorBrewer::brewer.pal(11, "RdBu"))


# Benchmarking examples ----

ctoi <- "Cancer cells"
cors.df <- all_cors.out %>%
  filter(ds_type == "Tumor datasets") %>%
  filter(!method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ungroup() %>%
  filter(dataset != "Tonsils") %>%
  filter(cor_method == "Spearman") %>%
  filter(celltype == ctoi)

methods_sorted <- dplyr::select(cors.df, method, cor_score) %>%
  group_by(method) %>%
  summarise(median = median(cor_score)) %>%
  arrange(-median) %>%
  pull(method)
cors.df$method <- factor(cors.df$method, levels = methods_sorted)

cors.df %>%
  separate(method, sep = "-", into = "method") %>%
  ggplot(., aes(x=method, y=cor_score, fill=is_xcell2)) +
  geom_boxplot(position = position_dodge(1), alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(6, 10, 2, 4)],
                              "#424242", "#00F5FF", "#FF3E96"), name = "Validation Dataset") +
                                scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Cancer-cells (Tumor Purity)", x = NULL, colour = NULL, fill = NULL)



cors.df %>%
  filter(dataset == "SC_NSCLC" & method == "FARDEEP-Rel")
truth <- read.table("../xCell2.0/Kassandara_data/cell_values/SC_NSCLC.tsv")
predi <- read.table("../xCell2.0/Kassandara_data/predicted_by_algorithms/fardeep_relative/SC_NSCLC_predicted_by_fardeep_relative.tsv")
cor(as.numeric(truth["Macrophages",]), as.numeric(predi["Macrophages", colnames(truth)]), method = "spearman")

data.frame(scores=as.numeric(predi["Macrophages", colnames(truth)]), true_fracs=as.numeric(truth["Macrophages",])) %>%
  ggplot(., aes(x=scores, y=true_fracs)) +
  geom_point(size=5) +
  stat_smooth(method = "lm") +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "True fraction", title = NULL, x = "Predicted scores", colour = NULL, fill = NULL)




ctoi <- c("CD8+ T-cells PD1 low", "CD8+ T-cells PD1 high")

all_cors.out %>%
  filter(ds_type == "Tumor datasets") %>%
  filter(!method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ungroup() %>%
  filter(dataset != "Tonsils") %>%
  filter(cor_method == "Spearman") %>%
  filter(celltype %in% ctoi) %>%
  mutate(expression = ifelse(celltype == "CD8+ T-cells PD1 low", "CD8 + PD1 Low", "CD8+ PD1 High")) %>%
  ggplot(., aes(x=expression, y=cor_score, fill=expression)) +
  geom_boxplot(position = position_dodge(1), alpha = 0.6, outlier.shape = NA,  alpha = 0.5) +
  geom_point(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
  scale_color_manual(values=RColorBrewer::brewer.pal(10, "Paired")[c(10,2,4,5,6,8)], name = "Validation Dataset") +
                                scale_fill_manual( values = c("CD8 + PD1 Low"="tomato", "CD8+ PD1 High"="#43CD80"), guide = FALSE) +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "xCell 2.0 perfornace on \nCD8+ PD1 low/high", x = NULL, colour = NULL, fill = NULL)


tumor_sigs <- readRDS("../xCell2.0/xcell2_tumor_refsigs_16_2_23_500genes_gubbs05.rds")
getcd8scores <- function(ds, ct){
  mix_in <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", ds, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)
  scores <- score_xCell2(mix = mix_in, sigs = tumor_sigs)
  scores <- scores[ct,]

  truth <- read_tsv(paste0("../xCell2.0/Kassandara_data/cell_values/", ds, ".tsv")) %>%
    dplyr::rename(celltype = 1) %>% # Fix for some datasets
    filter(!endsWith(celltype, "l)")) %>%
    mutate(celltype = plyr::mapvalues(celltype, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE))
  truth <- data.frame(truth[,-1], row.names = truth$celltype, check.names = F)
  truth <- truth[ct,]

  samples2use <- intersect(names(truth), names(scores))
  scores <- as.numeric(scores[samples2use])
  truth <- as.numeric(truth[samples2use])

  return(list("scores" = scores, "truth" = truth))

}


cd8pd1.df <- all_cors.out %>%
  filter(celltype %in% ctoi & method != "Kassandara") %>%
  dplyr::select(celltype, dataset) %>%
  unique() %>%
  rowwise() %>%
  mutate(scores = list(getcd8scores(dataset, celltype)[[1]]),
         truth = list(getcd8scores(dataset, celltype)[[2]]))

cd8pd1.df %>%
  dplyr::select(-method) %>%
  unnest() %>%
  filter(celltype == ctoi[2]) %>%
  ggplot(., aes(x=scores, y=truth)) +
  geom_point(size=5, col=RColorBrewer::brewer.pal(10, "Paired")[6]) +
  stat_smooth(method = "lm") +
  theme_linedraw() +
  ggpubr::stat_cor(aes(label = ..r.label..), method = "spearman", size=5) +
  facet_wrap(~dataset, scales = "free") +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "True fraction", title = ctoi[2], x = "Predicted scores", colour = NULL, fill = NULL)


 lapply(pd1_ds, function(ds){
   truth <- read_tsv(paste0("../xCell2.0/Kassandara_data/cell_values/", ds, ".tsv")) %>%
     dplyr::rename(celltype = 1) %>% # Fix for some datasets
     filter(!endsWith(celltype, "l)")) %>%
     filter(celltype != "Respiratory_cells") %>%
     filter(celltype != "Tumor KI67+") %>%
     mutate(celltype = plyr::mapvalues(celltype, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE))
   truth <- data.frame(truth[,-1], row.names = truth$celltype, check.names = F)
   tumor_sigs <- readRDS("../xCell2.0/xcell2_tumor_refsigs_16_2_23_500genes_gubbs05.rds")
   mix_in <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", pd1_ds[1], "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

  scores <- score_xCell2(mix = mix_in, sigs = tumor_sigs)
  cor(as.numeric(truth["CD8+ T-cells PD1 high",]), as.numeric(scores["CD8+ T-cells PD1 high", colnames(truth)]), method = "spearman")

  data.frame(scores=as.numeric(scores["CD8+ T-cells PD1 high", colnames(truth)]), true_fracs=as.numeric(truth["CD8+ T-cells PD1 high",])) %>%
    ggplot(., aes(x=scores, y=true_fracs)) +
    geom_point(size=5) +
    stat_smooth(method = "lm") +
    theme_linedraw() +
    theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
          panel.grid.minor = element_line(colour = "white"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 40, hjust=1, face = "bold"),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA)) +
    labs(y = "True fraction", title = NULL, x = "Predicted scores", colour = NULL, fill = NULL)
})


# All benchmarking results ----

# median_cors <- all_cors.out %>%
#   filter(ds_type == "Blood datasets" & !method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
#   filter(cor_method == "Spearman") %>%
#   filter(!celltype %in% c("Lymphocytes", "Memory T-helpers", "cDC", "Tregs")) %>%
#   dplyr::select(method, celltype, cor_score) %>%
#   group_by(method, celltype) %>%
#   summarise(median_cor = median(cor_score, na.rm = T)) %>%
#   pivot_wider(names_from = method, values_from = median_cor)
#
# median_cors.df <- data.frame(median_cors[,-1], row.names = median_cors$celltype)
# methods_sorted <- names(sort(apply(median_cors.df, 2, function(x){median(x, na.rm = T)}), decreasing = T))
# median_cors.mat <- as.matrix(median_cors.df[,methods_sorted])
# colors <- c(RColorBrewer::brewer.pal(11, "RdBu")[1:4], RColorBrewer::brewer.pal(10, "RdBu")[5:9], "black")
# corrplot::corrplot(median_cors.mat, method = 'circle', col.lim = c(0, 1),  col=colors, is.corr = FALSE)
#

 all_cors.out <- readRDS("/home/almogangel/xCell2_git/Data/benchmarking_data/kass_benchmarking_cors_sref_20ds_simulations.rds")


# B val
median_cors <- all_cors.out %>%
  filter(ds_type == "Blood datasets" & !method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  filter(cor_method == "Spearman") %>%
  filter(!celltype %in% c("Tregs", "Th1 cells", "Th2 cells", "Th17 cells")) %>%
  dplyr::select(method, celltype, cor_score) %>%
  group_by(method, celltype) %>%
  summarise(median_cor = median(cor_score, na.rm = T))

df <- pivot_wider(median_cors, names_from = method, values_from = median_cor)
df <- data.frame(df[,-1], row.names = df$celltype, check.names = FALSE)
df[df < 0] <- NA
methods_sorted <- names(sort(apply(df, 2, function(x){median(x, na.rm = T)}), decreasing = F))

median_cors$method <- factor(median_cors$method, levels = methods_sorted)

median_cors %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=median_cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), alpha = 0.6, outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=celltype), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1), labels = as.character(seq(0,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(1,3,5,6,4,9,2,8)],
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
                                scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median)", title = "Blood Validation Datasets", x = NULL, colour = NULL, fill = NULL)


# T ref
median_cors <- all_cors.out %>%
  filter(ds_type == "Tumor datasets" & !method %in% c("Kassandara", "quanTIseq", "quanTIseq-T")) %>%
  filter(cor_method == "Spearman") %>%
  filter(!dataset %in% c("ccRCC_cytof_CD45+", "SC_NSCLC")) %>%
  filter(!celltype %in% c("Tregs", "Monocytes", "Neutrophils", "Lymphocytes", "CD8+ T-cells PD1 high", "CD8+ T-cells PD1 low")) %>%
  dplyr::select(method, celltype, cor_score) %>%
  group_by(method, celltype) %>%
  summarise(median_cor = median(cor_score, na.rm = T))

df <- pivot_wider(median_cors, names_from = method, values_from = median_cor)
df <- data.frame(df[,-1], row.names = df$celltype, check.names = FALSE)
df[df < 0] <- NA
methods_sorted <- names(sort(apply(df, 2, function(x){median(x, na.rm = T)}), decreasing = F))

median_cors$method <- factor(median_cors$method, levels = methods_sorted)

median_cors %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=median_cor)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), alpha = 0.6, outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=celltype), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1), labels = as.character(seq(0,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired")[c(1,3,5,6,4,9,2)],
                              "#424242", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
                                scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r (median)", title = "Tumor Validation Datasets", x = NULL, colour = NULL, fill = NULL)

