library(tidyverse)

xCell2results <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.resMat.cyto.res.20dec.rds")

xCell2results <- bind_rows(xCell2results) %>%
  filter(topGenesFrac == 1)

xCell2results$is.filt <- "no"


filtsigs <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.sigsFilt.cyto.res.20dec.rds")
valrefs <- unlist(lapply(1:nrow(vals.refs.res), function(i){
  paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
}))
names(filtsigs) <- valrefs


for (i in 1:length(valrefs)) {

  # file <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", valrefs[i], "_sigsFilt_20dec.rds")
  # if (file.exists(file)) {
  #   sigs.filt <- readRDS(file)
  # }else{
  #   next
  # }

  filt.data <- xCell2results[xCell2results$sig_name %in% filtsigs[[valrefs[i]]] & xCell2results$val_ref == valrefs[i], ]
  filt.data$is.filt <- "yes"

  xCell2results <- rbind(filt.data, xCell2results)

}


all_vals <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "GSE53655", "SDY311",
              "SDY420", "SDY67", "GSE64385", "GSE65133", "GSE106898", "GSE64655",
              "GSE59654", "GSE107990", "GSE20300", "GSE65135", "GSE77343", "GSE77344",
              "ccRCC_cytof_CD45+", "NSCLC_cytof", "WU_ccRCC_RCCTC", "GSE120444", "GSE115823", "GSE93722")

vals_type <- c("rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq", "array",
               "array", "rnaseq", "array", "array", "array", "rnaseq",
               "array", "array", "array", "array", "array", "array",
               "rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq", "rnaseq")

val2type <- setNames(vals_type, all_vals)


xCell2results %>%
  mutate(ref = str_extract(val_ref, "(kass_blood|kass_tumor|bp|lm22|ts_blood|sc_pan_cancer)$"),
   val = str_replace(val_ref, "(kass_blood|kass_tumor|bp|lm22|ts_blood|sc_pan_cancer)$", "")) %>%
  mutate(val = str_sub(val, 1, -2)) %>%
  mutate(val_type = val2type[val]) %>%
  group_by(ref, celltype, is.filt) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=ref, y=cor, fill=is.filt)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


xCell2results_processed <- xCell2results %>%
  group_by(val_ref, celltype, is.filt) %>%
  summarise(cor = mean(cor), .groups = "drop") %>%
  group_by(val_ref, celltype) %>%
  summarise(delta_cor = cor[is.filt == "yes"] - cor[is.filt == "no"], .groups = "drop") %>%
  mutate(index = paste0(val_ref, "_", celltype)) %>%
  arrange(delta_cor) %>%
  mutate(index = factor(index, levels = unique(index)))

xCell2results_processed <- xCell2results_processed %>%
  mutate(ref = str_extract(val_ref, "(kass_blood|kass_tumor|bp|lm22|ts_blood|sc_pan_cancer)$"),
         val = str_replace(val_ref, "(kass_blood|kass_tumor|bp|lm22|ts_blood|sc_pan_cancer)$", "")) %>%
  mutate(val = str_sub(val, 1, -2)) %>%
  mutate(val_type = val2type[val])


# Calculate the median delta_cor
median_delta_cor <- median(xCell2results_processed$delta_cor)

# Create the plot
ggplot(xCell2results_processed, aes(x = index, y = delta_cor)) +
  geom_col(aes(fill = delta_cor)) +
  scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "grey", midpoint = 0) +
  geom_hline(yintercept = median_delta_cor, linetype = "dashed", color = "blue") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  theme_minimal() +
  labs(x = "", y = "Delta Correlation", title = "Waterfall Plot of Delta Correlation")






bind_rows(xCell2results) %>%
  filter(val_ref == "GSE107011_ts_blood" & topGenesFrac == 1) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor))

bind_rows(xCell2results) %>%
  filter(val_ref == "GSE107011_ts_blood" & topGenesFrac == 1 & sig_name %in% signatures_filtered) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor))

celltypes2use <- bind_rows(xCell2results) %>%
  filter(val_ref == "GSE107011_ts_blood" & topGenesFrac == 1) %>%
  pull(celltype) %>%
  unique() %>%
  gsub("#", "", .)

simulations.ct2use <- simulations[celltypes2use]



all_cors <- parallel::mclapply(celltypes2use, function(ctoi){

  signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]
  simulations.ct2use.ctoi <- simulations.ct2use[[ctoi]]

  mixNames <- unique(colnames(simulations.ct2use.ctoi))
  mixNames <- mixNames[-1]

  sigsUsed <- names(signatures_ctoi)
  xx <- list()

    for(mixname in rev(mixNames)){

    print(length(sigsUsed))
    if (length(sigsUsed) < 3) {
      return(NA)
    }

    simulations.ct2use.ctoi.sub <- simulations.ct2use.ctoi[,colnames(simulations.ct2use.ctoi) == mixname]
    signatures_ctoi.sub <- signatures_ctoi[names(signatures_ctoi) %in% sigsUsed]

    x <- parallel::mclapply(signatures_ctoi.sub, function(sig){

      mean(apply(simulations.ct2use.ctoi.sub, 2, function(m){
        cor(m[sig], gep_mat[sig, ctoi])
      }))

    }, mc.cores = 5)

    x.sorted <- sort(unlist(x), decreasing = T)
    x.sorted.sub <- x.sorted[1:round(length(x.sorted)*0.9)]
    sigsUsed <- names(x.sorted.sub)

    if (length(sigsUsed) < 3) {
      xx[[mixname]] <- x.sorted
    }else{
      xx[[mixname]] <- x.sorted.sub
    }

    }

  return(xx)

}, mc.cores = 5)
names(all_cors) <- celltypes2use

all_cors_final <- lapply(celltypes2use, function(ctoi){

  z = xCell2results %>%
    filter(val_ref == "GSE107011_ts_blood" & topGenesFrac == 1 & celltype == paste0(ctoi, "#"))

  sapply(all_cors[[ctoi]], function(s){

    z %>%
      filter(sig_name %in% names(s)) %>%
      pull(cor) %>%
      mean(.)

  })

})
names(all_cors_final) <- celltypes2use

enframe(all_cors_final, name = "celltype") %>%
  unnest_longer(value, values_to = "cor", indices_to = "mix") %>%
  mutate(mix = factor(mix, levels = rev(mixNames))) %>%
  ggplot(., aes(x=mix, y=cor, col=celltype)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


xCell2results %>%
  filter(val_ref == "BG_blood_lm22" & topGenesFrac == 1) %>%
  filter(sig_name %in% names(sigs_sub)) %>%
  group_by(celltype) %>%
  summarise(cor = mean(cor))


bind_rows(xCell2results) %>%
  mutate(topGenesFrac = as.factor(topGenesFrac)) %>%
  group_by(topGenesFrac, val_ref, celltype) %>%
  summarise(cor = median(cor)) %>%
  group_by(topGenesFrac, val_ref) %>%
  summarise(cor = median(cor)) %>%
  ggplot(., aes(x=topGenesFrac, y=cor, fill=topGenesFrac)) +
  geom_boxplot()


bind_rows(xCell2results) %>%
  mutate(ngenes = as.numeric(ngenes)) %>%  # Convert to numeric
  mutate(gene_category = cut(ngenes,
                             breaks = c(3, 5, 10, 20, 50, 100, Inf),
                             labels = c("3-5", "5-10", "10-20", "20-50", "50-100", "100+"),
                             include.lowest = TRUE)) %>%
  group_by(val_ref, gene_category, celltype) %>%
  summarise(cor = mean(cor)) %>%
  group_by(gene_category, val_ref) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=gene_category, y=cor, fill=gene_category)) +
  geom_boxplot()



















all_cors <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.sigsCors.19dec.rds")

all_cors <- bind_rows(all_cors)

all_cors %>%
  group_by(val_ref, ngenes) %>%
  summarise(cor = mean(cor)) %>%
  mutate(ngenes = as.numeric(ngenes)) %>%  # Convert to numeric
  mutate(gene_category = cut(ngenes,
                             breaks = c(3, 5, 10, 20, 50, 100, 120, 140, 160, 180, Inf),
                             labels = c("3-5", "5-10", "10-20", "20-50", "50-100", "100-120", "120-140", "140-160", "160-180", "180+"),
                             include.lowest = TRUE)) %>%
  group_by(val_ref, gene_category) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=gene_category, y=cor, fill=gene_category)) +
  geom_boxplot()

all_cors %>%
  group_by(val_ref, diff) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=diff, y=cor, fill=diff)) +
  geom_boxplot()

all_cors %>%
  group_by(val_ref, prob) %>%
  summarise(cor = mean(cor)) %>%
  ggplot(., aes(x=prob, y=cor, fill=prob)) +
  geom_boxplot()

all_cors %>%
  mutate(prob_diff = paste0(prob, "-", diff)) %>%
  group_by(val_ref, prob_diff) %>%
  summarise(cor = mean(cor)) %>%
  group_by(prob_diff) %>%
  summarise(cor = median(cor)) %>%
  arrange(-cor) %>% View()

  ggplot(., aes(x=prob_diff, y=cor, fill=prob_diff)) +
  geom_boxplot()







  library(pheatmap)


  data <- all_cors

  data <- data %>%
    filter(method %in% c("xCell2", "dtangle"))

  data <- data %>%
    mutate(refMethod = paste0(ref, "-", method)) %>%
    group_by(refMethod, val) %>%
    summarise(cors_list = list(cor),
              n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
    #            n_ct_samples = list(rep(1, length(cor)))) %>% # Rhos are weighted by number of samples per cell type
    rowwise() %>%
    mutate(cor = combineRhos(rhos = cors_list, sample_sizes = sqrt(n_ct_samples), use_median = FALSE, summarize = TRUE)) %>%
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


  data_matrix[data_matrix<0] <- 0


  # Generate the heatmap
  pheatmap(data_matrix,
           annotation_col = column_annotation,
           color = colors_pos,
           breaks = seq(0, 1, 0.02),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           labels_col = gsub(".*-", "", colnames(data_matrix)),
           gaps_col = gap_indices,
           legend_breaks = seq(-1, 1, by = 0.2))

  library(pheatmap)


  data <- all_cors

  data <- data %>%
    filter(method %in% c("xCell2", "dtangle"))

  data <- data %>%
    mutate(refMethod = paste0(ref, "-", method)) %>%
    group_by(refMethod, val) %>%
    summarise(cors_list = list(cor),
              n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
    #            n_ct_samples = list(rep(1, length(cor)))) %>% # Rhos are weighted by number of samples per cell type
    rowwise() %>%
    mutate(cor = combineRhos(rhos = cors_list, sample_sizes = sqrt(n_ct_samples), use_median = FALSE, summarize = TRUE)) %>%
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


  data_matrix[data_matrix<0] <- 0


  # Generate the heatmap
  pheatmap(data_matrix,
           annotation_col = column_annotation,
           color = colors_pos,
           breaks = seq(0, 1, 0.02),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           labels_col = gsub(".*-", "", colnames(data_matrix)),
           gaps_col = gap_indices,
           legend_breaks = seq(-1, 1, by = 0.2))

