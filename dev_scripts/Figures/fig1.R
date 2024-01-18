############################################################
#
############################################################

library(tidyverse)
library(xCell2)
library(parallel)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(grid)



setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
# sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
# sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))


# Load references matrices
refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
refsRDSList <- lapply(refList, function(ref_type){
  refs <- lapply(ref_type, function(ref){
    # Load reference
    ref.in <- readRDS(paste0("references/", ref, "_ref.rds"))
    ref.in
  })
  names(refs) <- ref_type
  refs
})

vals.refs.res <- refval.tbl %>%
  # Get number of samples in the validation dataset
  mutate(n_val_samples = ncol(cyto.vals$truth[[val_type]][[val_dataset]])) %>%
  filter(n_shared_celltypes > 2) %>%
  mutate(method = "xCell2", .before = everything())

# vals.refs.res.sc <- sc.refval.tbl %>%
#   # Get number of samples in the validation dataset
#   mutate(n_val_samples = ncol(sc.vals$truth[[val_type]][[val_dataset]])) %>%
#   filter(n_shared_celltypes > 2) %>%
#   mutate(method = "xCell2", .before = everything())

# Load other methods results
files <- list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations", pattern = ".cyto.res.rds", full.names = TRUE)
files <- files[!files %in% c("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.fig1.cyto.res.rds")]
cyto.Res <- lapply(files, function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)


# xCell 2.0 settings
set.seed(123)
thisseed <- 123
cores2use <- 60
# signature settings
# topGenesFrac = 1
load_sigs <- FALSE
sigs_suffix <- "18jan"
scores_results <- FALSE
minpbcells <- 30
minpbgroups <- 10
# simulations settings
# sim_method <- c("ref_multi", "ref_thin", "ref_mix_thin")
sim_method <- "mix_thin"
samplesFrac = 0.1
# xCell2Analysis
tranform <- TRUE
spillover <- FALSE
nSims <- 10



xCell2results <- lapply(1:nrow(vals.refs.res), function(i){

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
  valType <- vals.refs.res[i,]$val_type
  valDataset <- vals.refs.res[i,]$val_dataset
  refName <- vals.refs.res[i,]$ref_name
  valDataType <- vals.refs.res[i,]$val_data_type[[1]]


  res_mat_file <-  paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_resXGB_", sigs_suffix, ".rds")
  if (file.exists(res_mat_file)) {
    res_mat <-  readRDS(res_mat_file)
    print(paste0("reading ", val_ref, " skipping..."))
    return(res_mat)
  }



    # Load signatures?
  sig_file <-  paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_topGF", topGenesFrac, "_", sigs_suffix, ".rds")
  if (load_sigs) {
    sigsFile <- sig_file
  }else{
    sigsFile <- NULL
  }



  # Load filtering data
  if (valType == "tumor") {
    filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/tumor_rnaseq_filtering_data.rds")
  }else{
    if (valDataType == "rnaseq") {
      filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/blood_rnaseq_filtering_data.rds")
    }else{
      filtering_data <- readRDS("/bigdata/almogangel/xCell2/data/blood_array_filtering_data.rds")
    }
  }

  # Remove current validation from filtering data
  filt_datasets <- gsub(x = colnames(filtering_data$mixture), pattern = "#.*", replacement = "")
  filtering_data$mixture <- filtering_data$mixture[,filt_datasets != valDataset]
  filtering_data$truth <- filtering_data$truth[,filt_datasets != valDataset]
  filtering_data$truth <- filtering_data$truth[!apply(filtering_data$truth, 1, function(x){all(is.na(x))}),]



    sigs <-  xCell2::xCell2Train(ref = ref.in, labels = labels, mix = mix.in, ref_type, filtering_data = filtering_data, ref_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                                 return_sigs_filt = TRUE, sigsFile = sigsFile, minPBcells = minpbcells, minPBsamples = minpbgroups, seed = thisseed, samples_frac = samplesFrac, top_genes_frac = topGenesFrac,
                                 nCores = cores2use, simMethod = sim_method, mix = mix, sim_fracs = simFracs, ct_sims = nSims)


    # Return signatures correlations
    if (TRUE) {
      mix_ranked <- singscore::rankGenes(mix)
      truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
      sigs_celltypes <- unique(gsub("#.*", "", names(sigs)))
      celltypes <- intersect(sigs_celltypes, rownames(truth_mat))
      all_ctoi_res <- lapply(celltypes, function(ctoi){
        signatures_ctoi <- sigs[startsWith(names(sigs), paste0(ctoi, "#"))]
        scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
          suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
        })
        fracs <- truth_mat[ctoi,colnames(mix_ranked)]
        c <- apply(scores, 2, function(x){
          cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
        })
        tibble(sig_name = names(c), n_samples = sum(!is.na(fracs)), cor = c) %>%
          separate(sig_name, into = c("celltype", "prob", "diff", "ngenes"), sep = "_", remove = FALSE) %>%
          return(.)
      })
      names(all_ctoi_res) <- celltypes

      bind_rows(all_ctoi_res) %>%
        mutate(val_ref = val_ref,
               topGenesFrac = topGenesFrac)
    }


  })

  topGenesFrac.list %>%
    bind_rows() %>%
    return(.)

  #   xcell2object_file <-  paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_xCell2object_", sigs_suffix, ".rds")
  #   saveRDS(sigs, xcell2object_file)
  #
  #   if (!load_sigs) {
  #     saveRDS(sigs, sig_file)
  #     print(paste0(val_ref, " sigs are ready."))
  #     return(NA)
  #   }
  #
  #
  # # xCell2Analysis
  # res_mat <- xCell2::xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = 0.5)
  # res_mat <- res_mat[vals.refs.res[i,]$shared_celltypes[[1]], ]
  #
  # print(paste0(val_ref, " res mat ready."))
  #
  # return(res_mat)
})

# saveRDS(xCell2results, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.fig1.cyto.res.7dec.rds")
xCell2results <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.fig1.cyto.res.9dec.rds")





allLevel1plots <- list()
allLevel2plots <- list()
for (i in 1:nrow(vals.refs.res)) {

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
  valType <- vals.refs.res[i,]$val_type
  valDataset <- vals.refs.res[i,]$val_dataset
  refName <- vals.refs.res[i,]$ref_name


  # Load other methods results
  cyto.Res.tmp <- cyto.Res %>%
    filter(ref_name == refName, val_dataset == valDataset & method != "xCell2")

  xcell2.cyto.Res.tmp <- cyto.Res.tmp[1,]
  xcell2.cyto.Res.tmp$method <- "xCell2"

  xcell2.cyto.Res.tmp$res[[1]]  <- xCell2results[[i]]
  cyto.Res.tmp <- rbind(cyto.Res.tmp, xcell2.cyto.Res.tmp)



  # Calculate correlation with ground truth
  truth_mat <- cyto.vals$truth[[valType]][[valDataset]]


  level1plots <- lapply(1:nrow(cyto.Res.tmp), function(j){

    res <- cyto.Res.tmp[j,]$res[[1]]
    celltypes <- intersect(rownames(res), rownames(truth_mat))
    method <- cyto.Res.tmp[j,]$method[[1]]

    df <- lapply(celltypes, function(ct){

      truth <- truth_mat[ct,]
      res <- res[ct,]

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

    df <- df %>%
      left_join(cor_results, by = "celltype")

    # Plot correlations - level 1
    df$truth <- df$truth*100
    df$prediction <- df$prediction*100

    p <- df %>%
      drop_na() %>%
      ggplot(., aes(x = truth, y = prediction, color = celltype)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      facet_wrap(~celltype, scales = "free", strip.position = "top") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = NA, colour = "black", size = 1),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "gray95"),
        legend.position = "none"
      ) +
      labs(title = method, subtitle = paste0("Ref: ", refName, ", Val: ", valDataset),
           x = "Truth (%)", y = "prediction (%)") +
      geom_text(data = cor_results, aes(label = label, x = -Inf, y = Inf),
                hjust = 0, vjust = 1.2, size = 4, color = "black")

    p
  })
  names(level1plots) <- cyto.Res.tmp$method

  if (!refName[[1]] %in% names(allLevel1plots)) {
    allLevel1plots[[refName]] <- list()
  }

  if (!valDataset %in% names(allLevel1plots[[refName]])) {
    allLevel1plots[[refName]][[valDataset]] <- list()
  }

  allLevel1plots[[refName]][[valDataset]] <- level1plots


  getCors <- function(res, truth = truth_mat){

    celltypes <- intersect(rownames(res), rownames(truth))

    df <- lapply(celltypes, function(ct){

      truth <- truth[ct,]
      res <- res[ct,]

      samples <- intersect(names(res), names(truth))

      tibble(celltype = ct, truth = truth[samples], prediction = res[samples])

    }) %>%
      bind_rows()

    df %>%
      group_by(celltype) %>%
      summarize(
        cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
        p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
        n = sum(!is.na(truth) & !is.na(prediction))
      ) %>%
      mutate(cor = ifelse(is.na(cor), 0, cor))

  }


  cyto.Res.tmp <- cyto.Res.tmp %>%
    rowwise() %>%
    mutate(cors = list(getCors(res)))


  df <- cyto.Res.tmp %>%
    select(method, cors) %>%
    unnest(cols = c(cors))

  celltypes_per_method <- df %>%
    group_by(method) %>%
    summarize(celltypes = list(celltype), .groups = 'keep')

  shared_celltypes <- Reduce(intersect, celltypes_per_method$celltypes)

  # Filter the original data to keep only rows with shared cell types
  df <- df %>%
    filter(celltype %in% shared_celltypes)


  # Reordering the factor levels based on median of correlation for each method
  df <- df %>%
    mutate(method = factor(method, levels = names(sort(tapply(cor, method, median, na.rm = TRUE), decreasing = TRUE))))

  # Adding a column to specify color based on method
  df <- df %>%
    mutate(isxcell = ifelse(method == "xCell2", "tomato", "gray"))

  mincor <- min(df$cor)-0.1

  # Defining distinct colors for cell types
  num_colors <- length(unique(df$celltype))
  #color_palette <- brewer.pal(num_colors, "Set2")
  color_palette <- Polychrome::createPalette(num_colors,  c("#ff0000", "#00ff00", "#0000ff"))

  # Plotting
  p2 <- ggplot(df, aes(x = method, y = cor)) +
    geom_boxplot(aes(fill = isxcell), outlier.shape = NA) +
    geom_point(aes(color = celltype), position = position_dodge(width = 0.75), size = 3) +
    scale_fill_manual(values = c("gray", "tomato")) +
    scale_color_manual(values = as.character(color_palette)) +
    theme_minimal() +
    labs(x = "", y = "Spearman Rho", color = "") +
    scale_y_continuous(limits = c(mincor, 1), breaks = seq(-1, 1, 0.2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = ifelse(levels(df$method) == "xCell2", "bold", "plain"))) +
    theme(legend.position = "right") +
    theme(legend.key.width = unit(1, "cm")) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10, face = "bold")) +
    guides(fill = FALSE)

  # Reordering the factor levels based on median of correlation for each celltype
  df <- df %>%
    mutate(celltype = factor(celltype, levels = names(sort(tapply(cor, celltype, median, na.rm = TRUE), decreasing = FALSE))))

  mincor <- min(df$cor)-0.1
  color_palette[levels(df$method) == "xCell2"] <- "red"

  p3 <- ggplot(df, aes(x = celltype, y = cor)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +
    geom_point(aes(color = method, shape = isxcell), position = position_dodge(width = 0.75), size = 3) +
    scale_color_manual(values = as.character(color_palette)) +
    scale_shape_manual(values = c("gray" = 16, "tomato" = 8)) +  # Assigns shapes based on isxcell, using backticks for TRUE and FALSE
    theme_minimal() +
    labs(x = "", y = "Spearman Rho", color = "") +
    scale_y_continuous(limits = c(mincor, 1), breaks = seq(-1, 1, 0.2)) +
    theme(axis.text.x = element_text(hjust = 1)) +
    theme(legend.position = "right") +
    theme(legend.key.width = unit(1, "cm")) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10, face = "bold")) +
    guides(shape = FALSE) +
    coord_flip()

  refvalname <- paste0("Ref: ", refName, ", Val: ", valDataset)

  p23 <- gridExtra::grid.arrange(p2, p3, ncol = 1,
                                 top = textGrob(refvalname,
                                                gp = gpar(fontface = "bold", fontsize = 14)))



  if (!refvalname %in% names(allLevel2plots)) {
    allLevel2plots[[refvalname]] <- list
  }
  allLevel2plots[[refvalname]] <- p23
}


# Save level 1 plots
pdf("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/level1plots8dec.pdf",
    width = 12, height = 16)  # Adjust 'width' and 'height' as needed
plot.new()

for (ref in names(allLevel1plots)) {

  # Create a title page
  title <- ref
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
  text(x = 0.5, y = 0.5, labels = title, cex = 5, col = "black", font = 2)
  grid.newpage()  # Move to a new page after the title

  for(val in names(allLevel1plots[[ref]])){

    # Create a title page
    title <- val
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
    text(x = 0.5, y = 0.5, labels = title, cex = 5, col = "black")
    grid.newpage()  # Move to a new page after the title


    for (method in allLevel1plots[[ref]][[val]]) {
      grid.draw(method)
      # Add a new page after drawing each plot
      grid.newpage()
    }

  }

}
dev.off()


# Save level 2 plots
pdf("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/level2plots8dec.pdf",
    width = 12, height = 16)  # Adjust 'width' and 'height' as needed
for (plot in allLevel2plots) {
  grid.draw(plot)
  # Add a new page after drawing each plot
  grid.newpage()
}
dev.off()


# Level 3 plot

# This function calculate the mean Rho value using Fisher's Z transformation
combineRhos <- function(rhos, sample_sizes = NULL, use_median = TRUE, summarize = FALSE){

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

  #weights <- sample_sizes/max(sample_sizes)


  rhos[rhos == 1] <- 0.999999999
  rhos[rhos == -1] <- -0.999999999

  # Fisher's Z Transformation
  z_values <- 0.5 * log((1 + rhos) / (1 - rhos))

  if (!summarize) {
    z_weighted <- sample_sizes * z_values / mean(sample_sizes)
    rho_weighted <- (exp(2 * z_weighted) - 1) / (exp(2 * z_weighted) + 1)

    return(rho_weighted)
  }

  if (use_median) {
    weighted_median <- function(values, weights) {
      if (length(values) != length(weights)) {
        stop("values and weights must have the same length")
      }

      # Sort values and weights by values
      order_index <- order(values)
      values <- values[order_index]
      weights <- weights[order_index]

      # Calculate the cumulative sum of weights
      cum_weights <- cumsum(weights)
      total_weight <- sum(weights)

      # Find the index where the cumulative weight exceeds half of the total weight
      median_index <- which(cum_weights >= total_weight / 2)[1]

      # Return the corresponding value
      return(values[median_index])
    }
    z_median <- weighted_median(z_values, sample_sizes)
    rho_weighted_median <- (exp(2 * z_median) - 1) / (exp(2 * z_median) + 1)

    return(rho_weighted_median)
  }

  # Variance of Z values
  # var_z <- 1.06 / (sample_sizes - 3)
  # var_z <- (1 / (sample_sizes - 3)) + (0.57722 / (2 * (sample_sizes - 1)^2))
  # var_z <- 1 / (sample_sizes - 3)

  # Weighted Average of Z values
  # weights <- 1 / var_z
  weights <- sample_sizes
  z_mean <- sum(weights * z_values) / sum(weights)

  # Back Transformation
  rho_weighted_mean <- (exp(2 * z_mean) - 1) / (exp(2 * z_mean) + 1)
  return(rho_weighted_mean)

}



all_cors <- mclapply(1:nrow(vals.refs.res), function(i){

  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
  valType <- vals.refs.res[i,]$val_type
  valDataset <- vals.refs.res[i,]$val_dataset
  refName <- vals.refs.res[i,]$ref_name


  # # Load other methods results
  cyto.Res.tmp <- cyto.Res %>%
    filter(ref_name == refName, val_dataset == valDataset & method != "xCell2")

  if (nrow(cyto.Res.tmp) != 7) {
    stop("!!!!")
  }

  xcell2.cyto.Res.tmp <- cyto.Res.tmp[1,]
  xcell2.cyto.Res.tmp$method <- "xCell2"



  res_mat_file <-  paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_resXGB_29dec.rds")
  if (file.exists(res_mat_file)) {
    res_mat <- readRDS(res_mat_file)
  }else{
    return(NA)
  }


  xcell2.cyto.Res.tmp$res[[1]]  <- res_mat
  cyto.Res.tmp <- rbind(cyto.Res.tmp, xcell2.cyto.Res.tmp)


  # Calculate correlation with ground truth
  truth_mat <- cyto.vals$truth[[valType]][[valDataset]]


  shared_celltypes <- Reduce(intersect, lapply(cyto.Res.tmp$res, rownames))
  shared_samples <- Reduce(intersect, lapply(cyto.Res.tmp$res, colnames))


  getCors <- function(res, truth = truth_mat, shared_cts = shared_celltypes, shared_samp = shared_samples){

    celltypes <- intersect(rownames(res), rownames(truth))
    celltypes <- celltypes[celltypes %in% shared_celltypes]

    samples <- intersect(colnames(res), colnames(truth))
    samples <- samples[samples %in% shared_samp]


    df <- lapply(celltypes, function(ct){

      truth <- truth[ct,samples]
      res <- res[ct,samples]


      tibble(celltype = ct, truth = truth, prediction = res)

    }) %>%
      bind_rows()

    df %>%
      group_by(celltype) %>%
      dplyr::summarize(
        cor = cor(truth, prediction, method = "spearman", use = "pairwise.complete.obs"),
        p_value = cor.test(truth, prediction, method = "spearman", exact = FALSE)$p.value,
        n = sum(!is.na(truth) & !is.na(prediction))
      ) %>%
      mutate(cor = ifelse(is.na(cor), 0, cor))

  }


  df <- cyto.Res.tmp %>%
    rowwise() %>%
    mutate(cors = list(getCors(res))) %>%
    dplyr::select(method, cors) %>%
    unnest(cols = c(cors)) %>%
    mutate(ref = refName,
           val = valDataset)


  celltypes_per_method <- df %>%
    group_by(method) %>%
    dplyr::summarize(celltypes = list(celltype), .groups = 'keep')

  shared_celltypes <- Reduce(intersect, celltypes_per_method$celltypes)

  # Filter the original data to keep only rows with shared cell types
  df <- df %>%
    filter(celltype %in% shared_celltypes)


  return(df)
}, mc.cores = 5) %>%
  bind_rows()

# all_cors <- all_cors %>%
#   filter(n >= 15)

# all_cors <- all_cors[!is.na(all_cors)] %>%
#   bind_rows()

# microarray_vals <- c("SDY311", "SDY420", "SDY67", "GSE64385", "GSE65133", "GSE106898", "GSE64655", "GSE59654", "GSE107990", "GSE20300", "GSE65135", "GSE77343", "GSE77344", "GSE93722")
# all_cors <- all_cors %>%
#   filter(!val %in% microarray_vals)


tmp <- all_cors %>%
  #filter(n<20) %>%
  filter(!val %in% c("GSE65133")) %>%
  group_by(method, ref) %>%
  dplyr::summarise(cors_list = list(cor),
            n_ct_samples = list(n)) %>% # Rhos are weighted by number of samples per cell type
#            n_ct_samples = list(rep(1, length(cor)))) %>% # Rhos are weighted by number of samples per cell type
  rowwise() %>%
  #mutate(ref_rho = combineRhos(rhos = cors_list, sample_sizes = sqrt(n_ct_samples), use_median = FALSE, summarize = TRUE),
  mutate(ref_rho = list(combineRhos(rhos = cors_list, sample_sizes = sqrt(n_ct_samples), use_median = FALSE, summarize = FALSE)),
         n_val_cts = length(cors_list))
# %>%
#   group_by(method, ref) %>%
#   summarise(cors_list2 = list(ref_val_rho),
#             n_val_cts_combied = list(n_val_cts)) %>% # Rhos are weighted by number of cell types per validation dataset
# #            n_val_cts_combied = list(rep(1, length(n_val_cts)))) %>% # Rhos are weighted by number of cell types per validation dataset
#   rowwise() %>%
#   mutate(ref_rho = combineRhos(rhos = cors_list2, sample_sizes = n_val_cts_combied))


# Custom labeling function
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

if (min(tmp$ref_rho) < 0) {
  minRho <- round(min(tmp$ref_rho)-0.1, 1)
}else{
  minRho <- 0
}

tmp %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
  ggplot(., aes(x = tidytext::reorder_within(method, ref_rho, ref, sep = "#"), y = ref_rho)) +
  geom_bar(aes(fill=is_xcell2), stat="identity") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(limits = c(minRho, 1), breaks = seq(minRho, 1, by = 0.1)) +
  coord_flip() +
  facet_wrap(~ ref, scales = "free", ncol = 1) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels = bold_xCell2_labels) +
  labs(x="", y="Mean Spearman Rho") +
  guides(fill=FALSE)

# per ref
ds.info <- read_tsv("/bigdata/almogangel/xCell2_data/datasets.txt")

tmp <- all_cors %>%
  left_join(ds.info[,c(2,6,7)], by = c("val" = "Name")) %>%
  #filter(n >= 10) %>%
  #group_by(method, TruthMethod, DataType) %>%
  group_by(method, val) %>%
  summarise(cors_list = list(cor),
            n_val_samples = list(n)) %>%
  rowwise() %>%
  mutate(val_rho = combineRhos(rhos = cors_list, sample_sizes = n_val_samples))

# Custom labeling function
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

if (min(tmp$val_rho) < 0) {
  minRho <- round(min(tmp$val_rho)-0.1, 1)
}else{
  minRho <- 0
}

all_vals <- unique(tmp$val)
tmp %>%
  #mutate(TruthMethodType = paste0(TruthMethod, "-", DataType)) %>%
  filter(val %in% all_vals[1:8]) %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
  ggplot(., aes(x = method, y = val_rho)) +
  geom_bar(aes(fill=is_xcell2), stat="identity") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(limits = c(minRho, 1), breaks = seq(minRho, 1, by = 0.1)) +
  coord_flip() +
  facet_wrap(~ val, scales = "free", ncol = 1) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels = bold_xCell2_labels) +
  labs(x="", y="Mean Spearman Rho", title = "BG_blood") +
  guides(fill=FALSE)

p2 <- tmp %>%
  filter(val %in% all_vals[9:17]) %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
  ggplot(., aes(x = tidytext::reorder_within(method, val_rho, val, sep = "#"), y = val_rho)) +
  geom_bar(aes(fill=is_xcell2), stat="identity") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(limits = c(minRho, 1), breaks = seq(minRho, 1, by = 0.1)) +
  coord_flip() +
  facet_wrap(~ val, scales = "free", ncol = 1) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels = bold_xCell2_labels) +
  labs(x="", y="Mean Spearman Rho") +
  guides(fill=FALSE)


p3 <- tmp %>%
  filter(val %in% all_vals[18:length(all_vals)]) %>%
  mutate(is_xcell2 = ifelse(startsWith(method, "xCell2"), "yes", "no")) %>%
  ggplot(., aes(x = tidytext::reorder_within(method, val_rho, val, sep = "#"), y = val_rho)) +
  geom_bar(aes(fill=is_xcell2), stat="identity") +
  scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
  theme_minimal() +
  scale_y_continuous(limits = c(minRho, 1), breaks = seq(minRho, 1, by = 0.1)) +
  coord_flip() +
  facet_wrap(~ val, scales = "free", ncol = 1) +
  tidytext::scale_x_reordered() +
  scale_x_discrete(labels = bold_xCell2_labels) +
  labs(x="", y="Mean Spearman Rho") +
  guides(fill=FALSE)


p1 + p2 +p3
