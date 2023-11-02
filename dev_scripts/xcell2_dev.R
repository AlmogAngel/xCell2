library(tidyverse)
library(xCell2)

setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")


# Read reference-validation pairs

# "/bigdata/almogangel/xCell2/dev_scripts/prep_ref_val_pairs.R"
refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto_ref_val.rds")
sc.refval.tbl <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc_ref_val.rds")
# Load validation data
cyto.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/cyto.vals.rds")
sc.vals <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/ref_val_pairs/sc.vals.rds")

refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood", tumor = "sc_pan_cancer"))



getxCell2Res <- function(ref_val_table, vals, makeSigs = FALSE, useTransformation = TRUE, useSpillOver = FALSE, Params = list(),
                         saveSigs = FALSE, spillAlpha = 0, minPBc = 30, minPBg = 10, wGenes = TRUE, sigsSuff, seed2use){


  runxCell2 <- function(vals, refsRDSList, shared_celltypes, valType, valName, refType, refName, make_sigs = makeSigs, tranform = useTransformation, spillover = useSpillOver,
                        params2use = Params, save_sigs = saveSigs, spill_alpha = spillAlpha, minpbcells = minPBc, minpbgroups = minPBg, weight_genes = wGenes, sigs_suffix = sigsSuff, thisseed = seed2use){

    val_ref <- paste0(valName, "_", refName)
    print(val_ref)

    mix.in <- vals$mixtures[[valType]][[valName]]
    ref.in <- refsRDSList[[refType]][[refName]]$ref
    labels <- refsRDSList[[refType]][[refName]]$labels
    lineage_file <- refsRDSList[[refType]][[refName]]$lineage_file

    if (refType == "sc") {
      shared_cleaned_genes <- xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = TRUE, use_protein_coding = TRUE, n_var_genes = 5000)
      ref <- shared_cleaned_genes$ref
      mix <- shared_cleaned_genes$mix
    }else{
      shared_cleaned_genes <- xCell2CleanGenes(ref = ref.in, mix = mix.in, top_var_genes = FALSE, use_protein_coding = FALSE)
      ref <- shared_cleaned_genes$ref
      mix <- shared_cleaned_genes$mix
    }

    refType <- ifelse(refType == "rna_seq", "rnaseq", refType)


    if (make_sigs) {
      if (save_sigs) {


        sig_file <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
        if (!file.exists(sig_file)) {
          sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = TRUE,
                              sigsFile = NULL, rf_params = list(), minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes)
          saveRDS(sigs, sig_file)
          print("sigs done.")
        }else{
          sigs <- readRDS(sig_file)
          print("sigs loaded.")
        }

        mix_ranked <- singscore::rankGenes(mix)

        sigs.scores <- sapply(sigs, function(sig){
          singscore::simpleScore(mix_ranked, upSet = sig[[1]], centerScore = FALSE)$TotalScore
        })
        sigs.scores <- t(sigs.scores)
        colnames(sigs.scores) <- colnames(mix_ranked)

        return(sigs.scores)
      }else{
        sigsFile <- NULL
      }
    }else{
      sigsFile <- paste0("/bigdata/almogangel/xCell2_data/dev_data/sigs/", val_ref, "_sigs_", sigs_suffix, ".rds")
    }

    sigs <- xCell2Train(ref = ref, labels = labels, data_type = refType, lineage_file = lineage_file, return_sigs = FALSE,
                        sigsFile = sigsFile, rf_params = params2use, minPBcells = minpbcells, minPBsamples = minpbgroups, weightGenes = weight_genes, seed = thisseed)
    res <- xCell2Analysis(mix = mix, xcell2sigs = sigs, tranform = tranform, spillover = spillover, spillover_alpha = spill_alpha)
    res <- res[shared_celltypes, ]

    return(res)
  }


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

  # valType=x[1,]$val_type[[1]]; valName=x[1,]$val_dataset[[1]]; refType=x[1,]$ref_type[[1]]; refName=x[1,]$ref_name[[1]]; shared_celltypes=x[1,]$shared_celltypes[[1]]

  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
    filter(n_shared_celltypes > 2) %>%
    # Run xCell2
    mutate(res = list(runxCell2(vals, refsRDSList, shared_celltypes, valType = val_type, valName = val_dataset, refType = ref_type, refName = ref_name))) %>%
    mutate(method = "xCell2", .before = everything())

  return(vals.refs.res)


}

# Function to validate the all results use the same cell types:
validateSharedCTs <- function(res_combined, remove_default_refs = TRUE){


  if(remove_default_refs){
    res_combined <- res_combined %>%
      filter(ref_type != "default")
  }else{
    default.methods <- res_combined %>%
      filter(ref_type == "default") %>%
      pull(method) %>%
      unique()

    res_combined <- res_combined %>%
      filter(method %in% default.methods)
  }

  res_combined %>%
    # Group by reference and validation data
    group_by(ref_tissue, ref_type, ref_name, val_type, val_dataset) %>%
    # For each group, get a list of the intersection of cell type names across all methods
    mutate(common_cell_types = list(
      purrr::reduce(map(res, ~rownames(.x)), intersect)
    )) %>%
    # Update the 'res' column by keeping only rows (cell types) that are in 'common_cell_types'
    mutate(res = map2(res, common_cell_types, ~.x[rownames(.x) %in% .y, ])) %>%
    ungroup() %>%
    return()

}

plotCorrelations <- function(res_combined, vals, level, cor_method = "spearman", isSigsScore, useMedian = FALSE){



  # This function calculate the mean Rho value using Fisher's Z transformation
  combineRhos <- function(rhos, sample_sizes = NULL, use_median = useMedian){

    if (length(rhos) == 1) {
      return(rhos)
    }

    if (length(sample_sizes) != length(rhos)) {
      sample_sizes <- rep(sample_sizes, length(rhos))
    }

    rhos[rhos == 1] <- 0.999999999
    rhos[rhos == -1] <- -0.999999999

    # Fisher's Z Transformation
    z_values <- 0.5 * log((1 + rhos) / (1 - rhos))

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
    var_z <- 1.06 / (sample_sizes - 3)

    # Weighted Average of Z values
    weights <- 1 / var_z
    z_mean <- sum(weights * z_values) / sum(weights)

    # Back Transformation
    rho_weighted_mean <- (exp(2 * z_mean) - 1) / (exp(2 * z_mean) + 1)
    return(rho_weighted_mean)

  }


  # Level 1 plots
  if (level == 1) {
    getPlots <- function(res_mat, method, refName, valList, valType, valDataset){

      ref_val_combo <- paste0(method, ": ", "Ref - ", refName, ", Val - ", valDataset)

      truth_mat <- valList$truth[[valType]][[valDataset]]
      celltypes <- intersect(rownames(res_mat), rownames(truth_mat))

      df <- lapply(celltypes, function(ct){

        truth <- truth_mat[ct,]
        res <- res_mat[ct,]

        samples <- intersect(names(res), names(truth))
        tibble("celltype" = ct, "prediction" = res[samples], "truth" = truth[samples])

      }) %>%
        do.call(rbind, .)


      p <- ggpubr::ggscatter(df, x = "truth", y = "prediction",
                             add = "reg.line",
                             add.params = list(color = "blue", fill = "lightgray"),
                             conf.int = TRUE,
                             title = ref_val_combo)


      p2 <- ggpubr::facet(p, facet.by = "celltype", scales = "free") +
        ggpubr::stat_cor(method = "spearman", alternative = "greater")



      return(p2)
    }

    plots.list <- res_combined %>%
      rowwise() %>%
      mutate(plots = list(getPlots(res_mat = res, method = method, refName = ref_name,
                                   valList = vals, valType = val_type, valDataset = val_dataset))) %>%
      pull(plots)

    names <- res_combined %>%
      mutate(method_ref_val_comb = paste0(method, "##", ref_name, "##", val_dataset)) %>%
      pull(method_ref_val_comb)
    names(plots.list) <- names

    return(plots.list)
  }


  # Calculate p-values for levels 2,3 and 4
  # i = 35
  # valList=vals; valType=x[i,]$val_type[[1]]; valDataset=x[i,]$val_dataset[[1]]; valDataset=x[i,]$val_dataset[[1]]; res_mat=x[i,]$res[[1]]; corMethod=cor_method
  getCors <- function(res_mat, valList, valType, valDataset, corMethod = cor_method, is_sigs_score = isSigsScore){

    truth_mat <- valList$truth[[valType]][[valDataset]]

    if (is_sigs_score) {
      res_celltypes <- gsub("#.*", "", rownames(res_mat))
      celltypes <- intersect(res_celltypes, rownames(truth_mat))
      res_mat <- sapply(celltypes, function(ctoi){
        colMeans(res_mat[res_celltypes == ctoi,])
      })
      res_mat <- t(res_mat)
    }

    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))

    cor.tests <- lapply(celltypes, function(ct){

      truth <- truth_mat[ct,]

      if(any(is.na(truth))){
        stop("NAs in truth!!!!")
      }

      res <- res_mat[ct,]

      samples <- intersect(names(res), names(truth))

      cor.res <- suppressWarnings(cor(res[samples], truth[samples], method = "spearman"))
      cor.res <- ifelse(is.na(cor.res), 0, cor.res)
      cor.res

    })
    names(cor.tests) <- celltypes

    return(cor.tests)
  }

  res_combined.cors <- res_combined %>%
    rowwise() %>%
    mutate(cors = list(getCors(res, valList = vals, valType = val_type, valDataset = val_dataset))) %>%
    unnest_longer(cors, indices_to = "celltype")


  methods_palette <- c("BayesPrism" = "#424242", "CIBERSORTx" = "#424242", "DeconRNASeq" = "#424242",
                       "dtangle" = "#424242", "EPIC" = "#424242", "MCPcounter" = "#424242",
                       "quanTIseq" = "#424242", "xCell2" = "#8B0000")


  # Level 2,3 and 4 plots
  if(level != 1){


    getPlots <- function(data, refName, valDataset = NULL){

      if (!is.null(valDataset)) {

        # Level 2 plot:

        plot_title <- paste0("Ref: ", refName, ", Val: ", valDataset)

        # Cell types colors
        n_celltypes <- length(unique(data$celltype))
        celltypes_palette <- if(n_celltypes > 12){
          scales::hue_pal()(n_celltypes)
        }else{
          c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
            "#424242", "#8B1C62", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")
        }


        # Rhos plot
        levels <- data %>%
          group_by(method) %>%
          summarise(median_rho = median(cors)) %>%
          arrange(-median_rho) %>%
          pull(method)


        p1 <- data %>%
          mutate(method = factor(method, levels=levels)) %>%
          ggplot(., aes(x=method, y=cors)) +
          geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
          geom_jitter(aes(col=celltype), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 2) +
          scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +
          scale_fill_manual(values = methods_palette) +
          scale_color_manual(values = celltypes_palette) +
          theme_linedraw() +
          theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
                panel.grid.major.x = element_blank(),
                panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
                panel.grid.minor = element_line(colour = "white"),
                axis.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.key = element_rect(fill = NA),
                legend.background = element_rect(fill = NA),) +
          guides(fill=FALSE) +
          labs(x = "", y = "Spearman Rho", x = NULL, colour = "Cell Type", fill = "Method", title = plot_title)

        return(p1)

        # # Correlation plot
        # levels <- data %>%
        #   group_by(method) %>%
        #   summarise(median_cor = median(cor)) %>%
        #   arrange(-median_cor) %>%
        #   pull(method)
        #
        # p2 <- data %>%
        #   mutate(method = factor(method, levels=levels)) %>%
        #   ggplot(., aes(x=method, y=cor)) +
        #   geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
        #   geom_jitter(aes(col=celltype), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 2) +
        #   scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
        #   scale_fill_manual(values = methods_palette) +
        #   scale_color_manual(values = celltypes_palette) +
        #   geom_hline(yintercept = 0, color = "red", linetype = 2) +
        #   theme_linedraw() +
        #   theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        #         panel.grid.major.x = element_blank(),
        #         panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        #         panel.grid.minor = element_line(colour = "white"),
        #         axis.title = element_text(size = 16, face = "bold"),
        #         axis.text.x = element_text(size = 12),
        #         axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        #         legend.title = element_text(size = 12, face = "bold"),
        #         legend.text = element_text(size = 10),
        #         legend.key = element_rect(fill = NA),
        #         legend.background = element_rect(fill = NA),) +
        #   guides(fill=FALSE) +
        #   guides(color=FALSE) +
        #   labs(x = "", y = "r", x = NULL, colour = "Cell Type", fill = "Method")
        #
        #
        # # Get the legend from the second plot
        # g <- ggplotGrob(p1)
        # legend <- g$grobs[[which(g$layout$name == "guide-box")]]
        # p1 <- p1 + guides(color=FALSE)
        #
        # # Arrange the plots and the legend
        # p_final <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, p2, ncol = 1), legend, widths = c(10, 7))
        #
        # p_final <- ggplot() +
        #   theme_void() +
        #   annotation_custom(p_final)


      }else{

        # Level 3 plot

        plot_title <- refName

        colors <- c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                    "#424242", "#8B1C62", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")


        # P-value plot
        data_tmp <- data %>%
          group_by(method, val_type, val_dataset, n_val_samples) %>%
          summarise(cors_list = list(cors)) %>%
          rowwise() %>%
          mutate(ref_rho = combineRhos(rhos = cors_list, sample_sizes = n_val_samples))


        levels <- data_tmp %>%
          group_by(method) %>%
          summarise(median_rho = median(ref_rho)) %>%
          arrange(-median_rho) %>%
          pull(method)


        p_final <- data_tmp %>%
          mutate(method = factor(method, levels=levels)) %>%
          ggplot(., aes(x=method, y=ref_rho)) +
          geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
          geom_jitter(aes(col=val_dataset), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 3) +
          scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +
          scale_fill_manual(values = methods_palette) +
          scale_color_manual(values = colors) +
          theme_linedraw() +
          theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
                panel.grid.major.x = element_blank(),
                panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
                panel.grid.minor = element_line(colour = "white"),
                axis.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10),
                legend.key = element_rect(fill = NA),
                legend.background = element_rect(fill = NA),) +
          guides(fill=FALSE) +
          labs(x = "", y = "Mean Spearman Rho", x = NULL, colour = "Validation Data", fill = "Method", title = plot_title)

      }


      return(p_final)
    }



    if (level == 2) {
      # By reference-validation combination
      plots.list <- res_combined.cors %>%
        group_by(ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        nest() %>%
        rowwise() %>%
        mutate(plots = list(getPlots(data = data, refName = ref_name, valDataset = val_dataset))) %>%
        pull(plots)

      names <- res_combined.cors %>%
        group_by(ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        nest() %>%
        mutate(ref_val_comb = paste0(ref_name, "##", val_dataset)) %>%
        pull(ref_val_comb)
      names(plots.list) <- names

      return(plots.list)
    }

    if (level == 3) {
      # By reference
      plots.list <- res_combined.cors %>%
        group_by(ref_tissue, ref_type, ref_name) %>%
        nest() %>%
        rowwise() %>%
        mutate(plots = list(getPlots(data = data, refName = ref_name))) %>%
        pull(plots)

      names <- res_combined.cors %>%
        group_by(ref_tissue, ref_type, ref_name) %>%
        nest()  %>%
        pull(ref_name)
      names(plots.list) <- names

      return(plots.list)
    }


    if (level == 4) {
      # All references


      tmp <- res_combined.cors %>%
        select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, cors) %>%
        group_by(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        summarise(cors_list = list(cors)) %>%
        rowwise() %>%
        mutate(ref_val_rho = combineRhos(rhos = cors_list, sample_sizes = n_val_samples)) %>%
        group_by(method, ref_tissue, ref_type, ref_name) %>%
        summarise(cors_list2 = list(ref_val_rho),
                  sample_sizes = list(n_val_samples)) %>%
        rowwise() %>%
        mutate(ref_rho = combineRhos(rhos = cors_list2, sample_sizes = sample_sizes))


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
        ggplot(., aes(x = tidytext::reorder_within(method, ref_rho, ref_name, sep = "#"), y = ref_rho)) +
        geom_bar(aes(fill=is_xcell2), stat="identity") +
        scale_fill_manual(values = c("yes"="tomato", "no"="gray")) +
        theme_minimal() +
        scale_y_continuous(limits = c(minRho, 1), breaks = seq(minRho, 1, by = 0.1)) +
        coord_flip() +
        facet_wrap(~ ref_name, scales = "free", ncol = 1) +
        tidytext::scale_x_reordered() +
        scale_x_discrete(labels = bold_xCell2_labels) +
        labs(x="", y="Mean Spearman Rho") +
        guides(fill=FALSE) %>%
        return()


    }

  }

}


# (1) Signatures weighting analysis


# Make signatures
print("cyto no weights")
xcell2.cyto.ps30c10g.sigsres.noweight <- getxCell2Res(ref_val_table = refval.tbl, vals = cyto.vals, makeSigs = TRUE, saveSigs = TRUE, minPBc = 30, minPBg = 10, wGenes = FALSE, sigsSuff = "noweight")
saveRDS(xcell2.cyto.ps30c10g.sigsres.noweight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.sigsres.noweight.rds")
print("cyto weights")
xcell2.cyto.ps30c10g.sigsres.weight <- getxCell2Res(ref_val_table = refval.tbl, vals = cyto.vals, makeSigs = TRUE, saveSigs = TRUE, minPBc = 30, minPBg = 10, wGenes = TRUE, sigsSuff = "weight")
saveRDS(xcell2.cyto.ps30c10g.sigsres.weight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.sigsres.weight.rds")
print("sc no weights")
xcell2.sc.ps30c10g.sigsres.noweight <- getxCell2Res(ref_val_table = sc.refval.tbl, vals = sc.vals, makeSigs = TRUE, saveSigs = TRUE, minPBc = 30, minPBg = 10, wGenes = FALSE, sigsSuff = "noweight")
saveRDS(xcell2.sc.ps30c10g.sigsres.noweight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.sigsres.noweight.rds")
print("sc weights")
xcell2.sc.ps30c10g.sigsres.weight <- getxCell2Res(ref_val_table = sc.refval.tbl, vals = sc.vals, makeSigs = TRUE, saveSigs = TRUE, minPBc = 30, minPBg = 10, wGenes = TRUE, sigsSuff = "weight")
saveRDS(xcell2.sc.ps30c10g.sigsres.weight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.sigsres.weight.rds")


# Load signatures
x1 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.sigsres.noweight.rds") %>%
  mutate(method = "xCell2 No Weights (mean score)")

x2 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.sigsres.weight.rds") %>%
  mutate(method = "xCell2 Weights (mean score)")

cyto.Res <- x1 %>%
  bind_rows(x2, .)

plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 4, cor_method = "spearman", isSigsScore = TRUE)




# (2) Random Forest analysis

# Make signatures
print("cyto no weights")
xcell2.cyto.ps30c10g.noweight <- getxCell2Res(ref_val_table = refval.tbl, vals = cyto.vals, makeSigs = FALSE, saveSigs = FALSE, useTransformation = TRUE, sigsSuff = "noweight")
saveRDS(xcell2.cyto.ps30c10g.noweight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.nospill.noweight.rds")
print("cyto weights")
xcell2.cyto.ps30c10g.weight <- getxCell2Res(ref_val_table = refval.tbl, vals = cyto.vals, makeSigs = FALSE, saveSigs = FALSE, useTransformation = TRUE, sigsSuff = "weight")
saveRDS(xcell2.cyto.ps30c10g.weight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.nospill.weight.rds")
print("sc no weights")
xcell2.sc.ps30c10g.noweight <- getxCell2Res(ref_val_table = sc.refval.tbl, vals = sc.vals, makeSigs = FALSE, saveSigs = FALSE, useTransformation = TRUE, sigsSuff = "noweight")
saveRDS(xcell2.sc.ps30c10g.noweight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.nospill.noweight.rds")
print("sc weights")
xcell2.sc.ps30c10g.weight <- getxCell2Res(ref_val_table = sc.refval.tbl, vals = sc.vals, makeSigs = FALSE, saveSigs = FALSE, useTransformation = TRUE, sigsSuff = "weight")
saveRDS(xcell2.sc.ps30c10g.weight, "/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.nospill.weight.rds")


# Load results cytometry
cyto.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".cyto.res.rds", full.names = TRUE), function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)

x1 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.nospill.noweight.rds") %>%
  mutate(method = "xCell2 No Weights")

x2 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.cyto.ps30c10g.nospill.weight.rds") %>%
  mutate(method = "xCell2 Weights")


cyto.Res <- cyto.Res %>%
  filter(method != "xCell2") %>%
  bind_rows(x1, .) %>%
  bind_rows(x2, .)

cyto.Res <- validateSharedCTs(res_combined = cyto.Res, remove_default_refs = TRUE)

plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 4, cor_method = "spearman", isSigsScore = FALSE, useMedian = TRUE)


# Load results single-cell

sc.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".sc.res.rds", full.names = TRUE), function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)


x1 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.nospill.noweight.rds") %>%
  mutate(method = "xCell2 No Weights (mean score)")

x2 <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2.sc.ps30c10g.nospill.weight.rds") %>%
  mutate(method = "xCell2 Weights (mean score)")


sc.Res <- sc.Res %>%
  filter(method != "xCell2") %>%
  bind_rows(x1, .) %>%
  bind_rows(x2, .)

sc.Res <- validateSharedCTs(res_combined = sc.Res, remove_default_refs = TRUE)

plotCorrelations(res_combined = sc.Res, vals = sc.vals, level = 4, cor_method = "spearman", isSigsScore = FALSE)




# (3) Random Forest tuning

