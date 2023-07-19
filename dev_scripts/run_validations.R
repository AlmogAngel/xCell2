library(tidyverse)
library(xCell2)
library(EPIC)
library(BayesPrism)
library(MCPcounter)
library(dtangle)
library(DeconRNASeq)
library(ggh4x)


setwd("/bigdata/almogangel/xCell2_data/benchmarking_data/")

# Load cell types labels conversion file
celltype_conversion <- read_tsv("celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

# Set references to use
refList <- list(rna_seq = c(blood = "kass_blood", tumor = "kass_tumor", mixed = "bp"),
                array = c(mixed = "lm22"),
                sc = c(blood = "ts_blood",
                       tumor = "sc_pan_cancer"))



# Function to load validation mixtures and truth values
loadVals <- function(valList,
                     valMixDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/",
                     valTruthDir = "/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/",
                     conversion_mat = celltype_conversion){

  # load mixtures
  valMixtures <- lapply(valList, function(tissue){
    mixtures <- lapply(tissue, function(val){
      as.matrix(read.table(paste0(valMixDir, val, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE))
    })
    names(mixtures) <- tissue
    mixtures
  })

  # Load truth
  valTruth <- lapply(valList, function(tissue){
    truths <- lapply(tissue, function(val){
      truth <- as.matrix(read.table(paste0(valTruthDir, val, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1))
      # Change cell types labels
      rownames(truth) <- plyr::mapvalues(rownames(truth), conversion_mat$all_labels, conversion_mat$xCell2_labels, warn_missing = FALSE)
      rows <- rownames(truth)
      truth <- suppressWarnings(apply(truth, 2, as.numeric))
      rownames(truth) <- rows
      truth[is.na(truth)] <- 0
      truth <- truth[rowSums(truth) != 0,]
      rows <- rownames(truth)
      truth[!duplicated(rows),]
    })
    names(truths) <- tissue
    truths
  })

  return(list(mixtures = valMixtures, truth = valTruth))
}



combineRefVal <- function(vals, refList, splitDependencies){


  vals.tbl <- enframe(vals$mixtures, name = "val_type", value = "mixture") %>%
    unnest_longer(mixture, indices_to = "val_dataset") %>%
    left_join(unnest_longer(enframe(vals$truth, name = "val_type", value = "truth"), truth, indices_to = "val_dataset"), by = c("val_type", "val_dataset"))

  refs.tbl <- enframe(refList, name = "ref_type", value = "ref_name") %>%
    unnest_longer(ref_name, indices_to = "ref_tissue")

  types <- unique(vals.tbl$val_type)

  # Generate a tibble with validations-references combinations
  combined.tbl <- lapply(types, function(type){
    val_tmp <- vals.tbl[vals.tbl$val_type %in% type,]
    type <- if(type == "other") "mixed" else c(type, "mixed")
    ref_tmp <- refs.tbl[refs.tbl$ref_tissue %in% type,]
    crossing(val_tmp, ref_tmp)
  }) %>%
    do.call(rbind, .) %>%
    dplyr::select(-mixture, -truth)
  rm(vals.tbl, refs.tbl)


  # Load refs dependencies
  refs_dir <- "/bigdata/almogangel/xCell2_data/dev_data/"
  getDependencies <- function(lineage_file_checked){
    ont <- readr::read_tsv(lineage_file_checked, show_col_types = FALSE) %>%
      mutate_all(as.character)

    celltypes <- pull(ont[,2])
    celltypes <- gsub("_", "-", celltypes)
    dep_list <- vector(mode = "list", length = length(celltypes))
    names(dep_list) <- celltypes

    for (i in 1:nrow(ont)) {
      descendants <-  gsub("_", "-", strsplit(pull(ont[i,3]), ";")[[1]])
      descendants <- descendants[!is.na(descendants)]

      ancestors <-  gsub("_", "-", strsplit(pull(ont[i,4]), ";")[[1]])
      ancestors <- ancestors[!is.na(ancestors)]

      dep_list[[i]] <- list("descendants" = descendants, "ancestors" = ancestors)

    }

    return(dep_list)
  }
  refsDepList <- lapply(refList, function(ref_type){
    refsDep <- lapply(ref_type, function(ref){

      # Load reference
      ref <- readRDS(paste0("references/", ref, "_ref.rds"))

      # Get cell types dependencies
      dep_list <- getDependencies(ref$lineage_file)

      dep_cts <- sapply(names(dep_list), function(d){
        deps <- unname(unlist(dep_list[[d]]))
      })
      dep_cts[order(sapply(dep_cts, length), decreasing = T)]
    })
    names(refsDep) <- ref_type
    refsDep
  })


  # Split references with dependent cell types for deconvolution methods
  if (!splitDependencies) {

    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      mutate(refsDeps = NA) %>%
      # I non deconvolution method celltype_classes = shared_celltypes
      mutate(celltype_classes = list(shared_celltypes)) %>%
      return()

  }else{

    splitDepCellTypes <- function(refsDeps){

      refsDeps <- refsDeps[order(sapply(refsDeps, length), decreasing = T)]

      classes <- list("1" = c())
      for (ct in names(refsDeps)) {

        classes <- classes[order(sapply(classes, length), decreasing = F)]
        deps <- refsDeps[[ct]]

        for (class in 1:length(classes)) {
          if (!any(deps %in% classes[[class]])) {
            classes[[class]] <- c(classes[[class]], ct)
            break
          }
        }

        if (!ct %in% unname(unlist(classes))) {
          classes[[paste0(length(classes)+1)]] <- ct
        }
      }

      one_element_list <- unname(lapply(classes, length) == 1)
      if (sum(one_element_list) > 0) {
        classes[!one_element_list][[1]] <- unname(c(unlist(classes[!one_element_list][1]), unlist(classes[one_element_list])))
        classes <- classes[!one_element_list]
      }


      return(classes)
    }

    combined.tbl %>%
      rowwise() %>%
      # Get number of samples in the validation dataset
      mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
      # Get shared cell types between the validation and reference
      mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(refsDepList[[ref_type]][[ref_name]])))) %>%
      mutate(n_shared_celltypes = length(shared_celltypes)) %>%
      # Must have at least 3 cell types
      filter(n_shared_celltypes > 2) %>%
      # Create a list of dependent cell types - filter cell types that are shared between the reference and validation
      mutate(refsDeps = list(refsDepList[[ref_type]][[ref_name]][shared_celltypes])) %>%
      # Split dependent cell types in each CIBERSORTx run
      mutate(celltype_classes = list(splitDepCellTypes(refsDeps))) %>%
      unnest(celltype_classes) %>%
      return()

  }

}

# Function to validate the all results use the same cell types
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

# Multi-levels plot function:
# level 1 - Per method/validation/reference
# level 2 - Per validation/reference
# level 3 - Per reference
# level 4 - All references
plotCorrelations <- function(res_combined, vals, level, cor_method = "spearman"){


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
        ggpubr::stat_cor(method = "spearman")



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
  getCorTests <- function(res_mat, valList, valType, valDataset, corMethod = cor_method){

    truth_mat <- valList$truth[[valType]][[valDataset]]
    celltypes <- intersect(rownames(res_mat), rownames(truth_mat))

    cor.tests <- lapply(celltypes, function(ct){

      truth <- truth_mat[ct,]

      if(any(is.na(truth))){
        print("NAs in truth!!!!")
      }

      res <- res_mat[ct,]

      samples <- intersect(names(res), names(truth))

      cor.test(res[samples], truth[samples],
               alternative = "greater",
               method = corMethod,
               exact = TRUE)
    })
    names(cor.tests) <- celltypes

    return(cor.tests)
  }

  res_combined.pvals <- res_combined %>%
    rowwise() %>%
    mutate(cor_test = list(getCorTests(res, valList = vals, valType = val_type, valDataset = val_dataset))) %>%
    unnest_longer(cor_test, indices_to = "celltype") %>%
    rowwise() %>%
    mutate(cor = cor_test$estimate,
           pvalue = cor_test$p.value) %>%
    mutate(cor = replace_na(cor, 0),
           pvalue = replace_na(pvalue, 1)) %>%
    mutate(pvalue = ifelse(pvalue == 0, .Machine$double.eps, pvalue)) %>%  # Replace zero with smallest positive number
    mutate(tras_p_value =  -log(pvalue)) %>%
    ungroup()

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


        # P-value plot
        levels <- data %>%
          group_by(method) %>%
          summarise(median_p = median(tras_p_value)) %>%
          arrange(-median_p) %>%
          pull(method)

        max_y <- round((max(data$tras_p_value)/10)+0.5)*10

        p1 <- data %>%
          mutate(method = factor(method, levels=levels)) %>%
          ggplot(., aes(x=method, y=tras_p_value)) +
          geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
          geom_jitter(aes(col=celltype), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 2) +
          scale_y_continuous(limits = c(0, max_y), breaks = seq(0,max_y,10), labels = as.character(seq(0,max_y,10))) +
          scale_fill_manual(values = methods_palette) +
          scale_color_manual(values = celltypes_palette) +
          geom_hline(yintercept = -log(0.05), color = "red", linetype = 2) +
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
          labs(x = "", y = "-log(p-value)", x = NULL, colour = "Cell Type", fill = "Method", title = plot_title)

        # Correlation plot
        levels <- data %>%
          group_by(method) %>%
          summarise(median_cor = median(cor)) %>%
          arrange(-median_cor) %>%
          pull(method)

        p2 <- data %>%
          mutate(method = factor(method, levels=levels)) %>%
          ggplot(., aes(x=method, y=cor)) +
          geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
          geom_jitter(aes(col=celltype), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 2) +
          scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
          scale_fill_manual(values = methods_palette) +
          scale_color_manual(values = celltypes_palette) +
          geom_hline(yintercept = 0, color = "red", linetype = 2) +
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
          guides(color=FALSE) +
          labs(x = "", y = "r", x = NULL, colour = "Cell Type", fill = "Method")


        # Get the legend from the second plot
        g <- ggplotGrob(p1)
        legend <- g$grobs[[which(g$layout$name == "guide-box")]]
        p1 <- p1 + guides(color=FALSE)

        # Arrange the plots and the legend
        p_final <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, p2, ncol = 1), legend, widths = c(10, 7))

        p_final <- ggplot() +
          theme_void() +
          annotation_custom(p_final)


      }else{

        # Level 3 plot

        plot_title <- refName

        colors <- c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                    "#424242", "#8B1C62", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")


        # P-value plot
        data_tmp <- data %>%
          group_by(method, val_type, val_dataset, n_val_samples) %>%
          # summarise(ref_p = -log(metap::sumlog(pvalue)$p)) # combined p-values
          summarise(ref_p = -log(median(pvalue)))


        levels <- data_tmp %>%
          group_by(method) %>%
          summarise(median_p = median(ref_p)) %>%
          arrange(-median_p) %>%
          pull(method)

        max_y <- round((max(data_tmp$ref_p)/10)+0.5)*10 +10


        p_final <- data_tmp %>%
          mutate(method = factor(method, levels=levels)) %>%
          ggplot(., aes(x=method, y=ref_p)) +
          geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
          geom_jitter(aes(col=val_dataset), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 3) +
          scale_y_continuous(limits = c(0, max_y), breaks = seq(0,max_y,20), labels = as.character(seq(0,max_y,20))) +
          scale_fill_manual(values = methods_palette) +
          scale_color_manual(values = colors) +
          geom_hline(yintercept = -log(0.05), color = "red", linetype = 2) +
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
          labs(x = "", y = "-log(combined(p-value))", x = NULL, colour = "Validation Data", fill = "Method", title = plot_title)

      }


      return(p_final)
    }



    if (level == 2) {
      # By reference-validation combination
      plots.list <- res_combined.pvals %>%
        group_by(ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        nest() %>%
        rowwise() %>%
        mutate(plots = list(getPlots(data = data, refName = ref_name, valDataset = val_dataset))) %>%
        pull(plots)

      names <- res_combined.pvals %>%
        group_by(ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        nest() %>%
        mutate(ref_val_comb = paste0(ref_name, "##", val_dataset)) %>%
        pull(ref_val_comb)
      names(plots.list) <- names

      return(plots.list)
    }

    if (level == 3) {
      # By reference
      plots.list <- res_combined.pvals %>%
        group_by(ref_tissue, ref_type, ref_name) %>%
        nest() %>%
        rowwise() %>%
        mutate(plots = list(getPlots(data = data, refName = ref_name))) %>%
        pull(plots)

      names <- res_combined.pvals %>%
        group_by(ref_tissue, ref_type, ref_name) %>%
        nest()  %>%
        pull(ref_name)
      names(plots.list) <- names

      return(plots.list)
    }


    if (level == 4) {
      # All references

      tmp <- res_combined.pvals %>%
        group_by(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples) %>%
        # summarise(ref_val_p = ifelse(n() >1, metap::sumlog(pvalue)$p, pvalue)) %>% # combined p-values
        summarise(ref_val_p = median(pvalue)) %>%
        group_by(method, ref_tissue, ref_type, ref_name) %>%
        # summarise(ref_p = ifelse(n() >1, -log(metap::sumlog(ref_val_p)$p), ref_val_p))
        summarise(ref_p = median(ref_val_p))

      colors <- c(RColorBrewer::brewer.pal(11, "Paired")[c(2,4,6,8,10)],
                  "#424242", "#8B1C62", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")


      levels <- tmp %>%
        group_by(method) %>%
        summarise(median_ref_p = median(ref_p)) %>%
        arrange(-median_ref_p) %>%
        pull(method)

      max_y <- round((max(tmp$ref_p)/10)+0.5)*10 +10

      tmp %>%
        mutate(method = factor(method, levels=levels)) %>%
        ggplot(., aes(x=method, y=ref_p)) +
        geom_boxplot(aes(fill=method), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
        geom_jitter(aes(col=ref_name), position = position_jitterdodge(jitter.width = .1, dodge.width = .5), size = 3) +
        scale_y_continuous(limits = c(0, max_y), breaks = seq(0,max_y,20), labels = as.character(seq(0,max_y,20))) +
        scale_fill_manual(values = methods_palette) +
        scale_color_manual(values = colors) +
        geom_hline(yintercept = -log(0.05), color = "red", linetype = 2) +
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
        labs(x = "", y = "(-log(combined(combined(p-value)))", x = NULL, colour = "Reference", fill = "Method", title = "All Refrences") %>%
        return()

    }

  }

}

# Method functions -----------------------------------------------


# Run xCell2
getxCell2Res <- function(ref_val_table, vals, celltype_conversion){

  runxCell2 <- function(vals, shared_celltypes, valType, valName, sigs, refName){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    res <- xCell2Analysis(mix, sigs)
    res <- res[rownames(res) %in% shared_celltypes, ]

    return(res)
  }


  # Load xCell2 signatures
  sigsList <- sapply(unique(ref_val_table$ref_name), function(ref){
    readRDS(paste0("references/xcell2_sigs/", ref, "_sigs.rds"))
  })


  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    # Get number of samples in the validation dataset
    mutate(n_val_samples = ncol(vals$truth[[val_type]][[val_dataset]])) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), unique(gsub("#.*", "", names(sigsList[[ref_name]]@signatures)))))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    filter(n_shared_celltypes > 2) %>%
    # Run xCell2
    mutate(res = list(runxCell2(vals, shared_celltypes, valType = val_type, valName = val_dataset, sigs = sigsList[[ref_name]], refName = ref_name))) %>%
    mutate(method = "xCell2", .before = everything())

  return(vals.refs.res)


}

# Run CIBERSORTx
getCIBERSORTxRes <- function(ref_val_table, vals, celltype_conversion){


  runCIBERSORTx <- function(vals, valType, valName, refName, refType, celltypes2use , dir = "/bigdata/almogangel/CIBERSORTx_docker"){

    print(paste0(valName, "_", refName))
    single_cell <- ifelse(refType == "sc", TRUE, FALSE)


    # Subset cell types from the signature matrix
    sigmat_tmp <- read.csv(paste0(dir, "/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat_tmp <- cbind("NAME" = rownames(sigmat_tmp), sigmat_tmp[,celltypes2use])
    sigmat_tmp_file <- paste0(dir, "/sigmat-tmp.txt")
    write.table(sigmat_tmp, file = sigmat_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

    if (single_cell) {
      # Subset reference file
      ref_tmp <- read.csv(paste0(dir, "/", refName, "_ref.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
      index <- colnames(ref_tmp) %in% celltypes2use
      names <- colnames(ref_tmp)[index]
      ref_tmp <- cbind(rownames(ref_tmp), ref_tmp[,index])
      colnames(ref_tmp) <- c("genes", names)
      ref_tmp_file <- paste0(dir, "/ref-tmp.txt")
      write.table(ref_tmp, file = ref_tmp_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)
    }

    # Make mixture file
    mix <- vals$mixtures[[valType]][[valName]]
    mix_tmp <- cbind("genes" = rownames(mix), mix)
    mix_file <- paste0(dir, "/mix-tmp.txt")
    write.table(mix_tmp, file = mix_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)


    # Run CIBERSORTx
    token <- "b72da36961922443b75a1b65beef27c0"

    # Make results directory
    results_dir <- paste0(dir, "/results")
    if (!dir.exists(results_dir)){
      dir.create(results_dir)
    }

    # Clean old results
    if(length(list.files(results_dir)) > 0){
      system(paste0("rm -f ", results_dir, "/*"))
    }


    if (single_cell) {
      cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                    token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell , " --rmbatchSmode ", single_cell,
                    " --refsample ", ref_tmp_file, " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
    }else{
      cmd <- paste0("docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/fractions --username almog.angel@campus.technion.ac.il --token ",
                    token, " --sigmatrix ", sigmat_tmp_file,  " --mixture ", mix_file, " --single_cell ", single_cell ," --rmbatchBmode ", !single_cell,
                    " --verbose TRUE 1> ", results_dir, "/cibersortx.stdout 2> ", results_dir, "/cibersortx.stderr")
    }



    # Run Docker via shell
    system(cmd, wait = TRUE)


    # Load results
    res_file <- ifelse("CIBERSORTx_Adjusted.txt" %in% list.files(results_dir), "CIBERSORTx_Adjusted.txt", "CIBERSORTx_Results.txt")
    cibersortx_out <- t(read.table(paste0(results_dir, "/", res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE))
    cibersortx_out <- cibersortx_out[!rownames(cibersortx_out) %in% c("P-value", "Correlation", "RMSE"),]


    return(cibersortx_out)
  }


  # Add CIBERSORT default reference (LM22)
  lm22 <- read.table("/bigdata/almogangel/xCell2_data/benchmarking_data/references/sigmats/CIBERSORT_LM22_sigmat.txt", header = T, sep = "\t", check.names = F)

  ref_val_table <- ref_val_table %>%
    dplyr::select(val_type, val_dataset, n_val_samples) %>%
    unique() %>%
    mutate(ref_type = "default",
           ref_name = "CIBERSORT_LM22",
           ref_tissue = "mixed") %>%
    rowwise() %>%
    dplyr::select(-n_val_samples, everything()) %>%
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), colnames(lm22)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least 3 cell types
    filter(n_shared_celltypes > 2) %>%
    mutate(refsDeps = NA) %>%
    mutate(celltype_classes = list(shared_celltypes)) %>%
    rbind(., ref_val_table)

# Run CIBERSORTx
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runCIBERSORTx(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "CIBERSORTx", .before = everything())


  return(vals.refs.res)

}

# Run EPIC
getEPICRes <- function(ref_val_table, vals, celltype_conversion){

  runEPIC <- function(vals, default_refs, valType, valName, refsRDSList, refName, refType, refTissue, celltypes2use, dir = "references"){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    if (refType == "default") {

      res <- t(EPIC(bulk=mix, reference=default_refs[[refTissue]], withOtherCells=FALSE, scaleExprs=FALSE)$cellFractions)
      return(res)

    }


    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigGenes <- rownames(sigmat)

    # Subset GEP
    gep <- read.csv(paste0(dir, "/gep/", refName, "_gep.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    gep <- gep[,celltypes2use]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- ref.in$ref[,celltypeIndex]
    colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label


    # Generate reference for EPIC
    ref.raw.var <- sapply(unique(colnames(ref.raw)), function(ct){
      if(sum(colnames(ref.raw) == ct) > 1){
        apply(ref.raw[,colnames(ref.raw) == ct], 1, sd)
      }else{
        rep(0, length(ref.raw[,colnames(ref.raw) == ct]))
      }
    })
    ref.raw.var <- ref.raw.var[rownames(gep), colnames(gep)]

    epic_ref <- list("sigGenes" = sigGenes,
                     "refProfiles" = as.matrix(gep),
                     "refProfiles.var" = ref.raw.var)

    # Run EPIC
    res <- t(EPIC(bulk=mix, reference=epic_ref, withOtherCells=FALSE, scaleExprs=FALSE)$cellFractions)
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


  # Add EPIC default references (BRef/TRef)
  bref <- EPIC::BRef
  tref <- EPIC::TRef
  colnames(bref$refProfiles) <- plyr::mapvalues(colnames(bref$refProfiles), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(bref$refProfiles.var) <- plyr::mapvalues(colnames(bref$refProfiles.var), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(tref$refProfiles) <- plyr::mapvalues(colnames(tref$refProfiles), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  colnames(tref$refProfiles.var) <- plyr::mapvalues(colnames(tref$refProfiles.var), from = celltype_conversion$all_labels, to = celltype_conversion$xCell2_labels, warn_missing = FALSE)
  default_refs <- list("tumor" = tref, "blood" = bref)

  ref_val_table <- ref_val_table %>%
    dplyr::select(val_type, val_dataset, n_val_samples) %>%
    unique() %>%
    mutate(ref_type = "default",
           ref_name = ifelse(val_type == "tumor", "TRef", "BRef"),
           ref_tissue = ifelse(val_type == "tumor", "tumor", "blood")) %>%
    rowwise() %>%
    dplyr::select(-n_val_samples, everything()) %>%
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), colnames(default_refs[[val_type]]$refProfiles)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    # Must have at least 3 cell types
    filter(n_shared_celltypes > 2) %>%
    mutate(refsDeps = NA) %>%
    mutate(celltype_classes = list(shared_celltypes)) %>%
    rbind(., ref_val_table)


  # Run EPIC
  vals.refs.res <- ref_val_table %>%
    rowwise() %>%
    mutate(res = list(runEPIC(vals, default_refs, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, refTissue = ref_tissue, celltypes2use = celltype_classes))) %>%
    mutate(method = "EPIC", .before = everything())


  return(vals.refs.res)


}

# Run BayesPrism
getBayesPrismRes <- function(ref_val_table, vals, celltype_conversion){

  runBayesPrism <- function(vals, valType, valName, refsRDSList, refName, refType, celltypes2use, CPUs = 30){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- as.matrix(t(ref.in$ref[,celltypeIndex]))

    type <- ifelse(refType == "sc", "count.matrix", "GEP")
    ref.filtered <- cleanup.genes(input=ref.raw,
                                  input.type=type,
                                  species="hs",
                                  gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"))


    ref.filtered.pc <- select.gene.type(ref.filtered,
                                         gene.type = "protein_coding")

    labels <- ref.in$labels[celltypeIndex,]$label

    if (refType == "sc") {
      diff.exp.stat <- get.exp.stat(sc.dat=ref.raw[,colSums(ref.raw>0)>3],# filter genes to reduce memory use
                                    cell.type.labels=labels,
                                    cell.state.labels=labels,
                                    psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                    cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                    n.cores=CPUs) #number of threads


      ref.filtered.pc <- select.marker(sc.dat=ref.filtered.pc,
                                       stat=diff.exp.stat,
                                       pval.max=0.01,
                                       lfc.min=0.1)
    }else{

      # Subset marker genes from the signature matrix
      markers <- read.csv(paste0("references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
      shared_markers <- intersect(unique(markers$marker), colnames(ref.filtered.pc))
      ref.filtered.pc <- ref.filtered.pc[,shared_markers]
    }


    tumorKey <- NULL
    if ("malignant cell" %in% labels) {
      tumorKey <- "malignant cell"
    }

    myPrism <- new.prism(
      reference=ref.filtered.pc,
      mixture=t(mix),
      input.type=type,
      cell.type.labels = labels,
      cell.state.labels =labels,
      key=tumorKey)

    bp.res <- run.prism(prism = myPrism, n.cores = CPUs)

    res <- t(get.fraction (bp=bp.res,
                           which.theta="final",
                           state.or.type="type"))

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



  vals.refs.res <- ref_val_table %>%
    # Run BayesPrism
    rowwise() %>%
    mutate(res = list(runBayesPrism(vals, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "BayesPrism", .before = everything())



  return(vals.refs.res)



}

# Run MCPcounter
getMCPcounterRes <- function(ref_val_table, vals, celltype_conversion){

  runMCPcounter <- function(vals, valType, valName, markers, celltypes2use, refName){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    markers_tmp <- markers[celltypes2use]

    markers_tmp <- enframe(markers_tmp, value = "HUGO symbols", name = "Cell population") %>%
      unnest(`HUGO symbols`) %>%
      dplyr::select(`HUGO symbols`, `Cell population`) %>%
      as.data.frame()

    res <- MCPcounter.estimate(expression = mix, featuresType = "HUGO_symbols", genes = markers_tmp)

    return(res)
  }

  # Get marker genes

  refs <- unique(ref_val_table$ref_name)
  refs_markers <- sapply(refs, function(refName){
    markers <- read.csv(paste0("references/markers/", refName, "_markers.txt"), sep = "\t", header = T, check.names = F)
    split(markers$marker, markers$label)
  })


  vals.refs.res <- ref_val_table %>%
    # Add marker genes
    mutate(markers = list(refs_markers[[ref_name]])) %>%
    # Get shared cell types between the validation and reference
    mutate(shared_celltypes = list(intersect(rownames(vals$truth[[val_type]][[val_dataset]]), names(markers)))) %>%
    mutate(n_shared_celltypes = length(shared_celltypes)) %>%
    filter(n_shared_celltypes > 2) %>%
    # Run MCPcounter
    mutate(res = list(runMCPcounter(vals, valType = val_type, valName = val_dataset, markers, celltypes2use = shared_celltypes, refName = ref_name))) %>%
    mutate(method = "MCPcounter", .before = everything())


  return(vals.refs.res)

}

# Run dtangle
getdtangleRes <- function(ref_val_table, vals, celltype_conversion){

  rundtangle <- function(vals, valType, celltypes2use, refsRDSList, valName, refName, refType){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset cell types from the raw reference matrix
    ref.in <- refsRDSList[[refType]][[refName]]
    celltypeIndex <- ref.in$labels$label %in% celltypes2use
    ref.raw <- ref.in$ref[,celltypeIndex]
    colnames(ref.raw) <- ref.in$labels[celltypeIndex,]$label

    # Normalize to CPM
    if (refType == "sc") {
      seurat_object <- Seurat::CreateSeuratObject(counts = ref.raw)
      seurat_object <- Seurat::NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
      ref.cpm <- seurat_object@assays$RNA@data
      colnames(ref.cpm) <- colnames(ref.raw)
      ref.raw <- ref.cpm
    }


    celltypes <- unique(colnames(ref.raw))
    pure_samples_list <- lapply(celltypes, function(ct){
      which(colnames(ref.raw) == ct)
    })
    names(pure_samples_list) <- celltypes

    shared_genes <- intersect(rownames(mix), rownames(ref.raw))
    mix <- mix[shared_genes,]
    ref.raw <- as.matrix(ref.raw[shared_genes,])

    if (refType == "sc") {
      y <- cbind(ref.raw, mix)
      y <- limma::normalizeBetweenArrays(y)
      ref.raw <- y[,1:ncol(ref.raw)]
      mix <- y[,(ncol(ref.raw)+1):ncol(y)]
    }


    markerMethod <- ifelse(refType == "sc", "ratio", "p.value")
    res <- dtangle(Y = t(mix), references = t(ref.raw), pure_samples = pure_samples_list, marker_method = markerMethod, data_type = "rna-seq")
    res <- t(res$estimates)

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


  vals.refs.res <- ref_val_table %>%
    # Run dtangle
    rowwise() %>%
    mutate(res = list(rundtangle(vals, valType = val_type, valName = val_dataset, refsRDSList, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "dtangle", .before = everything())




  return(vals.refs.res)

}

# Run DeconRNASeq
getDeconRNASeqRes <- function(ref_val_table, vals, celltype_conversion){

  runDeconRNASeq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]

    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,celltypes2use]

    # Run DeconRNASeq
    res <- DeconRNASeq(data.frame(mix), sigmat, checksig=FALSE, known.prop = FALSE,  fig = FALSE)
    res <- as.matrix(res$out.all)
    rownames(res) <- colnames(mix)
    res <- t(res)

    return(res)
  }


  vals.refs.res <- ref_val_table %>%
    # Run DeconRNASeq
    rowwise() %>%
    mutate(res = list(runDeconRNASeq(vals, valType = val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "DeconRNASeq", .before = everything())


  return(vals.refs.res)

}



# Run quanTIseq
# Workaround from: https://github.com/icbi-lab/quanTIseq/tree/master/quantiseq/deconvolution
# as used in "Twelve Years of Cellular Deconvolution: Applications, Benchmark, Methodology, and Challenges"
# TODO: check this workaround
getquanTIseqRes <- function(ref_val_table, vals, celltype_conversion){

  source("/bigdata/almogangel/xCell2/dev_scripts/quantiseq_code.R")
  runquanTIseq <- function(vals, valType, valName, refName, refType, celltypes2use, dir = "references"){

    print(paste0(valName, "_", refName))
    mix <- vals$mixtures[[valType]][[valName]]


    # Subset sigGenes from the signature matrix
    sigmat <- read.csv(paste0(dir, "/sigmats/", refName, "_sigmat.txt"), sep = "\t", header = T, check.names = F, row.names = 1)
    sigmat <- sigmat[,celltypes2use]

    # Run quanTIseq
    res <- t(quanTIseq(sigmat, mix, scaling=rep(1, ncol(sigmat)), method="lsei"))
    res <- res[rownames(res) != "Other",]


    return(res)
  }


  vals.refs.res <- ref_val_table %>%
    # Run quanTIseq
    rowwise() %>%
    mutate(res = list(runquanTIseq(vals, valType =  val_type, valName = val_dataset, refName = ref_name, refType = ref_type, celltypes2use = celltype_classes))) %>%
    mutate(method = "quanTIseq", .before = everything())


  return(vals.refs.res)

}


# -------------- Run cytometry/other validations --------------

print("Running Cytometry/Other Validations...")


cyto.vals.list <- list(blood = c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "GSE53655", "GSE60424"),
                      tumor = c("ccRCC_cytof_CD45+", "NSCLC_cytof", "GSE121127", "WU_ccRCC_RCCTC"),
                      other = c("GSE120444", "GSE115823"))

cyto.vals <- loadVals(cyto.vals.list)

refval.tbl.nodeps <- combineRefVal(vals = cyto.vals, refList, splitDependencies = TRUE)
refval.tbl <- combineRefVal(vals = cyto.vals, refList, splitDependencies = FALSE)


print("Running CIBERSORTx...")
cbrx.cyto.res <- getCIBERSORTxRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(cbrx.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.cyto.res.rds")

print("Running EPIC...")
epic.cyto.res <- getEPICRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(epic.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.cyto.res.rds")

print("Running BayesPrism...")
bp.cyto.res <- getBayesPrismRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(bp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/bp.cyto.res.rds")

print("Running MCPcounter")
mcp.cyto.res <- getMCPcounterRes(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(mcp.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mcp.cyto.res.rds")

print("Running dtangle")
dtan.cyto.res <- getdtangleRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(dtan.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.cyto.res.rds")

print("Running DeconRNASeq")
decon.cyto.res <- getDeconRNASeqRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(decon.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/decon.cyto.res.rds")

print("Running quanTIseq")
quanti.cyto.res <- getquanTIseqRes(ref_val_table = refval.tbl.nodeps, vals = cyto.vals, celltype_conversion)
saveRDS(quanti.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/quanti.cyto.res.rds")

print("Running xCell2...")
xcell2.cyto.res <- getxCell2Res(ref_val_table = refval.tbl, vals = cyto.vals, celltype_conversion)
saveRDS(xcell2.cyto.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.cyto.res.rds")



# Combine and plot results -----------------

# Load results
cyto.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".cyto.res.rds", full.names = TRUE), function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    # merge all splitted results of deconvolution method
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)


# Validate shared cell types and remove default results
cyto.Res <- validateSharedCTs(res_combined = cyto.Res, remove_default_refs = TRUE)


# Save plots to PDF
cyto.level1.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 1, cor_method = "spearman")
cyto.level2.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 2, cor_method = "spearman")
cyto.level3.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 3, cor_method = "spearman")
cyto.level4.plot <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 4, cor_method = "spearman")

all_plots <- c(list("allRefs" = cyto.level4.plot), cyto.level3.plots, cyto.level2.plots, cyto.level1.plots)

pdf("xCell2_cytoVal_plots.pdf")
for (p in all_plots) {
  print(p)
}
dev.off()


# # Default reference analysis
#
# # Load results
# cyto.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".cyto.res.rds", full.names = TRUE), function(f){
#   readRDS(f) %>%
#     dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
#     # merge all splitted results of deconvolution method
#     group_by(across(-ncol(.))) %>%
#     summarise(res = list(do.call(rbind, res)), .groups = 'drop')
# }) %>%
#   do.call(rbind, .)
#
#
# # Validate shared cell types and remove default results
# cyto.Res <- validateSharedCTs(res_combined = cyto.Res, remove_default_refs = FALSE)
#
# default.cyto.level1.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 1, cor_method = "spearman")
# default.cyto.level2.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 2, cor_method = "spearman")
# default.cyto.level3.plots <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 3, cor_method = "spearman")
# default.cyto.level4.plot <- plotCorrelations(res_combined = cyto.Res, vals = cyto.vals, level = 4, cor_method = "spearman")
#





# -------------- Run single-cell validations --------------

print("Running Single-cell Validations...")


sc.vals.list <- list(blood = c("sc_pbmc"),
                     tumor = c("SC_glioblastoma", "SC_GSE84133", "SC_HNSCC", "SC_lymphomas", "SC_melanoma", "SC_NSCLC"))

sc.vals <- loadVals(sc.vals.list)

sc.refval.tbl.nodeps <- combineRefVal(vals = sc.vals, refList, splitDependencies = TRUE)
sc.refval.tbl <- combineRefVal(vals = sc.vals, refList, splitDependencies = FALSE)


print("Running CIBERSORTx...")
cbrx.sc.res <- getCIBERSORTxRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(cbrx.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/cbrx.sc.res.rds")

print("Running EPIC...")
epic.sc.res <- getEPICRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(epic.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/epic.sc.res.rds")

print("Running BayesPrism...")
bp.sc.res <- getBayesPrismRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(bp.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/bp.sc.res.rds")

print("Running MCPcounter")
mcp.sc.res <- getMCPcounterRes(ref_val_table = sc.refval.tbl, vals = sc.vals, celltype_conversion)
saveRDS(mcp.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/mcp.sc.res.rds")

print("Running dtangle")
dtan.sc.res <- getdtangleRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(dtan.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/dtan.sc.res.rds")

print("Running DeconRNASeq")
decon.sc.res <- getDeconRNASeqRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(decon.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/decon.sc.res.rds")

print("Running quanTIseq")
quanti.sc.res <- getquanTIseqRes(ref_val_table = sc.refval.tbl.nodeps, vals = sc.vals, celltype_conversion)
saveRDS(quanti.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/quanti.sc.res.rds")

print("Running xCell2...")
xcell2.sc.res <- getxCell2Res(ref_val_table = sc.refval.tbl, vals = sc.vals, celltype_conversion)
saveRDS(xcell2.sc.res, "/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.sc.res.rds")



# Combine and plot results -----------------


# Load results
sc.Res <- lapply(list.files("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/", pattern = ".sc.res.rds", full.names = TRUE), function(f){
  readRDS(f) %>%
    dplyr::select(method, ref_tissue, ref_type, ref_name, val_type, val_dataset, n_val_samples, n_shared_celltypes, res) %>%
    # merge all splitted results of deconvolution method
    group_by(across(-ncol(.))) %>%
    summarise(res = list(do.call(rbind, res)), .groups = 'drop')
}) %>%
  do.call(rbind, .)


# Validate shared cell types and remove default results
sc.Res <- validateSharedCTs(sc.Res, remove_default_refs = TRUE)


# Save plots to PDF
sc.level1.plots <- plotCorrelations(res_combined = sc.Res, vals = sc.vals, level = 1, cor_method = "spearman")
sc.level2.plots <- plotCorrelations(res_combined = sc.Res, vals = sc.vals, level = 2, cor_method = "spearman")
sc.level3.plots <- plotCorrelations(res_combined = sc.Res, vals = sc.vals, level = 3, cor_method = "spearman")
sc.level4.plot <- plotCorrelations(res_combined = sc.Res, vals = sc.vals, level = 4, cor_method = "spearman")

all_plots <- c(list("allRefs" = sc.level4.plot), sc.level3.plots, sc.level2.plots, sc.level1.plots)

pdf("xCell2_scVal_plots.pdf")
for (p in all_plots) {
  print(p)
}
dev.off()



