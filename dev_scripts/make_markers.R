library(tidyverse)

dir <- "/bigdata/almogangel/xCell2_data/benchmarking_data/references/"


# ref_files <- list.files(dir, pattern = "*_ref.rds")
ref_files <- list("bulk" = c("bp_ref.rds", "kass_blood_ref.rds", "kass_tumor_ref.rds", "lm22_ref.rds"), "sc" = c("ts_blood_ref.rds", "sc_pan_cancer_ref.rds"))

# ref_files <- list("sc" = c("ts_blood_ref.rds", "sc_pan_cancer_ref.rds"))
# ref_files <- list("bulk" = c("bp_ref.rds"))


# Make markers from bulk references using method from dtangle/ABIS  ------------
# Credit: https://github.com/gjhunt/dtangle
find_markers <- function(Y, references = NULL, pure_samples = NULL, data_type = NULL,
                         gamma = NULL, marker_method = "ratio") {


  #' Row-binds \code{Y} with \code{references} and generates \code{pure_samples}.
  #' @inheritParams dtangle
  combine_Y_refs <- function(Y, references, pure_samples) {

    if (is.null(pure_samples)) {
      pure_samples <- lapply(1:nrow(references), identity)
      names(pure_samples) <- rownames(references)
    }

    if (is.null(colnames(Y)) & !is.null(colnames(references))) {
      colnames(Y) <- colnames(references)
    }

    if (!is.null(references) & is.null(colnames(references)))
      colnames(references) <- colnames(Y)

    if (!is.null(references))
      Y <- as.matrix(rbind(as.matrix(references), as.matrix(Y)))

    if (is.null(colnames(Y)))
      colnames(Y) <- 1:ncol(Y)

    return(list(Y = Y, pure_samples = pure_samples))
  }

  cmbd <- combine_Y_refs(Y, references, pure_samples)
  Y <- cmbd$Y
  pure_samples <- cmbd$pure_samples

  if (any(lengths(pure_samples) == 1) & marker_method == "p.value") {
    message("Can't use p.value method.")
    marker_method <- "ratio"
  }
  if (is.null(gamma))
    gamma <- get_gamma(data_type)
  K <- length(pure_samples)
  N <- dim(Y)[2]
  pure <- unlist(pure_samples)
  C <- array(0, c(K, N))
  colnames(C) <- colnames(Y)
  if (marker_method == "ratio") {
    avg_exp_fn <- function(x) {
      colMeans(2^(Y[x, , drop = FALSE]))/gamma
    }
    eta_hats <- t(sapply(pure_samples, avg_exp_fn))
    C <- t(sapply(1:K, function(i) {
      eta_hats[i, ]/apply(eta_hats[-i, , drop = FALSE], 2, sum)
    }))
  } else if (marker_method == "regression") {
    for (i in 1:K) {
      X <- as.numeric(pure %in% pure_samples[[i]])
      Yp <- as.matrix(Y[pure, ])
      m <- stats::lm(Yp ~ 1 + X)
      cfdf <- data.frame(t(stats::coef(m)))
      C[i, ] <- cfdf$X
    }
  } else if (marker_method == "diff") {
    for (i in 1:K) {
      C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, stats::median)
    }
    less_second <- function(x) {
      x - sort(x, decreasing = TRUE)[2]
    }
    C <- apply(C, 2, less_second)
  } else if (marker_method == "p.value") {
    for (i in 1:K) {
      C[i, ] <- apply(Y[pure_samples[[i]], , drop = FALSE], 2, mean)
    }
    calc_pvals <- function(i) {
      x <- C[, i]
      top <- which(x == max(x))[1]
      second <- order(x, decreasing = TRUE)[2]
      pvs <- rep(NA, length(x))
      for (j in 1:length(x)) {
        pvs[j] <- tryCatch({
          x1 <- Y[pure_samples[[j]], i]
          x2 <- Y[pure_samples[[second]], i]
          n1 <- length(x1)
          n2 <- length(x2)
          sd1 <- stats::sd(x1)
          sd2 <- stats::sd(x2)
          sp <- sqrt(((n1 - 1) * sd1 + (n2 - 1) * sd2)/(n1 + n2 - 2))
          t.value <- (mean(x1) - mean(x2))/(sp * sqrt((1/n1) + (1/n2)))
          tmp <- stats::pt(abs(t.value), df = n1 + n2 - 2)
          tmp
        })
      }
      pvs[-top] <- 0
      return(pvs)
    }
    C <- sapply(1:ncol(C), calc_pvals)
  } else {
    stop("Marker method not found.")
  }
  pick_top <- function(x) {
    m <- which(x == max(x, na.rm = TRUE))
    if (length(m) != 1)
      return(c(NA, NaN))
    return(c(m, x[m]))
  }
  M <- apply(C, 2, pick_top)
  M <- data.frame(t(M))
  colnames(M) <- c("top", "value")
  M$rn <- 1:N
  rownames(M) <- colnames(Y)
  M$Cell.Type <- names(pure_samples)[M$top]
  if (marker_method == "p.value") {
    diffmm <- find_markers(Y = Y, pure_samples = pure_samples, data_type = data_type,
                           gamma = gamma, marker_method = "diff")$M
    M$diff <- diffmm$value
    iM <- M[stats::complete.cases(M), ]
    sM <- iM[order(iM$top, -iM$value, -iM$diff), ]
  } else {
    iM <- M[stats::complete.cases(M), ]
    sM <- iM[order(iM$top, -iM$value), ]
  }
  L <- lapply(1:K, function(i) {
    vals <- sM[sM$top == i, "rn"]
    names(vals) <- rownames(sM[sM$top == i, ])
    return(vals)
  })
  V <- lapply(1:K, function(i) {
    vals <- sM[sM$top == i, "value"]
    names(vals) <- rownames(sM[sM$top == i, ])
    return(vals)
  })
  names(L) <- names(pure_samples)
  names(V) <- names(pure_samples)
  return(list(L = L, V = V, M = M, sM = sM))
}

lapply(names(ref_files), function(ref_type){
  markerMethod <- ifelse(ref_type == "bulk", "p.value", "ratio")

  lapply(ref_files[[ref_type]], function(ref_file){

    ref.name <- gsub("_ref.rds", "", ref_file)
    print(ref.name)


    ref.in <- readRDS(paste0(dir, ref_file))
    ref <- ref.in$ref
    colnames(ref) <- ref.in$labels$label
    ref <- as.matrix(ref)

    celltypes <- unique(colnames(ref))
    pure_samples_list <- sapply(celltypes, function(ct){
      which(colnames(ref) == ct)
    })


    if (any(lengths(pure_samples_list) == 1)) {
      markerMethod <- "ratio"

      # # Remove cell types with a single sample -  cannot use p.value method
      # single_sample_celltypes <- names(which(lengths(pure_samples_list) == 1))
      # ref <- ref[,!colnames(ref) %in% single_sample_celltypes]
      # celltypes <- unique(colnames(ref))
      # pure_samples_list <- sapply(celltypes, function(ct){
      #   which(colnames(ref) == ct)
      # })
    }

    gma <- ifelse(ref.name == "lm22", dtangle:::gma$ma_gene, dtangle:::gma$rna_seq)

    markers.out <- find_markers(Y = t(ref), pure_samples = pure_samples_list, marker_method = markerMethod, gamma = gma)
    markers.df <- data.frame(label = markers.out$M$Cell.Type, marker = rownames(markers.out$M))
    write.table(markers.df,
                paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", ref.name, "_markers.txt"),
                row.names = F, quote = F, sep = "\t", col.names = T)
  })



})




# # Make markers from bulk references using method from BayesPrism: ------------
# # Credit: https://github.com/Danko-Lab/BayesPrism
# library(BayesPrism)
# lapply(ref_files_sc, function(ref_file){
#
#   ref.name <- gsub("_ref.rds", "", ref_file)
#   print(ref.name)
#
#
#   ref.in <- readRDS(paste0(dir, ref_file))
#   ref.raw <- t(as.matrix(ref.in$ref))
#
#   ref.filtered <- cleanup.genes(input=ref.raw,
#                                 input.type="count.matrix",
#                                 species="hs",
#                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"))
#
#
#   ref.filtered.pc <-  select.gene.type(ref.filtered,
#                                        gene.type = "protein_coding")
#
#   labels <- ref.in$labels$label
#
#   diff.exp.stat <- get.exp.stat(sc.dat=ref.raw[,colSums(ref.raw>0)>3],# filter genes to reduce memory use
#                                   cell.type.labels=labels,
#                                   cell.state.labels=labels,
#                                   psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
#                                   cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
#                                   n.cores=30) #number of threads
#
#
#   markers.tbl <- enframe(diff.exp.stat, name = "label") %>%
#     unnest_longer(value, indices_to = "marker") %>%
#     group_by(label)
#
#
#   markers.df <- markers.tbl %>%
#     filter(value$pval.up.min <= 0.01 & value$min.lfc >= 0.1) %>%
#     select(label, marker) %>%
#     as.data.frame()
#
#   write.table(markers.df,
#               paste0("/bigdata/almogangel/xCell2_data/benchmarking_data/references/markers/", ref.name, "_markers.txt"),
#               row.names = F, quote = F, sep = "\t", col.names = T)
# })


