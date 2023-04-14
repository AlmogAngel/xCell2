library(tidyverse)

# Load Kassandra's blood reference
blood_ref <- all_models_expr[,!is.na(all_models_annot$Blood_model_annot)]
blood_labels <- all_models_annot[!is.na(all_models_annot$Blood_model_annot),]
all(colnames(blood_ref) == blood_labels[!is.na(blood_labels$Blood_model_annot),]$Sample)

blood_labels <- blood_labels %>%
  mutate(label = plyr::mapvalues(Blood_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(Dataset, Sample, ont, label)

# Add lab data
laboratory_data_expressions <- read_tsv("../xCell2.0/Kassandara_data/laboratory_data_expressions.tsv")
laboratory_data_expressions <- data.frame(laboratory_data_expressions[,-1], row.names = laboratory_data_expressions$Gene, check.names = F)
laboratory_data_annotation <- read_tsv("../xCell2.0/Kassandara_data/laboratory_data_annotation.tsv")
colnames(laboratory_data_annotation)[1] <- "Sample"
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)
laboratory_data_expressions <- laboratory_data_expressions[,laboratory_data_annotation$Sample]
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)

lab_labels <- laboratory_data_annotation %>%
  mutate(label = plyr::mapvalues(Cell_type, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  mutate(Dataset = "Kassandra_lab") %>%
  dplyr::select(Dataset, Sample, ont, label)

blood_labels <- rbind(blood_labels, lab_labels)
blood_ref <- cbind(blood_ref, laboratory_data_expressions)

# -------------------------------------------------------------------------------------------------

xcell2_blood_ref

# Create simulations
getMixtures <- function(ref, labels, ct, ct_data, dep_list, max_control_type = 5,
                           mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02))){

  getControlsMeanExpression <- function(ref, labels, ct, dep_list, max_control_type, mixture_fractions){

    controls2use <- names(dep_list)[!names(dep_list) %in% c(ct, dep_list[[ct]])]

    controls <- sapply(1:length(mixture_fractions), function(x){
      labels %>%
        filter(label %in% controls2use) %>%
        group_by(Dataset) %>%
        slice_sample(n=1) %>% # One cell type per dataset
        group_by(label) %>%
        slice_sample(n=1) %>% # One sample per datasets
        ungroup() %>%
        slice_sample(n=max_control_type) %>%
        pull(Sample) %>%
        ref[,.] %>%
        as.matrix() %>%
        Rfast::rowmeans()
    })

    controls_fracs <- controls %*% diag(1-mixture_fractions)
    return(controls_fracs)

  }

  mixSmaples <- function(ref, labels, ct, dep_list, max_control_type, ct_mean_expression = x[1,]$ct_mean_expression[[1]], mixture_fractions){

    # Generate a matrix of CTOI fractions:
    m <- matrix(rep(ct_mean_expression, length(mixture_fractions)), byrow = FALSE, ncol = length(mixture_fractions)) %*% diag(mixture_fractions)
    # Get random controls fractions
    c <- getControlsMeanExpression(ref, labels, ct, dep_list, max_control_type, mixture_fractions)
    # Add controls
    m <- m + c

    rownames(m) <- rownames(ref)
    colnames(m) <- as.character(mixture_fractions)

    return(m)

  }



  ct_data %>%
    group_by(Dataset) %>%
    summarise(samples = list(Sample)) %>%
    rowwise() %>%
    mutate(ct_mean_expression = list(Rfast::rowmeans(as.matrix(ref[,samples])))) %>%
    rowwise() %>%
    mutate(mixtures = list(mixSmaples(ref, labels, ct, dep_list, max_control_type, ct_mean_expression, mixture_fractions))) %>%
    pull(mixtures) %>%
    return()



}

getSigsCor <- function(xcell2_ref, ct, mixture_ranked){

  print(ct)
  sigs <- xcell2_blood_ref@all_signatures[startsWith(names(xcell2_blood_ref@all_signatures), paste0(ct, "#"))]

  # Score every ranked mixtures of CTOI
  cors <- sapply(mixture_ranked, function(ranked_mix){

    # Score
    scores <- sapply(sigs, simplify = TRUE, function(sig){
      singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
    })
    if (is.list(scores)) {
      sigs <- sigs[-which(lengths(scores) == 0)]
      scores <- sapply(sigs, simplify = TRUE, function(sig){
        singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
      })
    }
    colnames(scores) <- names(sigs)

    # Correlation
    fracs <- as.numeric(colnames(ranked_mix))
    apply(scores, 2, function(sig_scores){
      cor(sig_scores, fracs, method = "spearman")
    })

  })

  median_cors <- Rfast::rowMedians(cors)
  names(median_cors) <- rownames(cors)

  return(median_cors)
}

simulations <- blood_labels %>%
  group_by(ont, label) %>%
  nest() %>%
  rowwise() %>%
  mutate(mixture = list(getMixtures(ref = blood_ref, labels = blood_labels, ct = label, ct_data = data, dep_list = dep_list))) %>%
  mutate(mixture_ranked = list(lapply(mixture, function(mix){singscore::rankGenes(mix)})))

simulations_final <- simulations %>%
  mutate(signatures_cor = list(getSigsCor(xcell2_ref = xcell2_blood_ref, ct = label, mixture_ranked = mixture_ranked)))

# saveRDS(simulations_final, "../xCell2.0/tmp_R_data/simulations_maxtype5.rds")

filtered_sigs <- simulations_final %>%
  mutate(sig_filtered = list(names(signatures_cor)[signatures_cor >= quantile(signatures_cor, 0.9)])) %>%
  pull(sig_filtered) %>%
  unlist()










