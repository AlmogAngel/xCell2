

clean_names = function(names) {
  names %>%
    gsub("^X", "", .) %>%
    gsub("\\.|\\-|\\s|\\+", "_", .) %>%
    gsub(",", "_", .)
}



predict_response_lightgbm <- function(data,
                                    return_shap = FALSE,
                                     num_threads = 1){

  lightgbm::setLGBMthreads(1)


  # Split to train/test
  trainIndex <- caret::createDataPartition(y = data$label, p = 0.75, list = FALSE, times = 1)

  train_set <- data[trainIndex, ]
  test_set <- data[-trainIndex, ]

  train_x <- train_set[,-ncol(train_set)]
  train_y <- as.numeric(train_set$label)-1

  test_x <- test_set[,-ncol(test_set)]
  test_y <- as.numeric(test_set$label)-1


  # >>> Tune learning rate:
  dtrain <- lightgbm::lgb.Dataset(data = data.matrix(train_x), label = train_y, categorical_feature = which(sapply(train_x, is.factor)))

  initial_learning_rates <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 3, 5)

  cv_results <- parallel::mclapply(initial_learning_rates, function(lr) {
    params <- list(objective = "binary",
                   metric = "auc",
                   learning_rate = lr)

    lightgbm::lgb.cv(params = params,
                     data = dtrain,
                     nfold = 10,
                     nrounds = 1000,
                     early_stopping_rounds = 10,
                     verbose = -1)

  }, mc.cores = num_threads)

  # Extract AUC for each learning rate
  auc_values <- sapply(cv_results, function(res) {
    max(unlist(res$record_evals$valid$auc$eval))
  })

  # Find the best learning rate
  best_initial_lr <- initial_learning_rates[which.max(auc_values)]

  # Create a sequence of learning rates around the initial best value
  refinement_factor <- 0.8   # Here we use a relative adjustment of ±80% of the initial best learning rate
  learning_rates <- unique(round(seq(from = best_initial_lr * (1 - refinement_factor),
                        to = best_initial_lr * (1 + refinement_factor),
                        length.out = 50), 2))
  learning_rates <- learning_rates[learning_rates > 0]

  cv_results <- parallel::mclapply(learning_rates, function(lr) {
    params <- list(objective = "binary",
                   metric = "auc",
                   learning_rate = lr)

    lightgbm::lgb.cv(params = params,
                     data = dtrain,
                     nfold = 10,
                     nrounds = 1000,
                     early_stopping_rounds = 10,
                     verbose = -1)
  }, mc.cores = num_threads)

  auc_values <- sapply(cv_results, function(res) {
    max(unlist(res$record_evals$valid$auc$eval))
  })

  # Find the refined best learning rate
  best_lr <- learning_rates[which.max(auc_values)]


  # >>> Tune other hyperparameters - stage 1:
  dtrain <- lightgbm::lgb.Dataset(data = data.matrix(train_x), label = train_y, categorical_feature = which(sapply(train_x, is.factor)), params = list(feature_pre_filter = FALSE))

  tune_grid_stage1 <- expand.grid(
    num_leaves = c(5, 10, 15, 31, 60),
    lambda_l1 = c(0, 1, 5, 10),
    lambda_l2 = c(0, 5, 10, 20)
  )


  tune_lgb_stage1 <- function(params) {
    lightgbm::lgb.cv(
      params = list(
        objective = "binary",
        metric = "auc",
        learning_rate = best_lr,
        num_leaves = params$num_leaves,
        lambda_l1 = params$lambda_l1,
        lambda_l2 = params$lambda_l2
      ),
      data = dtrain,
      nfold = 5,
      nrounds = 1000,
      early_stopping_rounds = 10,
      verbose = -1
    )
  }

  results_stage1 <- parallel::mclapply(1:nrow(tune_grid_stage1), function(i) {
    params <- tune_grid_stage1[i, ]
    cv_result <- tune_lgb_stage1(params)
    auc <- max(unlist(cv_result$record_evals$valid$auc$eval))
    list(params = params, auc = auc)
  }, mc.cores = num_threads)

  # Find the best hyperparameters from Stage 1
  best_result_stage1 <- results_stage1[[which.max(sapply(results_stage1, function(res) res$auc))]]
  best_params_stage1 <- best_result_stage1$params

  # # >>> Tune other hyperparameters - stage 2:
  # tune_grid_stage2 <- expand.grid(
  #   min_split_gain = c(0.5, 1, 1.5, 3),
  #   max_bin = c(50, 100, 200, 500),
  #   min_sum_hessian_in_leaf = c(1, 3, 5, 10),
  #   extra_trees = c(TRUE, FALSE),
  #   path_smooth = c(0, 0.5, 1, 10, 50),
  #   boost_from_average = c(TRUE, FALSE),
  #   early_stopping_rounds = 30
  # )
  #
  # tune_lgb_stage2 <- function(params) {
  #   lightgbm::lgb.cv(
  #     params = list(
  #       objective = "binary",
  #       metric = "auc",
  #       learning_rate = best_lr,  # Use the best learning rate found earlier
  #       num_leaves = best_params_stage1$num_leaves,
  #       lambda_l1 = best_params_stage1$lambda_l1,
  #       lambda_l2 = best_params_stage1$lambda_l2,
  #       min_split_gain = params$min_split_gain,
  #       min_sum_hessian_in_leaf = params$min_sum_hessian_in_leaf,
  #       extra_trees = params$extra_trees,
  #       path_smooth = params$path_smooth,
  #       boost_from_average = params$boost_from_average
  #     ),
  #     data = dtrain,
  #     nfold = 5,
  #     nrounds = 10000,
  #     early_stopping_rounds = params$early_stopping_rounds,
  #     verbose = -1
  #   )
  # }
  #
  # results_stage2 <- parallel::mclapply(1:nrow(tune_grid_stage2), function(i) {
  #   params <- tune_grid_stage2[i, ]
  #   dtrain = lightgbm::lgb.Dataset(data = train_x, label = train_y, params = list(max_bin = params$max_bin))
  #   cv_result <- tune_lgb_stage2(params)
  #   auc <- max(unlist(cv_result$record_evals$valid$auc$eval))
  #   list(params = params, auc = auc)
  # }, mc.cores = num_threads)
  #
  # # Find the best hyperparameters from Stage 2
  # best_result_stage2 <- results_stage2[[which.max(sapply(results_stage2, function(res) res$auc))]]
  # best_params_stage2 <- best_result_stage2$params


  # Fit best model
  # best_params <- cbind(best_params_stage1, best_params_stage2, "learning_rate" = best_lr)
  best_params <- cbind(best_params_stage1, "learning_rate" = best_lr)

  train_params <- list(
    boosting = "dart", # Use DART boosting
    drop_rate = 0.01,
    skip_drop = 0.03, # DART-specific parameter
    feature_fraction = 0.5,
    feature_pre_filter = FALSE,
    max_bin = 300
  )

  dtrain <- lightgbm::lgb.Dataset(data = data.matrix(train_x), label = train_y, categorical_feature = which(sapply(train_x, is.factor)))
  dtest <- lightgbm::lgb.Dataset(data = data.matrix(test_x), label = test_y, categorical_feature = which(sapply(test_x, is.factor)))




  lightgbm::setLGBMthreads(num_threads)
  # dtrain <- lightgbm::lgb.Dataset(data = train_x, label = train_y, params = list(max_bin = best_params_stage2$max_bin))


  best_model <- lightgbm::lgb.train(
    data = dtrain,
    params = list(
      objective = "binary",
      metric = "auc",
      learning_rate = best_params$learning_rate,
      num_leaves = best_params$num_leaves,
      lambda_l1 = best_params$lambda_l1,
      lambda_l2 = best_params$lambda_l2,

      boosting = "dart", # Use DART boosting
      drop_rate = 0.01,
      skip_drop = 0.03, # DART-specific parameter
      feature_fraction = 0.5,
      feature_pre_filter = FALSE,
      max_bin = 300,

      extra_trees=TRUE,
      force_col_wise=TRUE

      # min_data_in_leaf = best_params$min_data_in_leaf,
      # feature_fraction = best_params$feature_fraction,
      # bagging_fraction = best_params$bagging_fraction,
      # min_split_gain = best_params$min_split_gain,
      # min_sum_hessian_in_leaf = best_params$min_sum_hessian_in_leaf,
      # extra_trees = best_params$extra_trees,
      # path_smooth = best_params$path_smooth,
      # boost_from_average = best_params$boost_from_average
    ),
    valids = list(val = dtest),
    nrounds = 10000,
    verbose = -1)

  gbm_params <- list(
    metric = "auc",  #   "multi_logloss"
    objective = "binary",  #   "multiclass"
    # num_class = 4,
    learning_rate = best_lr,
    num_iterations = 2500,
    num_leaves = 400,
    max_depth = 9,
    min_data_in_leaf = 50,
    min_gain_to_split = 0.0001,
    extra_trees=TRUE,
    force_col_wise=TRUE)

  best_model <- lightgbm::lgb.train(
    data = dtrain,
    params = gbm_params,
    valids = list(val = dtest),
    nrounds = 10000,
    verbose = -1)




  if (return_shap) {
    shap_values <- shapviz::shapviz(best_model, X_pred = as.matrix(test_x))

    return(shap_values)
    # shapviz::sv_importance(shap_values, kind = "beeswarm")



  }


  test_pred <- predict(best_model, data.matrix(test_x))
  test_roc <- pROC::roc(test_y, test_pred)
  test_auc <- pROC::auc(test_roc)

  auc <- as.numeric(gsub("Area under the curve: ", "", test_auc))
  return(auc)

}



create_shap_and_suv_data <- function(dir, data, num_of_models, case) {

  lightgbm::setLGBMthreads(1)
  data <-  na.omit(data)

  models_shap_fetures_res <- c()
  models_suv_res <- c()
  for (i in 1:num_of_models)
  {
    # Split to train/test
    trainIndex = caret::createDataPartition(y = data$label, p = 0.75, list = FALSE, times = 1)
    train_set <- data[trainIndex, ]
    test_set <- data[-trainIndex, ]
    label_idx = grep("label", colnames(train_set))
    train_x = data.matrix(train_set[, -label_idx])
    train_y = train_set[, label_idx]
    test_x = data.matrix(test_set[, -label_idx])
    test_y = test_set[, label_idx]
    train_y <- ifelse(train_y == 0, 0, 1)
    test_y <- ifelse(test_y == 0, 0, 1)

    dtrain = lightgbm::lgb.Dataset(data = train_x, label = train_y)
    # lr <- choose_lr(dtrain)

    train_params <- list(
      boosting = "dart",
      drop_rate = 0.01,
      skip_drop = 0.03,
      feature_fraction = 0.5,
      feature_pre_filter = FALSE
    )

    gbm_params <- list(
      metric = "auc",
      objective = "binary",
      learning_rate = 0.2,#lr,
      num_iterations = 2500,
      num_leaves = 400,
      max_depth = 9,
      min_data_in_leaf = 50,
      min_gain_to_split = 0.0001,
      extra_trees=TRUE,
      force_col_wise=TRUE)

    dtrain <- lightgbm::lgb.Dataset(data = train_x, label = train_y, params=train_params)     # for catagorial data categorical_feature=c(ncol(data) -1,ncol(data) -2,ncol(data) -3)
    dtest <- lightgbm::lgb.Dataset(data = test_x, label = test_y)

    model <- lightgbm::lgb.train(
      data = dtrain,
      params = gbm_params,
      valids = list(val = dtest),
      verbose = -1)

    print(paste0("iter: ", i, ", score: " ,model$best_score))
    shp <- shapviz::shapviz(model, X_pred = as.matrix(test_x))
    models_shap_fetures_res[[i]] <- colMeans(abs(get_shap_values(shp))) # abs mean
    models_suv_res[[i]] <- predict(model, test_x) # for suv
  }

  saveRDS(models_shap_fetures_res, paste0(dir,"shap_grid_res_100_", case, ".rds"))
  saveRDS(models_suv_res, paste0(dir, "suv_grid_res_100_", case, ".rds"))

  #plot_shap(models_shap_fetures_res, case)
  #plot_suv(models_suv_res)
}

plot_suv <- function(res, pred_col){
  surv.in <- readRDS("/bigdata/loainaom/overallSurvival.rds")
  metadata <- readRDS("/bigdata/almogangel/xCell2/dev_scripts/loai/totaldata_mini.rds")

  res <- do.call(rbind, lapply(seq_along(models_suv_res_pd), function(i) {
    data.frame(
      Model = i,
      Sample = names(models_suv_res_pd[[i]]),
      Value = unlist(models_suv_res_pd[[i]])
    )
  }))
  mean_pred_res <- aggregate(Value ~ Sample, data = res, mean)
  rownames(mean_pred_res) <- mean_pred_res$Sample

  common_rows <- Reduce(intersect, list(rownames(surv.in), rownames(mean_pred_res)))

  surv.in <- surv.in[common_rows, ]
  mean_pred_res <- mean_pred_values_pd[common_rows, ]
  metadata <- metadata[common_rows, ]

  surv.in$pred_col <- mean_pred_res$Value
  surv.in$cancerType <- metadata$Cancer_Type

  cancer_types <- c("Melanoma","Urothelial Carcinoma","NSCLC")

  km_fits <- list()
  km_plots <-list()
  for (cancer_type in cancer_types) {
    surv_in_subset <- surv.in[surv.in$cancerType == cancer_type, ]
    quantiles <- quantile(surv_in_subset[[pred_col]], probs = 0.5)
    surv_in_subset$group <- ifelse(surv_in_subset[[pred_col]] <= quantiles[[1]], "Low 50%", "High 50%")

    # Filter out NA group
    surv.in_filtered <- surv_in_subset[!is.na(surv_in_subset$group), ]

    # Create a survival object
    surv_obj <- Surv(time = surv.in_filtered$survival_time, event = surv.in_filtered$survival_status)

    # Fit the Kaplan-Meier model
    km_fit <- survfit(surv_obj ~ group, data = surv.in_filtered)

    # Plot the Kaplan-Meier curve
    km_plot <- ggsurvplot(km_fit, data = surv.in_filtered, pval = TRUE, conf.int = TRUE,
                          xlab = "Time",
                          ylab = "Survival Probability",
                          #legend.labs = c("Low 50%", "High 50%"),
                          legend.title = "",
                          risk.table = TRUE,
                          risk.table.col = "strata",
                          title = paste0("Cancer Type: ", cancer_types[i], ", Prediction: ", pred_col),
                          ggtheme = theme_minimal())

    # Store the plot in the list
    km_plots[[cancer_type]] <- km_plot$plot
  }

  # Combine the Kaplan-Meier plots into one figure
  combined_plot <- plot_grid(plotlist = km_plots, labels = "", ncol = 3)

  # Print the combined plot
  print(combined_plot)
}

plot_shap <- function(res, case) {
  shap_values_df <- do.call(rbind, lapply(seq_along(res), function(i) {
    data.frame(
      Model = i,
      Feature = names(res[[i]]),
      Value = unlist(res[[i]])
    )
  }))

  # Calculate the median SHAP value for each feature
  feature_medians <- aggregate(Value ~ Feature, data = shap_values_df, median)

  # Select the top 10 features with the highest median SHAP values
  top10_features <- feature_medians[order(-feature_medians$Value), "Feature"][1:10]

  # Filter the data frame to include only the top 10 features
  top10_shap_values_df <- shap_values_df[shap_values_df$Feature %in% top10_features, ]

  # Create the box plot
  ggplot(top10_shap_values_df, aes(x = reorder(Feature, Value, FUN = median), y = Value)) +
    geom_boxplot() +
    coord_flip() +  # Flip coordinates for a horizontal box plot
    theme_minimal() +
    labs(title = paste0("Top 10 Features with Highest Median SHAP Values for ", case),
         x = "Feature",
         y = "mean|SHAP Value|") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_aucs_res <- function(auc_pd, auc_cr, auc_response) {

  # Flatten the list into a data frame
  methods_auc_df <- do.call(rbind, lapply(seq_along(method2use_all), function(i) {
    data.frame(Value = auc_pd[[i]], method = method2use_all[i], Type = "PD")
  }))  %>%
    bind_rows(do.call(rbind, lapply(seq_along(method2use_all), function(i) {
      data.frame(Value = auc_cr[[i]], method = method2use_all[i], Type = "CR")
    }))) %>%
    bind_rows(do.call(rbind, lapply(seq_along(method2use_all), function(i) {
      data.frame(Value = auc_response[[i]], method = method2use_all[i], Type = "Response")
    })))


  methods_auc_df <- methods_auc_df %>%
    mutate(method = factor(method, levels = method2use_all),
           is_xcell2 = ifelse(method == "xcell2", "yes", "no"))

  # Perform statistical tests for each grid group
  #xcell2_vs_cbrx <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "cbrx")),
  #                                        method = "wilcox.test",
  #                                        group.by = "grid")
  #xcell2_vs_metadata <- ggpubr::compare_means(auc ~ method, data = grid_res %>% filter(method %in% c("xcell2", "metadata")),
  #                                            method = "wilcox.test",
  #                                            group.by = "grid")

  # Plot the data with p-values
  ggplot(methods_auc_df, aes(x = method, y = Value, fill = Type, alpha=is_xcell2)) +
    geom_boxplot(width = .5, position = position_dodge(width = 0.75), show.legend = TRUE) +
    scale_fill_manual(values = c("PD" = "#A86060", "CR" = "#60A860", "Response"= "#6060A8")) +
    scale_alpha_manual(values=c("no"=0.4, "yes"=1)) +
    guides(alpha = "none") +
    theme_minimal() +
    labs(y = "AUC", x="Method") +
    scale_y_continuous(limits = c(0.58, 0.80), breaks = seq(0.6, 0.80, by = 0.05)) +
    scale_x_discrete(labels = c(
      "xcell2" = "xCell2",
      "xcell1" = "xCell1",
      "dtangle" = "dtangle",
      "cbrx" = "CIBERSORTx",
      "bayesprism" = "BayesPRISM",
      "decon" = "DeconRNASeq",
      "epic" = "EPIC",
      "mcpcounter" = "MCPcounter",
      "metadata" = "Clinical Data",
      "tide" = "TIDE",
      "impers" = "IMPERS",
      "easier" = "EASIER"
    )) +
    #  ggpubr::stat_compare_means(comparisons = list(c("xcell2", "cbrx"), c("xcell2", "dtangle")),
    #                             method = "wilcox.test",
    #                             label = "p.format",
    #                            label.y = c(0.76, 0.78)) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 20)
    )
}

plot_aucs_res_for_2 <- function(auc_pd, auc_response) {

  # Flatten the list into a data frame
  methods_auc_df <- do.call(rbind, lapply(seq_along(method2use_all), function(i) {
    data.frame(Value = auc_pd[[i]], method = method2use_all[i], Type = "PD")
  }))  %>% bind_rows(do.call(rbind, lapply(seq_along(method2use_all), function(i) {
    data.frame(Value = auc_response[[i]], method = method2use_all[i], Type = "Response")
  })))


  methods_auc_df <- methods_auc_df %>%
    mutate(method = factor(method, levels = method2use_all),
           is_xcell2 = ifelse(method == "xcell2", "yes", "no"))

  # Plot the data with p-values
  ggplot(methods_auc_df, aes(x = method, y = Value, fill = Type, alpha=is_xcell2)) +
    geom_boxplot(width = .5, position = position_dodge(width = 0.75), show.legend = TRUE) +
    scale_fill_manual(values = c("PD" = "#A86060", "Response"= "#6060A8")) +
    scale_alpha_manual(values=c("no"=0.4, "yes"=1)) +
    guides(alpha = "none") +
    theme_minimal() +
    labs(y = "AUC", x="Method") +
    scale_y_continuous(limits = c(0.58, 0.80), breaks = seq(0.6, 0.80, by = 0.05)) +
    scale_x_discrete(labels = c(
      "xcell2" = "xCell2",
      "xcell1" = "xCell1",
      "dtangle" = "dtangle",
      "cbrx" = "CIBERSORTx",
      "bayesprism" = "BayesPRISM",
      "decon" = "DeconRNASeq",
      "epic" = "EPIC",
      "mcpcounter" = "MCPcounter",
      "metadata" = "Clinical Data",
      "tide" = "TIDE",
      "impers" = "IMPERS",
      "easier" = "EASIER"
    )) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 20)
    )
}
