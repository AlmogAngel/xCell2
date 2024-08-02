

clean_names = function(names) {
  names %>%
    gsub("^X", "", .) %>%
    gsub("\\.|\\-|\\s|\\+", "_", .) %>%
    gsub(",", "_", .)
}

dummy_metadata <- function(metadata, cancers, label, remove_sd) {

  if (remove_sd) {
    metadata <- metadata[metadata$RECIST != "SD",]
  }

  metadata <- metadata[metadata$Cancer_Type %in% cancers,]

  label_column <- ifelse(label %in% c("PD", "SD", "PR", "CR"), "RECIST", "Response")
  columns <- if (length(cancers) == 1) c('treatment','Gender') else c('treatment', 'Cancer_Type', 'Gender')
  columns <- c(columns, label_column)
  metadata <- metadata[,columns]

  metadata <- na.omit(metadata)

  if (label_column == "Response") {
    labels <- factor(ifelse(metadata$Response == "Response", 1, 0))
  }else{
    labels <- factor(ifelse(metadata$RECIST == label, 1, 0))
  }

  columns <- columns[columns != label_column]
  metadata <- metadata[,columns]

  metadata_encoded <- fastDummies::dummy_cols(metadata, remove_first_dummy = TRUE, select_columns = columns, remove_selected_columns = TRUE)
  metadata_encoded$label <- labels
  rownames(metadata_encoded) <- rownames(metadata)

  return(metadata_encoded)
}

create_lightgbm_data <- function(scores = NULL,
                                metadata = NULL,
                                metadata.encoded = NULL,
                                genes = NULL,
                                use.survival = FALSE,
                                survival_time = NULL,
                                survival_status = NULL,
                                use.scores = c(),
                                use.genes = c(),
                                remove.cancertype = FALSE,
                                outcome = c(),
                                algo = c()) {


  data.out <-



  if (!is.null(metadata.encoded)) {
    if (use.scores == TRUE && use.genes == TRUE) {
      data = cbind(scores, metadata.encoded, genes)
    } else if (use.scores == TRUE && use.genes == FALSE) {
      data = cbind(scores, metadata.encoded)
    } else if (use.scores == FALSE && use.genes == TRUE) {
      data = cbind(genes, metadata.encoded)
    } else if (use.scores == FALSE && use.genes == FALSE) {
      data = as.data.frame(metadata.encoded)
    }
  } else {
    if (use.scores == TRUE && use.genes == TRUE) {
      data = cbind(scores, genes)
    } else if (use.scores == TRUE && use.genes == FALSE) {
      data = as.data.frame(scores)
    } else if (use.scores == FALSE && use.genes == TRUE) {
      data = as.data.frame(genes)
    }
  }

  if (use.survival == TRUE){

    data$survival_time = survival_time
    data$survival_status = survival_status
    data = data[!is.na(data$survival_time) & !is.na(data$survival_status),]

  } else {

    if(outcome=='Response') {
      data$outcome = as.numeric(metadata$Response=='Response')
    } else if (outcome == 'PD') {
      data$outcome = metadata$PD
    } else if (outcome == 'CR'){
      data$outcome = metadata$CR
    }

    data = data[!is.na(data$outcome),]
  }

  # if (!is.null(metadata.encoded)){
  #   if (remove.cancertype == FALSE){
  #     data = data %>% dplyr::select(-c('treatment.anti.CTLA4', 'Cancer_Type.NSCLC' ,'Gender.female'))
  #   } else {
  #     data = data %>% dplyr::select(-c('treatment.anti.CTLA4' ,'Gender.female'))
  #   }
  # }

  View(data)

  return(as.data.frame(data))

}


predict_response_lightgbm <- function(data,
                                    return_shap = FALSE,
                                     num_threads = 1){

  lightgbm::setLGBMthreads(1)


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


    # >>> Tune learning rate:
  dtrain = lightgbm::lgb.Dataset(data = train_x, label = train_y)

  initial_learning_rates <- c(0.001, 0.01, 0.1, 0.3, 0.5, 1, 2, 3, 5, 10)

  cv_results <-  parallel::mclapply(initial_learning_rates, function(lr) {
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

  # # Create a sequence of learning rates around the initial best value
  # refinement_factor <- 0.8   # Here we use a relative adjustment of ±80% of the initial best learning rate
  # learning_rates <- unique(round(seq(from = best_initial_lr * (1 - refinement_factor),
  #                       to = best_initial_lr * (1 + refinement_factor),
  #                       length.out = 50), 2))
  #
  # learning_rates <- learning_rates+0.001 # remove?
  #
  # cv_results <- parallel::mclapply(learning_rates, function(lr) {
  #   params <- list(objective = "binary",
  #                  metric = "auc",
  #                  learning_rate = lr)
  #
  #   lightgbm::lgb.cv(params = params,
  #                    data = dtrain,
  #                    nfold = 10,
  #                    nrounds = 1000,
  #                    early_stopping_rounds = 10,
  #                    verbose = -1)
  # }, mc.cores = num_threads)
  #
  # auc_values <- sapply(cv_results, function(res) {
  #   max(unlist(res$record_evals$valid$auc$eval))
  # })
  #
  # # Find the refined best learning rate
  # best_lr <- learning_rates[which.max(auc_values)]
  best_lr <- best_initial_lr

  # # >>> Tune other hyperparameters - stage 1:
  # dtrain = lightgbm::lgb.Dataset(data = train_x, label = train_y, params = list(feature_pre_filter = FALSE))
  #
  # tune_grid_stage1 <- expand.grid(
  #   num_leaves = c(6, 15, 31, 63, 100),
  #   max_depth = c(-1, 8, 15, 30),
  #   min_data_in_leaf = c(5, 10, 20, 40),
  #   feature_fraction = c(0.2, 0.6, 0.8, 1),
  #   bagging_fraction = c(0.2, 0.6, 0.8, 1),
  #   lambda_l1 = c(0, 0.5, 1),
  #   lambda_l2 = c(0, 2, 4, 6, 10)
  # )
  #
  #
  # tune_lgb_stage1 <- function(params) {
  #   lightgbm::lgb.cv(
  #     params = list(
  #       objective = "binary",
  #       metric = "auc",
  #       learning_rate = best_lr,
  #       num_leaves = params$num_leaves,
  #       max_depth = params$max_depth,
  #       min_data_in_leaf = params$min_data_in_leaf,
  #       feature_fraction = params$feature_fraction,
  #       bagging_fraction = params$bagging_fraction,
  #       lambda_l1 = params$lambda_l1,
  #       lambda_l2 = params$lambda_l2
  #     ),
  #     data = dtrain,
  #     nfold = 5,
  #     nrounds = 1000,
  #     early_stopping_rounds = 10,
  #     verbose = -1
  #   )
  # }
  #
  # results_stage1 <- parallel::mclapply(1:nrow(tune_grid_stage1), function(i) {
  #   params <- tune_grid_stage1[i, ]
  #   cv_result <- tune_lgb_stage1(params)
  #   auc <- max(unlist(cv_result$record_evals$valid$auc$eval))
  #   list(params = params, auc = auc)
  # }, mc.cores = num_threads)
  #
  # # Find the best hyperparameters from Stage 1
  # best_result_stage1 <- results_stage1[[which.max(sapply(results_stage1, function(res) res$auc))]]
  # best_params_stage1 <- best_result_stage1$params
  #
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
  #       max_depth = best_params_stage1$max_depth,
  #       min_data_in_leaf = best_params_stage1$min_data_in_leaf,
  #       feature_fraction = best_params_stage1$feature_fraction,
  #       bagging_fraction = best_params_stage1$bagging_fraction,
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
  #
  #
  # Fit best model
  #best_params <- cbind(best_params_stage1, best_params_stage2, "learning_rate" = best_lr)

  train_params <- list(
    boosting = "dart",            # Use DART boosting
    drop_rate = 0.01,
    skip_drop = 0.03, # DART-specific parameter
    feature_fraction = 0.5,
    feature_pre_filter = FALSE,
    max_bin = 300
  )
  dtrain = lightgbm::lgb.Dataset(data = train_x, label = train_y, params=train_params)



  lightgbm::setLGBMthreads(num_threads)
  # dtrain <- lightgbm::lgb.Dataset(data = train_x, label = train_y, params = list(max_bin = best_params_stage2$max_bin))
  dtest <- lightgbm::lgb.Dataset(data = test_x, label = test_y)


  # best_model <- lightgbm::lgb.train(
  #   data = dtrain,
  #   params = list(
  #     objective = "binary",
  #     metric = "auc",
  #     learning_rate = best_params$learning_rate,
  #     num_leaves = best_params$num_leaves,
  #     min_data_in_leaf = best_params$min_data_in_leaf,
  #     feature_fraction = best_params$feature_fraction,
  #     bagging_fraction = best_params$bagging_fraction,
  #     lambda_l1 = best_params$lambda_l1,
  #     lambda_l2 = best_params$lambda_l2,
  #     min_split_gain = best_params$min_split_gain,
  #     min_sum_hessian_in_leaf = best_params$min_sum_hessian_in_leaf,
  #     extra_trees = best_params$extra_trees,
  #     path_smooth = best_params$path_smooth,
  #     boost_from_average = best_params$boost_from_average
  #   ),
  #   valids = list(val = dtest),
  #   nrounds = 10000,
  #   early_stopping_rounds = 100,
  #   verbose = -1)

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





  # --------> Almog review Loai's code so far

  # SHAP
  # new_labels <- list()
  # shap_long_for_train <- SHAPforxgboost::shap.prep(xgb_model = best_model, X_train = train_x)
  # shap_long_for_test <- SHAPforxgboost::shap.prep(xgb_model = best_model, X_train = test_x)

  # shap_imp <- SHAPforxgboost::shap.importance(shap_values_long, top_n = 20)
  # shapplot_for_train <- SHAPforxgboost::shap.plot.summary.wrap1(best_model, train_x, top_n = 20)
  # shapplot_for_test <- SHAPforxgboost::shap.plot.summary.wrap1(best_model, test_x, top_n = 20)


  # train_pred = predict(best_model, train_x)
  # train_roc = pROC::roc(train_y, train_pred)
  # train_auc = pROC::auc(train_roc)
  # cat("Train AUC:", train_auc, "\n")


  if (return_shap) {
    shap_values <- shapviz::shapviz(best_model, X_pred = as.matrix(test_x))

    return(shap_values)
    # shapviz::sv_importance(shap_values, kind = "beeswarm")



  }



  test_pred = predict(best_model, test_x)
  test_roc = pROC::roc(test_y, test_pred)
  test_auc = pROC::auc(test_roc)
  # cat("Test AUC:", test_auc, "\n")
  #
  # train_df = data.frame(row.names = rownames(train_x), pred = train_pred)
  # test_df = data.frame(row.names = rownames(test_x), pred = test_pred)
  #
  # print(' - - - - - - - - - - - - - - - - - - - - - RUN ENDED - - - - - - - - - - - - - - - - - - - - - ')
  # return(list(test_df = test_df, train_df = train_df, best_model = best_model, best_params = best_params,
  #             test_auc = test_auc, train_auc = train_auc, algorithm = algo,
  #             shap_long_for_train = shap_long_for_train, shap_long_for_test = shap_long_for_test))

  auc <- as.numeric(gsub("Area under the curve: ", "", test_auc))
  return(auc)

}
