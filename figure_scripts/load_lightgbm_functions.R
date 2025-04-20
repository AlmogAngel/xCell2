library(tidyverse)
library(lightgbm)
library(caret)
library(pROC)
library(shapviz)


# Define hyperparameter bounds
bounds <- list(
  learning_rate       = c(0.001, 0.10),
  num_leaves          = c(4L, 64L),
  max_depth           = c(3L, 10L),
  min_data_in_leaf    = c(5L, 50L),
  feature_fraction    = c(0.5, 1.0),
  bagging_fraction    = c(0.5, 1.0),
  bagging_freq        = c(0L, 5L),
  lambda_l1           = c(0.0, 5.0),
  lambda_l2           = c(0.0, 5.0),
  drop_rate           = c(0.0, 0.5),
  skip_drop           = c(0.0, 1.0)
)

# Function to sample a random value from bounds
sample_value <- function(lb, ub) {
  if (is.integer(lb)) {
    # For integer parameters, sample from the sequence of integers
    return(as.integer(sample(seq(lb, ub), 1)))
  } else {
    # For continuous parameters, sample uniformly
    return(runif(1, lb, ub))
  }
}

# Function to generate a random set of hyperparameters from bounds
sample_params <- function(bounds) {
  params <- list()
  for (name in names(bounds)) {
    lb <- bounds[[name]][1]
    ub <- bounds[[name]][2]
    params[[name]] <- sample_value(lb, ub)
  }
  return(params)
}

# Define a function that runs one random search iteration
run_rs_iteration <- function(iter, train_matrix, train_label, categorical_features, bounds) {
  
  lightgbm::setLGBMthreads(1)
  
  
  dtrain <- lgb.Dataset(data = train_matrix,
                        label = train_label,
                        categorical_feature = categorical_features,
                        params = list(feature_pre_filter = FALSE))
  
  
  # Generate random hyperparameters
  hp <- sample_params(bounds)
  
  # Set up LightGBM parameters using the sampled hyperparameters
  params <- list(
    objective = "binary",
    metric    = "auc",
    boosting  = "dart",
    learning_rate    = hp$learning_rate,
    num_leaves       = hp$num_leaves,
    max_depth        = hp$max_depth,
    min_data_in_leaf = hp$min_data_in_leaf,
    feature_fraction = hp$feature_fraction,
    bagging_fraction = hp$bagging_fraction,
    bagging_freq     = hp$bagging_freq,
    lambda_l1        = hp$lambda_l1,
    lambda_l2        = hp$lambda_l2,
    drop_rate        = hp$drop_rate,
    skip_drop        = hp$skip_drop,
    force_row_wise   = TRUE  # avoid sparse format issues
  )
  
  # Perform cross-validation with LightGBM
  cv <- lgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    stratified = TRUE,
    verbose = -1
  )
  
  # Compute the best AUC from CV results
  best_auc <- max(unlist(cv$record_evals$valid$auc$eval))
  
  # Return the hyperparameters and the corresponding AUC
  return(tibble(
    learning_rate = hp$learning_rate,
    num_leaves = hp$num_leaves,
    max_depth = hp$max_depth,
    min_data_in_leaf = hp$min_data_in_leaf,
    feature_fraction = hp$feature_fraction,
    bagging_fraction = hp$bagging_fraction,
    bagging_freq = hp$bagging_freq,
    lambda_l1 = hp$lambda_l1,
    lambda_l2 = hp$lambda_l2,
    drop_rate = hp$drop_rate,
    skip_drop = hp$skip_drop,
    best_auc = best_auc
  ))
}


run_lightgbm <- function(data, outer_iterations, ncores){
  
  # For reproducibility
  set.seed(123)
  lightgbm::setLGBMthreads(1)
  
  # Set up outer loop parameters
  auc_values <- numeric(outer_iterations)
  model_predictions <- vector("list", outer_iterations)
  shap_values_list <- vector("list", outer_iterations)

  
  train_idx_list <- createDataPartition(y = data$Response, p = 0.75, list = TRUE, times = outer_iterations) # outer_iterations different random splits
  
  # Outer loop
  for(i in 1:length(train_idx_list)){
    
    # --- Data splitting (25% test, 75% train) ---
    train_data <- data[train_idx_list[[i]], ]
    test_data  <- data[-train_idx_list[[i]], ]
    
    # Create matrices for LightGBM (exclude Response)
    train_matrix <- data.matrix(train_data[, setdiff(names(train_data), "Response")])
    train_label <- train_data$Response
    test_matrix  <- data.matrix(test_data[, setdiff(names(test_data), "Response")])
    test_label <- test_data$Response
    
    # Create LightGBM datasets (pass categorical feature names)
    dtrain <- lgb.Dataset(data = train_matrix,
                          label = train_label,
                          categorical_feature = categorical_features,
                          params = list(feature_pre_filter = FALSE))
    
    
    # --- Tune learning rate using Bayesian Optimization ---
    
    
    # Run the random search iterations in parallel using mclapply
    results_list <- parallel::mclapply(1:(ncores*5), function(x){
      run_rs_iteration(iter = x, train_matrix, train_label, categorical_features, bounds)
    }, mc.cores = ncores)
    
    hp_tuning_results <- bind_rows(results_list) |> 
      arrange(-best_auc)
    
    best_params <- hp_tuning_results[1,] |> 
      select(-best_auc) |> 
      as.list()
    
    # --- Train final model ---
    
    final_params <- list(
      objective = "binary",
      metric    = "auc",
      boosting  = "dart"
    )
    final_params <- c(final_params, best_params)
    
    final_model <- suppressWarnings(lgb.train(
      params = final_params,
      data   = dtrain,
      nrounds = 1000,
      verbose = -1
    ))
    
    # --- Evaluate model on test set ---
    preds <- suppressMessages(predict(final_model, test_matrix))
    auc_value <- as.numeric(auc(test_data$Response, preds))
    auc_values[i] <- auc_value
    
    # Store the predictions and AUC for this iteration
    model_predictions[[i]] <- list(predictions = preds, auc = auc_value)
    
    # --- Compute and store SHAP values using shapviz ---
    # Compute SHAP values for the test set predictions
    shap_obj <- shapviz(final_model, X_pred = test_matrix)
    shap_values_list[[i]] <- shap_obj
    
    cat(sprintf("Iteration %d / %d, Test AUC: %.4f\n", i, outer_iterations, auc_value))
  }
  
  # --- Summarize overall performance ---
  median_auc_overall <- median(auc_values)
  cat("Median AUC across 100 models:", median_auc_overall, "\n")
  
  
  out <- list(auc_values = auc_values,
              model_predictions = model_predictions,
              shap_values_list = shap_values_list
              )

  return(out)
  
}





