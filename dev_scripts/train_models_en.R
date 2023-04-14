# Elastic Net

# ref<-ref_old
# labels<-labels_old
# ref <-ref[,test]
# labels<-labels[test,]


trainModels <- function(ref, labels, pure_ct_mat_test_test, signatures_collection_filtered, dep_list){
  
  celltypes <- unique(labels[,2])

  sig_models_list <- lapply(celltypes[celltypes != "T-cells"], function(ctoi){
    

    
    signature_filtered_ctoi <- signatures_collection_filtered[startsWith(names(signatures_collection_filtered), paste0(ctoi, "#"))]
    
    print(ctoi)
    # # Make mixture
    # mix <- makeMixture(ctoi, ref, labels, pure_ct_mat_test, dep_list)
    # mix_ranked <- singscore::rankGenes(mix$mixture)
    # scores_mat <- sapply(signature_filtered_ctoi, simplify = TRUE, function(sig){
    #   singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    # })
    
    # # Train Elastic Net models with different alphas
    # alphas <- seq(0, 1, 0.2)
    # models <- lapply(alphas, function(x){
    #   glmnet::cv.glmnet(scores_mat, mix$ctoi_fracs, alpha = x, family = "gaussian")
    # })
    # 
    # sig_models <- tibble(celltype = ctoi,
    #                      signatures = list(signature_filtered_ctoi),
    #                      n_sigs = length(signature_filtered_ctoi),
    #                      model = models,
    #                      alpha = alphas)
    # 
    # en_models <- sig_models %>%
    #   rowwise() %>%
    #   mutate(predictions = list(as.numeric(predict(model, newx = scores_mat)))) %>%
    #   mutate(spearman = cor(mix$ctoi_fracs, predictions, method = "spearman"),
    #          pearson = cor(mix$ctoi_fracs, predictions, method = "pearson")) %>%
    #   mutate(score = sum(spearman, pearson)) %>%
    #   arrange(desc(score))
    # 
    # en_models
    
    # Train RF models
    # randomForest::tuneRF(x=scores_mat, y=mix$ctoi_fracs, ntreeTry = 100, mtryStart = 22, stepFactor = 1.5)

    # ntree_values <- c(100, 1000)
    # models <- lapply(ntree_values, function(x){
    #   randomForest::randomForest(scores_mat, mix$ctoi_fracs, ntree = x, mtry = 22)
    # })
    
    mix_train_ranked <-  singscore::rankGenes(mix_train_list[[ctoi]]$mixture)
    scores_mat_train <- sapply(signature_filtered_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_train_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    
    mix_test_ranked <-  singscore::rankGenes(mix_test_list[[ctoi]]$mixture)
    scores_mat_test <- sapply(signature_filtered_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_test_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
    
    models <- list(randomForest::randomForest(scores_mat_train, mix_train_list[[ctoi]]$ctoi_fracs))

    sig_models <- tibble(celltype = ctoi,
                         signatures = list(signature_filtered_ctoi),
                         n_sigs = length(signature_filtered_ctoi),
                         model = models)

    rf_models <- sig_models %>%
      rowwise() %>%
      mutate(predictions = list(as.numeric(predict(model, newdata = scores_mat_test)))) %>%
      mutate(spearman = cor(mix_test_list[[ctoi]]$ctoi_fracs, predictions, method = "spearman"),
             pearson = cor(mix_test_list[[ctoi]]$ctoi_fracs, predictions, method = "pearson"))
    
    rf_models
    
    # rbind(en_models[1,], rf_models[1,])

  })
  
  
  saveRDS(sig_models_list, "/bigdata/almogangel/xCell2_data/sig_models_list_en.rds")
  
  names(sig_models_list) <- celltypes[celltypes != "T-cells"]
  models <- sig_models_list %>% 
    bind_rows() %>% 
    group_by(celltype) %>% 
    arrange(desc(score)) %>% 
    top_n(1)
  
  return(models)
  
}
