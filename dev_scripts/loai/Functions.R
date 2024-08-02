library(xCell)


################## Get Ensemble ID for each gene

GetEnsmblID = function(dataset) {    ## Get ensemble ID

  geneSymbols = as.vector(rownames(dataset))
  geneIDs = ensembldb::select(EnsDb.Hsapiens.v86, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
  geneIDs = geneIDs[!duplicated(geneIDs$SYMBOL),]
  rownames(geneIDs) = geneIDs$SYMBOL
  ind = intersect(rownames(geneIDs), rownames(dataset))
  dataset = dataset[ind,]
  rownames(dataset) = geneIDs$GENEID
  return(dataset)
}

################## Get the transcript length of each gene

GetGeneLength2 = function(counts){

  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genelength =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id', 'transcript_length','cds_length'), filters =  'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = FALSE)
  gene_canonical_transcript =  getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_is_canonical'), filters =  'ensembl_gene_id', values = rownames(counts), mart = ensembl, useCache = FALSE)
  gene_canonical_transcript_subset = gene_canonical_transcript[!is.na(gene_canonical_transcript$transcript_is_canonical),]
  genelength = merge(gene_canonical_transcript_subset, genelength, by = c("ensembl_gene_id", "ensembl_transcript_id"))
  return(genelength)

}

################## Get the length of each gene (different approach - less recommended)

# GetGeneLength <- function(dataset) {    ## Gene length
#
#   exons = exonsBy(EnsDb.Hsapiens.v86, by="gene")
#   exons = reduce(exons)
#   len = sum(width(exons))
#   insect = intersect(rownames(dataset),names(len))
#   genelength = len[insect]
#   #dataset = dataset[insect,]
#   return(genelength)
#
# }

################## Remove unwated characters from rownames or colnames

clean_names = function(names) {
  names %>%
    gsub("^X", "", .) %>%
    gsub("\\.|\\-|\\s|\\+", "_", .) %>%
    gsub(",", "_", .)
}

################## FUN

# Normalizedata <- function(dataset, genelength){  ## Normalize TPM
#
#   rpkm <- apply(X = subset(dataset),
#                 MARGIN = 2,
#                 FUN = function(x) {
#                   10^9 * x / genelength / sum(as.numeric(x))
#                 })
#
#   TPM <- apply(rpkm, 2, function(x) x / sum(as.numeric(x)) * 10^6) %>% as.data.frame()
#   return(TPM)
# }

################## Normalize data, from counts to TPM - code by Michael Love

countToTpm = function(counts = NULL, genelength = NULL) {

  x = counts / genelength
  tpm.mat = t( t(x) * 1e6 / colSums(x) )
  return(tpm.mat)

}

################## Remove genes according to threshold

CleanData = function(TPM = NULL, t = NULL, k = NULL) {

  TPM = TPM[rowSums(TPM[])>0,]
  thresh = TPM > t
  keep = rowSums(thresh) >= k
  table(keep)
  TPM = TPM[keep,]
  return(TPM)

}

################## Get the symbol of each gene

GetHugoSymbols = function(TPM = NULL, use.column = FALSE) {

  ensembl.genes = as.vector(rownames(TPM))
  geneIDs = ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  geneIDs = geneIDs[!duplicated(geneIDs$GENEID),]
  rownames(geneIDs) = geneIDs$GENEID

  if(nrow(TPM) != nrow(geneIDs)) {
    stop("Mismatch between TPM rows and geneIDs fetched.")
  }

  if (use.column == FALSE){
    rownames(TPM) = geneIDs$SYMBOL
  } else {
    TPM$row_name = geneIDs$SYMBOL
  }

  return(TPM)
}

################## Perform upper quantile normalization

upperQuantile = function(myTPM = NULL) {
  tpm.quartiles.expressed = apply(myTPM, 2, function(x){quantile(x[x>0], 0.75)})
  uqn.tpm.expressed = t(t(myTPM)/ tpm.quartiles.expressed)
  return(uqn.tpm.expressed)
}

################## Load the separate data sets (as a combined data frame or separated)

load_tpm_data = function(converted.to.TPM = NULL,
                         originally.TPM = NULL,
                         use.quantile.Norm = FALSE,
                         use.log2 = FALSE,
                         use.mini = FALSE,
                         new.data = TRUE,
                         return.dataframe = TRUE){


  if (new.data == TRUE){
    if(use.mini == FALSE){
      TPMs = dplyr::lst(originally.TPM[[1]],
                        originally.TPM[[2]],
                        converted.to.TPM[[1]],
                        originally.TPM[[3]],
                        originally.TPM[[4]],
                        originally.TPM[[5]],
                        originally.TPM[[6]],
                        originally.TPM[[7]],
                        originally.TPM[[8]],
                        converted.to.TPM[[2]],
                        converted.to.TPM[[3]],
                        originally.TPM[[9]],
                        originally.TPM[[10]],
                        originally.TPM[[11]],
                        originally.TPM[[12]],
                        originally.TPM[[13]],
                        originally.TPM[[14]],
                        originally.TPM[[15]],
                        converted.to.TPM[[4]],
                        originally.TPM[[16]],
                        originally.TPM[[17]])
    } else{
      TPMs = dplyr::lst(originally.TPM[[1]],
                        # originally.TPM[[2]],
                        converted.to.TPM[[1]],
                        originally.TPM[[3]],
                        originally.TPM[[4]],
                        originally.TPM[[5]],
                        originally.TPM[[6]], #
                        originally.TPM[[7]],
                        originally.TPM[[8]],
                        converted.to.TPM[[2]],
                        converted.to.TPM[[3]],
                        originally.TPM[[9]],
                        originally.TPM[[10]],
                        originally.TPM[[11]], # 16
                        originally.TPM[[12]],
                        originally.TPM[[13]],
                        originally.TPM[[14]],
                        originally.TPM[[15]],
                        converted.to.TPM[[4]],
                        originally.TPM[[16]],
                        originally.TPM[[17]])
    }

  } else {

    TPMs = dplyr::lst(converted.to.TPM[[1]],
                      originally.TPM[[1]],
                      converted.to.TPM[[2]],
                      originally.TPM[[2]],
                      converted.to.TPM[[3]],
                      converted.to.TPM[[4]],
                      originally.TPM[[3]],
                      originally.TPM[[4]],
                      originally.TPM[[5]],
                      converted.to.TPM[[5]],
                      converted.to.TPM[[6]],
                      originally.TPM[[6]],
                      originally.TPM[[7]] ,
                      originally.TPM[[8]],
                      converted.to.TPM[[7]],
                      converted.to.TPM[[8]],
                      converted.to.TPM[[9]],
                      originally.TPM[[9]],
                      converted.to.TPM[[10]],
                      originally.TPM[[10]],
                      originally.TPM[[11]])


  }

  if (use.quantile.Norm == TRUE){
    for (i in seq_along(TPMs)){
      TPMs[[i]] = upperQuantile(TPMs[[i]])
    }
  }

  if (use.log2 == TRUE){
    for (i in seq_along(TPMs)){
      TPMs[[i]] = log2(TPMs[[i]] + 1)
    }
  }

  TPMs = lapply(TPMs, function(x) {x$rowname <- rownames(x); return(x)})

  if (return.dataframe == TRUE){
    TPMs_combined = Reduce(function(x, y) merge(x, y, by="rowname", all=FALSE), TPMs)
    rownames(TPMs_combined) = TPMs_combined$rowname
    TPMs_combined$rowname = NULL
    return(TPMs_combined)
  } else {
    TPMs = lapply(TPMs, function(df) df[, -ncol(df), drop = FALSE])
    return(TPMs)
  }

}

################## FUN

MatchData = function(TPM = NULL, meta_data = NULL) {   ## Order / Match the data

  insect = intersect(colnames(TPM), rownames(meta_data))
  TPM = TPM[, insect]
  meta_data = meta_data[insect,]
  return(TPM)

}

################## FUN

MatchMeta = function(TPM, meta_data) {   ## Order / Match the Metadata

  insect = intersect(colnames(TPM), rownames(meta_data))
  TPM = TPM[, insect]
  meta_data = meta_data[insect,]
  return(meta_data)

}


################## Find the model with the closest AUC to the median of all AUC scroes

aucClosest = function(models_list = NULL, metric = NULL) {

  auc_scores = sapply(models_list, function(model) model[[14]])
  if (metric == "mean"){
    median_auc = mean(auc_scores)
  } else {
    median_auc = median(auc_scores)
  }

  differences = abs(auc_scores - median_auc)
  closest_index = which.min(differences)

  return(closest_index)
}

################## Run xCell1

RunxCell = function(transcriptomic_Data = NULL) {

  xCellscores = xCell::xCellAnalysis(transcriptomic_Data,
                              signatures = NULL,
                              genes = NULL,
                              spill  = NULL,
                              rnaseq = TRUE,
                              file.name = NULL,
                              scale = TRUE,
                              alpha = 0.5,
                              save.raw = FALSE,
                              parallel.sz = 8,
                              parallel.type = "SOCK",
                              cell.types.use = cells_2use)
  return (xCellscores)
}


################## Choose the target binary column and perform wilcoxon test

wilcoxonTest = function(scoresData = NULL, metadata = NULL, column = NULL, group1_indices = NULL) {

  # Determine indices based on column value
  if (!is.null(column) && length(column) > 0) {
    if (column %in% c('Response','Benefit')) {
      idx1 = which(metadata[[column]] == 'Response')
      idx2 = which(metadata[[column]] == 'NoResponse')
    } else if (column %in% c('PD','CR')) {
      idx1 = which(metadata[[column]] == 0)
      idx2 = which(metadata[[column]] == 1)
    }
  } else if (is.null(column) || length(column) == 0) {
    idx1 = group1_indices
    idx2 = setdiff(1:nrow(scoresData), idx1)
  }

  # Use apply() to compute the p-values
  p = apply(scoresData, 1, function(row) {
    wilcox.test(row[idx1], row[idx2], exact = FALSE)$p.value
  })

  lfc = log2(rowMeans(scoresData[, idx1]) / rowMeans(scoresData[, idx2]))
  p.adj = p.adjust(p, method = 'fdr')

  return(list(p, p.adj, lfc))
}


################## Visualize the scores for a cell type using boxplots

cellBox = function(scoresData = NULL ,metadata = NULL, column = NULL, celltype = NULL, k = NULL) {

  resp = data.frame(sample = colnames(scoresData),
                    response = metadata[[deparse(substitute(column))]])

  scoresData %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    right_join(resp, .) %>%
    ggplot(aes(x = response, y = {{celltype}})) +
    geom_boxplot(aes(fill = response)) +
    scale_fill_manual(values = c('red','lightblue')) +
    geom_jitter(color="black", size=0.8, alpha=0.9) +
    ylim(0, k)
}


################## Visualize the data using t-SNE

tsnePlot = function(data = NULL, column = NULL, metadata = NULL, p = NULL, seed = NULL, color_vec = NULL) {    ## t-SNE for scores

  if(!is.null(seed)){
    set.seed(seed)
  }

  tsne_result = Rtsne(data, perplexity = p)

  tsne_data = data.frame(
    x = tsne_result$Y[, 1],
    y = tsne_result$Y[, 2],
    group_info = as.factor(metadata[[column]]))

  tsne_data = na.omit(tsne_data)

  # Plot the t-SNE results using ggplot2
  tsne_plot = ggplot(tsne_data, aes(x = x, y = y, color = group_info)) +
    geom_point(size = 3) +
    theme_classic() +
    labs(color = column) +
    theme_bw() +
    ggtitle(paste('t-SNE of', column)) +
    theme(
      title = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12)
    )


  if(!is.null(color_vector)){
    tsne_plot = tsne_plot + scale_color_manual(values = color_vector)
  }

  return(tsne_plot)

}

################## Visualize the scores using Heatmap

scoresHeatmap = function(scoresData = NULL, metadata = NULL, column1 = NULL, column2 = NULL){     ## Heatmap cell

  scores_data_scaled = t(scale(t(scoresData)))
  scores_data_scaled[scores_data_scaled > 1] <- 2
  scores_data_scaled[scores_data_scaled < -1] <- -2
  return(pheatmap(scores_data_scaled ,annotation_col = metadata[,c(column1,column2)],
                  cluster_cols = TRUE, show_colnames = F, clustering_method = 'ward.D'))

}

################## Find the most variant genes

calcVariation = function(data){

  Var = apply(data, 1, var)
  selectedGenes1000 = names(Var[order(Var, decreasing = T)][1:1000])
  selectedGenes50 = names(Var[order(Var, decreasing = T)][1:50])
  listofvar = list(Var, selectedGenes1000, selectedGenes50)
  return(listofvar)

}

################## Visualize the transcriptomic data using Heatmap

geneticHeatmap = function(data = NULL ,metadata = NULL, column1 = NULL, column2 = NULL, VAR, n , show.rownames = TRUE){   ## Choose VAR[[2]] or Var[[3]]

  data_scaled = t(scale(t(data)))
  data_scaled[data_scaled > 1] = 2
  data_scaled[data_scaled < -1] = -2
  pheatmap(data_scaled[VAR[[n]],] ,annotation_col = meta_data[,c(column1,column2)],
           cluster_cols = TRUE, show_colnames = FALSE,cluster_rows = TRUE, show_rownames = show.rownames, clustering_method = 'ward.D')

}

################## Differential expression analysis using LIMMA

limmaNorm = function(metadata = NULL, dataset = NULL, column = NULL){

  d = DGEList(dataset)
  group = as.factor(metadata[[column]])
  d$samples$group = group
  keep.exprs = edgeR::filterByExpr(d, group = group)
  d = d[keep.exprs,, keep.lib.sizes = FALSE]
  d = calcNormFactors(d , method = 'TMM' )
  return(d)

}

################## Differential expression analysis using LIMMA

limmaAnalysis = function(d = NULL, design = NULL , contrast = NULL) {

  Voom = voomWithQualityWeights(d, design, plot = TRUE)
  vfit = lmFit(Voom, design = design) # Fit linear model for each gene
  vfit = contrasts.fit(vfit , contrasts = contrast) # compute statistics and fold changes for comparison interests.
  efit = eBayes(vfit) # computes more precise estimates by sharing gene information using empirical bayes moderation
  return(efit)

}

################## Differential expression analysis using LIMMA (for TPM data, not counts!)

limmaAnalysisTPM = function(data = NULL ,num = NULL, design = NULL, contrast = NULL) {

  lTPM = log2(data + num)
  print(dim(lTPM))
  print(dim(design))
  print(all(colnames(lTPM) == rownames(design)))
  vfit = lmFit(lTPM, design)
  vfit = contrasts.fit(vfit, contrasts = contrast)
  efit = eBayes(vfit, robust = TRUE, trend = TRUE)
  return(efit)

}

################################# ML #########################################

sco_meta = function(metadat,scores){

  INDx = intersect(metadat$samples,colnames(scores))
  Scores_new = scores[,INDx]
  rownames(metadat) = metadat$samples
  metadat = metadat[INDx,]
  ls = list(metadat, Scores_new)
  return(ls)
}

################## Add three model features for xCell1 results

MES = function(xcell.scores = NULL) { # Get the three additional features for xcell1
  ImmuneScore = apply(xcell.scores[c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells'),],2,sum)/1.5
  StromaScore = apply(xcell.scores[c('Adipocytes','Endothelial cells','Fibroblasts'),],2,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore
  xcell.scores = rbind(ImmuneScore,StromaScore,MicroenvironmentScore)
  return(xcell.scores)
}

################## Make dummy variables for the data

make_dumm = function(data){

  data$Gender = as.numeric(factor(data$Gender))-1

  data_dummy = dummy_cols(data)
  rownames(data_dummy) = rownames(data)
  data_dummy = data_dummy %>% select(-c(Cancer_Type,treatment))

  if (colnames(data)[1] == 'CR') {
    data_dummy[,1] = plyr::mapvalues(data_dummy[,1], c(1,0), c('Complete','nonComplete'))
  } else if (colnames(data)[1] == 'PD'){
    data_dummy[,1] = plyr::mapvalues(data_dummy[,1], c(1,0), c('Prog','nonProg'))
  }

  return(data_dummy)

}

################## Show distributions

showDistribution = function(data = NULL, column1 = NULL, column2 = NULL) {

  ggplot(data,
         aes(x = data[[column2]],
             group = data[[column1]],
             fill = data[[column1]])) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()

}

################## Waterfall plot - useful for evaluating ML classifiers (Waterfall2 is better)

waterfall = function(model,metric,meta,test) {

  d = test
  ind = intersect(rownames(d), rownames(meta))
  test_meta = meta[ind,]

  sel = rownames(test_meta)
  Recist_test = meta[sel,]$RECIST

  pred = predict(model, data.matrix(test[,-1]))
  df = data.frame(row.names = rownames(d),xcell = scale(pred), truth=d[[metric]], type=test_meta$Cancer_Type)
  df =  df[order(df$xcell),]

  ggplot(df[order(df$xcell),], aes(seq(nrow(df)), scale(xcell), fill = as.factor(Recist_test))) +
    geom_col() +
    labs(x = 'Test set', y = 'scaled xCell2 prediction scores')+ #ggtitle(paste('Test of',title,'VS all ')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())}


waterfall2 = function(metadata = NULL,
                      test.pred = NULL,
                      case = c(),
                      v.line = NULL,
                      y_label = 'xCell prediction scores',
                      train.recist = c(),
                      test.recist = c())
{

  test.pred$shit = 'shit'
  ind = intersect(rownames(test.pred), rownames(metadata))
  test_meta = metadata[ind,]
  test.pred = test.pred[ind,]
  pred = test.pred[,1]
  mytruth = as.factor(test_meta[[case]])
  cat('My truth is the number of real', case, '0 :', sum(mytruth == 0), '1 :', sum(mytruth == 1))
  print(table(test_meta$Response))

  df = data.frame(row.names = rownames(test.pred),
                  Algo_scaled = as.numeric(scale(pred)),
                  Algo = as.numeric(pred),
                  truth = mytruth,
                  type = test_meta$Cancer_Type,
                  RECIST = test_meta$RECIST)

  df = df[order(df$Algo_scaled),]

  total_above_zero = sum(df$Algo_scaled > 0)
  red_above_zero = sum(df$truth == 1 & df$Algo_scaled > 0)
  total_below_zero = sum(df$Algo_scaled <= 0)
  red_below_zero = sum(df$truth == 1 & df$Algo_scaled <= 0)

  count_df = data.frame(
    total = c(total_above_zero, total_below_zero),
    red = c(red_above_zero, red_below_zero),
    direction = c("Above", "Below")
  )

  # Adjusted the fill aesthetic to RECIST, and the color scale to accommodate the 4 RECIST values.
  plt = ggplot(df, aes(seq(nrow(df)), Algo_scaled, fill = RECIST)) +
    geom_col(show.legend = TRUE) +
    labs(x = 'Test set', y = y_label) +
    scale_fill_manual(values = c("CR" = "dodgerblue3", "PR" = "yellow2", "SD" = "springgreen3", "PD" = "firebrick1")) +
    annotate("text", x = 120, y = max(df$Algo_scaled), size = 6,
             label = paste("Above zero: ", count_df$total[1], " total bars, ", count_df$red[1], " red bars")) +
    annotate("text", x = 120, y = min(df$Algo_scaled), size = 6,
             label = paste("Below zero: ", count_df$total[2], " total bars, ", count_df$red[2], " red bars")) +
    geom_vline(aes(xintercept = v.line), color = "orange", size = 1) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 17),         # Axis titles
      axis.text = element_text(size = 15),          # Axis ticks
      legend.text = element_text(size = 15),        # Legend text
      legend.title = element_text(size = 17),       # Legend title
      plot.title = element_text(size = 20, hjust = 0.5)   # Main title
    )  + ggtitle('Waterfall Plot for xCell 2.0 Predictions')

  plot(plt)

  return(list(plt,df))
}

################## Plot feature importance of a model (regardless of the direction)

plot_top_feature_importance = function(feature_importance = NULL, n = 15){

  # Sort and get the top n features
  sorted_importance = sort(feature_importance, decreasing = TRUE)
  top_n_features = head(sorted_importance, n)

  # Create a dataframe for ggplot
  importance_df = data.frame(
    feature = names(top_n_features),
    importance = top_n_features
  )

  # Plotting using ggplot
  ggplot(importance_df, aes(x = reorder(feature, importance), y = importance)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    coord_flip() + # for horizontal bars
    labs(title = paste("Top", n, "Features Importance"),
         x = "",
         y = "Importance") +
    theme_minimal()
}

################## Insert a model and see the auc of the model

auc_xgb = function(model){

  # xgb_trn = roc(xCell2_MLList[[1]][,1],predict(model1, dtrain))
  # plot(xgb_trn, main = 'xgb train')
  # auc(xgb_trn)

  xgb_tst = roc(xCell2_MLList[[2]][,1], predict(model1, data.matrix(xCell2_MLList[[2]][,-1])))
  #plot(xgb_tst, main = 'xgb test')
  return(auc(xgb_tst))

}

################## Accept data, split it, and return the train test sets

create_data_for_learning = function(scores = NULL, metadata = NULL, target='RECIST',add.scores=TRUE,
                                    metadata.use = c('Gender','treatment','Cancer_Type'),
                                    recist.use=c(), case='CR', train.size=0.75, seed=0) {

  if (seed>0) {
    set.seed(seed)
  }

  metadata = metadata %>% dplyr::filter(RECIST %in% recist.use)
  scores = scores[rownames(metadata),]

  if (add.scores == FALSE) {
    scores = cbind(metadata[,target],metadata[,metadata.use])
  } else {
    scores = cbind(metadata[,target],metadata[,metadata.use],scores)
  }

  colnames(scores) = make.names(colnames(scores))
  scores[,1] = as.numeric(scores[,1]==case)

  scores = dummy_cols(scores)
  scores = scores[,!(colnames(scores) %in% metadata.use)]
  colnames(scores) = make.names(colnames(scores))

  scores[,1] = factor(scores[,1])
  IND = createDataPartition(y = scores[,1], p=train.size, list = FALSE)
  data_train = scores[IND, ]
  data_test = scores[-IND,]

  list(data_train,data_test,seed)

}

################## FUN

train_ici = function(data = NULL, cpus=32) {

  data_train = data[[1]]
  data_test = data[[2]]

  # Create an xgboost classification learner with probability prediction type
  lrn <- makeLearner("classif.xgboost", predict.type = "prob")
  # Set the learner's hyperparameters
  lrn$par.vals <- list(objective="binary:logistic", eval_metric="auc",
                       nrounds=100L, eta=0.1)

  params <- makeParamSet(makeDiscreteParam("booster", values = c("gbtree","gblinear")),
                         makeNumericParam("eta", lower = 0.1, upper = 0.11),
                         makeIntegerParam("max_depth", lower = 3L, upper = 10L),
                         makeNumericParam("min_child_weight", lower = 1L, upper = 10L),
                         makeNumericParam("subsample", lower = 0.5, upper = 1),
                         makeNumericParam("colsample_bytree", lower = 0.5, upper = 1))

  # Create a resampling description for cross-validation with stratified sampling
  # and 5 iterations
  rdesc <- makeResampleDesc("CV", stratify = T, iters = 5L)
  # Create a random search control for the parameter tuning process
  ctrl <- makeTuneControlRandom(maxit = 10L)

  # Start parallel processing with the specified number of CPUs
  parallelStartSocket(cpus = cpus)

  # Rename the first column of both the training and test sets as "target"
  colnames(data_train)[1] = 'target'
  colnames(data_test)[1] = 'target'

  # Create classification tasks for the training and test sets
  traintask <- makeClassifTask(data = data_train, target = "target")
  testtask <- makeClassifTask(data = data_test, target = "target")

  # Perform parameter tuning on the training set
  mytune <- tuneParams(learner = lrn, task = traintask, resampling = rdesc,
                       measures = acc, par.set = params, control = ctrl,
                       show.info = T)

  # Update the learner with the tuned hyperparameters
  lrn_tune <- setHyperPars(lrn, par.vals = mytune$x)

  # Train the model on the training set
  xgmodel <- train(learner = lrn_tune, task = traintask)
  # Make predictions on the test set
  xgpred <- predict(xgmodel, testtask)


  # Stop parallel processing
  parallelStop()

  bestparams = list(booster = "gbtree", objective = "binary:logistic",
                    eval_metric = "auc", eta = xgmodel$learner$par.vals$eta,
                    max_depth = xgmodel$learner$par.vals$max_depth,
                    min_child_weight = xgmodel$learner$par.vals$min_child_weight,
                    subsample = xgmodel$learner$par.vals$subsample,
                    colsample_bytree = xgmodel$learner$par.vals$colsample_bytree)


  dtrain = xgb.DMatrix(data = data.matrix(data_train[,-1]),label=as.numeric(as.character(data_train[,1])))
  xgb = xgb.train(params = bestparams, data = dtrain, nrounds = 1000, maximize = F)

  dtest = xgb.DMatrix(data = data.matrix(data_test[,-1]),label=as.numeric(as.character(data_test[,1])))


  # Plot the ROC curve using the test set labels and model predictions
  plot.roc(data_test[,1], predict(xgb, dtest))

  # Return the trained xgboost model
  return(xgb)

}


evaluate_model = function(model,data) {
  roc(model$data$truth,model$data$prob.1)
  plot.roc(data[[2]][,1],predict(xgb, data[[2]]))

}


################## FUN

customize_labels = function (p, font.title = NULL,
                             font.subtitle = NULL, font.caption = NULL,
                             font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

survivalCurve = function(mymodel = NULL,
                         model_num = NULL,
                         survival_df = NULL,
                         use.TME = TRUE,
                         insert.title = NULL,
                         color1 = NULL,
                         color2 = NULL,
                         split.method = NULL,
                         split.num = NULL) {

  model_predictionscores = mymodel[[model_num]][[1]]
  model_OS = survival_df[rownames(model_predictionscores),]
  rn = rownames(model_OS)[!grepl("^NA", rownames(model_OS))]
  model_OS = na.omit(model_OS)
  rownames(model_OS) = rn
  model_predictionscores = model_predictionscores[rownames(model_OS),]

  df = data.frame(survival_time = as.numeric(model_OS$survival_time),
                  survival_status = as.numeric(model_OS$survival_status),
                  prediction_scores = as.numeric(model_predictionscores),
                  row.names = rownames(model_OS))

  if (use.TME) {
    if (split.method == 'median') {
      median_prediction_score = median(df$prediction_scores)
      df$group = ifelse(df$prediction_scores >= median_prediction_score, 'above_median', 'below_median')

    } else if (split.method == 'tertile') {
      lower_third_threshold = quantile(df$prediction_scores, 0.33)
      upper_third_threshold = quantile(df$prediction_scores, 0.66)
      df$group = with(df, ifelse(prediction_scores <= lower_third_threshold, 'lower_third',
                                 ifelse(prediction_scores > upper_third_threshold, 'upper_third', NA)))
      df = na.omit(df)
    } else {
      num_samples = min(split.num, nrow(df))
      upper_samples = ceiling(num_samples / 2)
      lower_samples = num_samples - upper_samples

      df = df[order(-df$prediction_scores), ]

      # Assign 'upper_third' to the top samples and 'lower_third' to the bottom samples
      df$group = c(rep('upper_group', upper_samples), rep('lower_group', lower_samples), rep(NA, nrow(df) - num_samples))
      df = na.omit(df)
    }

    df$group = as.factor(df$group)
    fit = survfit(Surv(survival_time, survival_status) ~ group, data = df)

  } else {
    # Handle the case when use.TME is FALSE
    # [Add any specific logic if needed]
  }
  ggsurv = ggsurvplot(
    fit,
    data = df,
    risk.table = TRUE,
    pval = TRUE,
    palette = c(color1, color2),
    xlim = c(0, 36),
    break.time.by = 6,
    size = 2
  )

  ggsurv$plot = ggsurv$plot +
    ggtitle(insert.title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 23),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      axis.text.x = element_text(size = 19),
      axis.text.y = element_text(size = 19),
      legend.text = element_text(size = 19),
      axis.text.x.top = element_text(size = 10),
      axis.text.y.right = element_text(size = 10)
    )

  # Customize labels if needed
  # ggsurv$table = customize_labels(
  #   ggsurv$table,
  #   font.title    = c(19, "bold", "darkgreen"),
  #   font.x        = c(16, "bold", "black"),
  #   font.xtickslab = c(16, "bold", "brown")
  # )

  print(ggsurv)
  return(list(df = df, fit = fit, whole_plot = ggsurv, plot = ggsurv$plot, table = ggsurv$table))
}

