library(tidyverse)

totaldata_mini <- readRDS("/bigdata/loainaom/For_Almog/totaldata_mini.rds")
resIMPRES <- readRDS("/bigdata/loainaom/newRuns_0708/RESULTS/PD_oldgenes/new_xCell2/100Runs/newAlgos/auc_scores_pd_impres.rds")

# Set the working directory or specify the full path
path <- "/bigdata/loainaom/newRuns_0708/RESULTS/PD_oldgenes/new_xCell2/100Runs/"

# Manually read each file and assign to the appropriate variable
resxcell1 <- readRDS(paste0(path, "resxcell1.rds"))
resTIDE <- readRDS(paste0(path, "resTIDE.rds"))
resxcell2_allSig_new <- readRDS(paste0(path, "resxcell2_allSig_new.rds"))
resxcell2_dice <- readRDS(paste0(path, "resxcell2_dice.rds"))
resxcell2_kass <- readRDS(paste0(path, "resxcell2_kass.rds"))
resxcell2_pan <- readRDS(paste0(path, "resxcell2_pan.rds"))
resxcell2_bp <- readRDS(paste0(path, "resxcell2_bp.rds"))
resxcell2_lm22 <- readRDS(paste0(path, "resxcell2_lm22.rds"))
resCIBER <- readRDS(paste0(path, "resCIBER.rds"))
resxcell1tide <- readRDS(paste0(path, "resxcell1tide.rds"))
resxcell2Tide <- readRDS(paste0(path, "resxcell2tide_New.rds"))
resxcell2Easier <- readRDS(paste0(path, "resxcell2Easier.rds"))
rescibertide <- readRDS(paste0(path, "resCIBERTide.rds"))
resEASIER <- readRDS(paste0(path, "resEASIER.rds"))
resQuanTIseq <- readRDS(paste0(path, "resQuanTIseq.rds"))
resEasierTide <- readRDS(paste0(path, "resEasierTide.rds"))
resEasierQuanseq <- readRDS(paste0(path, "resEasierQuanseq.rds"))

xcell2new <- readRDS("/bigdata/loainaom/newRuns_0708/RESULTS/PD_oldgenes/new_xCell2/100Runs/newxCell2_foralmog.rds")


pd_lists = c(resxcell1,
             resTIDE,
             resxcell2_allSig_new,
             resxcell2_dice,
             resxcell2_kass,
             resxcell2_pan,
             resxcell2_bp,
             resxcell2_lm22,
             resCIBER,
             resxcell1tide,
             resxcell2Tide,
             resxcell2Easier,
             rescibertide,
             resEASIER,
             resQuanTIseq,
             resEasierTide,
             resEasierQuanseq)


pd_lists <- c(pd_lists, xcell2new)

names(pd_lists)[1:400] <- paste0("xcell1")
names(pd_lists)[401:600] <- paste0("tide")
names(pd_lists)[601:800] <- paste0("xcell2_allSig")
names(pd_lists)[801:1000] <- paste0("xcell2_dice")
names(pd_lists)[1001:1200] <- paste0("xcell2_kass")
names(pd_lists)[1201:1400] <- paste0("xcell2_pan")
names(pd_lists)[1401:1600] <- paste0("xcell2_bp")
names(pd_lists)[1601:1800] <- paste0("xcell2_lm22")
names(pd_lists)[1801:2000] <- paste0("cibersort")
names(pd_lists)[2001:2200] <- paste0("xcell1_tide")
names(pd_lists)[2201:2400] <- paste0("xcell2_tide")
names(pd_lists)[2401:2600] <- paste0("xcell2_easier")
names(pd_lists)[2601:2800] <- paste0("ciber_tide")
names(pd_lists)[2801:3000] <- paste0("EASIER")
names(pd_lists)[3001:3200] <- paste0("QuanTIseq")
names(pd_lists)[3201:3400] <- paste0("Easier_tide")
names(pd_lists)[3401:3600] <- paste0("Easier_Quanseq")
names(pd_lists)[3601:3700] <- paste0("xCell 2.0")

auc_score_saved = list()
auc_scores = data.frame(AUC = numeric(), Model = character())
for (i in 1:length(pd_lists)) {
  pred_scores_df = pd_lists[[i]][[1]]

  # Get the real label values based on the row names of the prediction data frame
  real_labels = totaldata_mini[rownames(pred_scores_df), "PD"]

  auc_score = pROC::auc(pROC::roc(response = real_labels, predictor = as.numeric(pred_scores_df$pred)))
  auc_score_saved[[i]] = auc_score



  if (names(pd_lists[i]) == 'xcell1'){

    if (i %% 4 == 1) {
      model_type = "Metadata"
    } else if (i %% 4 == 2) {
      model_type = "xCell 1.0"
    } else if (i %% 4 == 3) {
      model_type = "Genes"
    } else {
      model_type = "Genes + xCell 1.0"
    }

  } else if (names(pd_lists[i]) == 'tide') {

    if (i %% 2 == 1) {
      model_type = "TIDE"
    } else {
      model_type = "Genes + TIDE"
    }

  } else if (names(pd_lists[i]) == 'xcell2_allSig') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (MultiRef)"
    } else {
      model_type = "Genes + xCell 2.0 (MultiRef)"
    }

  } else if (names(pd_lists[i]) == 'xcell2_dice') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (Dice)"
    } else {
      model_type = "Genes + xCell 2.0 (Dice)"
    }

  } else if (names(pd_lists[i]) == 'xcell2_kass') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (Kass)"
    } else {
      model_type = "Genes + xCell 2.0 (Kass)"
    }

  } else if (names(pd_lists[i]) == 'xcell2_pan') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (Pan Cancer)"
    } else {
      model_type = "Genes + xCell 2.0 (Pan Cancer)"
    }

  } else if (names(pd_lists[i]) == 'xcell2_bp') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (Bp)"
    } else {
      model_type = "Genes + xCell 2.0 (Bp)"
    }

  } else if (names(pd_lists[i]) == 'xcell2_lm22') {

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 (LM22)"
    } else {
      model_type = "Genes + xCell 2.0 (LM22)"
    }

  } else if (names(pd_lists[i]) == 'cibersort') {

    if (i %% 2 == 1) {
      model_type = "CIBERSORTx"
    } else {
      model_type = "Genes + CIBERSORTx"
    }

  } else if (names(pd_lists[i]) == 'xcell1_tide'){

    if (i %% 2 == 1) {
      model_type = "xCell 1.0 + TIDE"
    } else {
      model_type = "Genes + xCell 1.0 + TIDE"
    }

  } else if (names(pd_lists[i]) == 'xcell2_tide'){

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 + TIDE"
    } else {
      model_type = "Genes + xCell 2.0 + TIDE"
    }

  } else if (names(pd_lists[i]) == 'xcell2_easier'){

    if (i %% 2 == 1) {
      model_type = "xCell 2.0 + EASIER"
    } else {
      model_type = "Genes + xCell 2.0 + EASIER"
    }

  } else if (names(pd_lists[i]) == 'ciber_tide') {

    if (i %% 2 == 1) {
      model_type = "CIBERSORTx + TIDE"
    } else {
      model_type = "Genes + CIBERSORTx + TIDE"
    }

  } else if (names(pd_lists[i]) == 'EASIER') {

    if (i %% 2 == 1) {
      model_type = "EASIER"
    } else {
      model_type = "Genes + EASIER"
    }

  } else if (names(pd_lists[i]) == 'QuanTIseq'){

    if (i %% 2 == 1) {
      model_type = "QuanTIseq"
    } else {
      model_type = "Genes + QuanTIseq"
    }

  } else if (names(pd_lists[i]) == 'Easier_tide'){

    if (i %% 2 == 1) {
      model_type = "EASIER + TIDE"
    } else {
      model_type = "Genes + EASIER + TIDE"
    }

  } else if (names(pd_lists[i]) == 'Easier_Quanseq'){

    if (i %% 2 == 1) {
      model_type = "EASIER + QuanTIseq"
    } else {
      model_type = "Genes + EASIER + QuanTIseq"
    }

  }

  if (names(pd_lists[i]) == "xCell 2.0") {
    model_type = "xCell 2.0"
  }

  # Append to auc_scores data frame
  auc_scores = rbind(auc_scores, data.frame(AUC = auc_score, Model = model_type))
}
auc_scores$Model = fct_relevel(auc_scores$Model, "Metadata", "xCell 1.0", "Genes", "Genes + xCell 1.0",
                               "TIDE", "Genes + TIDE",
                               "xCell 2.0 (MultiRef)", "Genes + xCell 2.0 (MultiRef)",
                               "xCell 2.0 (Dice)", "Genes + xCell 2.0 (Dice)",
                               "xCell 2.0 (Kass)", "Genes + xCell 2.0 (Kass)",
                               "xCell 2.0 (Pan Cancer)", "Genes + xCell 2.0 (Pan Cancer)",
                               "xCell 2.0 (Bp)", "Genes + xCell 2.0 (Bp)",
                               "xCell 2.0 (LM22)", "Genes + xCell 2.0 (LM22)",
                               "CIBERSORTx", "Genes + CIBERSORTx",
                               "xCell 1.0 + TIDE", "Genes + xCell 1.0 + TIDE",
                               "xCell 2.0 + TIDE", "Genes + xCell 2.0 + TIDE",
                               "xCell 2.0 + EASIER", "Genes + xCell 2.0 + EASIER",
                               "CIBERSORTx + TIDE", "Genes + CIBERSORTx + TIDE",
                               "EASIER", "Genes + EASIER",
                               "QuanTIseq", "Genes + QuanTIseq",
                               "EASIER + TIDE", "Genes + EASIER + TIDE",
                               "EASIER + QuanTIseq", "Genes + EASIER + QuanTIseq",
                               "xCell 2.0")
medians = with(auc_scores, tapply(AUC, Model, median, na.rm = TRUE))
auc_scores$Model = with(auc_scores, reorder(Model, AUC, FUN = median))
color_vec = c("Metadata"="gray",
              "xCell 1.0"="gray",
              "Genes"="gray",
              "Genes + xCell 1.0"="gray",
              "TIDE"="gray",
              "Genes + TIDE"="gray",
              "xCell 2.0 (MultiRef)"="gray",
              "Genes + xCell 2.0 (MultiRef)"="gray",
              "xCell 2.0 (Dice)"="gray",
              "Genes + xCell 2.0 (Dice)"="gray",
              "xCell 2.0 (Kass)"="gray",
              "Genes + xCell 2.0 (Kass)"="gray",
              "xCell 2.0 (Pan Cancer)"="gray",
              "Genes + xCell 2.0 (Pan Cancer)"="gray",
              "xCell 2.0 (Bp)"="gray",
              "Genes + xCell 2.0 (Bp)"="gray",
              "xCell 2.0 (LM22)"="gray",
              "Genes + xCell 2.0 (LM22)"="gray",
              "CIBERSORTx"="gray",
              "Genes + CIBERSORTx"="gray",
              "xCell 1.0 + TIDE"="gray",
              "Genes + xCell 1.0 + TIDE"="gray",
              "xCell 2.0 + TIDE"="gray",
              "Genes + xCell 2.0 + TIDE"="gray",
              "xCell 2.0 + EASIER" = "gray",
              "Genes + xCell 2.0 + EASIER" = "gray",
              "CIBERSORTx + TIDE"="gray",
              "Genes + CIBERSORTx + TIDE"="gray",
              "EASIER" = "gray",
              "Genes + EASIER" = "gray",
              "QuanTIseq" = "gray",
              "Genes + QuanTIseq" = "gray",
              "EASIER + TIDE" = "gray",
              "Genes + EASIER + TIDE" = "gray",
              "EASIER + QuanTIseq" = "gray",
              "Genes + EASIER + QuanTIseq" = "gray",
              "xCell 2.0" = "tomato")
resIMPRES_scores = as.numeric(resIMPRES) # Your provided AUC scores
impres_df = data.frame(AUC = resIMPRES_scores, Model = rep("IMPRES", length(resIMPRES_scores)))
auc_scores = rbind(auc_scores, impres_df)
auc_scores$Model = factor(auc_scores$Model, levels = unique(auc_scores$Model))
color_vec["IMPRES"] = "gray"
model_medians <- with(auc_scores, tapply(AUC, Model, median, na.rm = TRUE))
# Order the models by median AUC
ordered_models <- names(sort(model_medians))
# Set the levels of the Model factor in the auc_scores data frame to this ordering
auc_scores$Model <- factor(auc_scores$Model, levels = ordered_models)
# Update the color vector if necessary (make sure it includes all necessary model types)


# Plot the updated results
results <- ggplot(auc_scores, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot() +
  labs(title = "Benchmarking Predictive Models of Progressive Disease",
       x = "Model Type",
       y = "AUC Score") +
  theme_bw() +
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 19),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_fill_manual(values = color_vec)
results


d <- auc_scores %>%
  filter(Model %in% c("Metadata", "xCell 1.0", "TIDE", "CIBERSORTx", "IMPRES", "EASIER", "xCell 2.0"))


d %>%
  ggplot(., aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "AUC") +
  theme_minimal() +
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 19),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_fill_manual(values = color_vec) +
  ggpubr::stat_compare_means(comparisons = list(c("TIDE", "xCell 2.0")),
                             method = "wilcox.test", label = "p.signif", paired = FALSE)





# Survival analysis ----

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

# Define models and related parameters
models_data = list(  # 9 , 16, 17, 20, 23, 31, 39, 53, 55,64, 65, 68
  xCell2_new = list(data = xcell2new[[1]], col_suffix = "xCell 2.0", title = 'xCell2-New', model_num = 1)
)


# Create survival curves, rename columns, and filter data
list_of_1 = list()
list_of_0 = list()
for (model in names(models_data)) {
  survival_data = survivalCurve(
    mymodel = models_data[[model]]$data,
    model_num = models_data[[model]]$model_num,
    use.TME = ifelse(is.null(models_data[[model]]$use.TME), TRUE, models_data[[model]]$use.TME),
    split.method  = 'num',
    split.num = 413,
    survival_df = overallSurvival,
    color1 = 'brown1',
    color2 = "darkturquoise",
    insert.title = models_data[[model]]$title
  )

  colname_suffix = paste0("group_", models_data[[model]]$col_suffix)
  colnames(survival_data[[1]]) = paste0(colnames(survival_data[[1]]), "_", models_data[[model]]$col_suffix)

  list_of_1[[model]] <- survival_data[[1]] %>% dplyr::filter(get(colname_suffix) == 'upper_group') %>% dplyr::mutate(rowname = rownames(.))
  list_of_0[[model]] <- survival_data[[1]] %>% dplyr::filter(get(colname_suffix) == 'lower_group') %>% dplyr::mutate(rowname = rownames(.))
}
# Merge the data
ones = Reduce(function(x, y) merge(x, y, by="rowname", all=TRUE), list_of_1)
zeros = Reduce(function(x, y) merge(x, y, by="rowname", all=TRUE), list_of_0)
df = ones
# Compute the survival curves
fit_tide = survfit(Surv(survival_time_TIDE, survival_status_TIDE) ~ 1, data = df)
fit_ciber = survfit(Surv(survival_time_CIBERSORT, survival_status_CIBERSORT) ~ 1, data = df)
fit_xcell2 = survfit(Surv(`survival_time_xCell 2.0`, `survival_status_xCell 2.0`) ~ 1, data = df)
fit_xcell1 = survfit(Surv(`survival_time_xCell 1.0`, `survival_status_xCell 1.0`) ~ 1, data = df)
fit_baseline = survfit(Surv(survival_time_Metadata, survival_status_Metadata) ~ 1, data = df)
fit_genes = survfit(Surv(survival_time_Genes, survival_status_Genes) ~ 1, data = df)
fit_easier = survfit(Surv(survival_time_EASIER, survival_status_EASIER) ~ 1, data = df)
fit_quan = survfit(Surv(survival_time_quanTIseq, survival_status_quanTIseq) ~ 1, data = df)
p = ggsurvplot_combine(
  list(fit_tide, fit_ciber,fit_xcell2,fit_genes,fit_xcell1,fit_baseline, fit_easier, fit_quan),
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.labs = c("TIDE", "CIBERSORTx","xCell 2.0","Genes","xCell 1.0","Baseline","EASIER","quanTIseq"),
  palette = c("forestgreen", "darkorange2","black", "blue3", "gold2", "firebrick3","dodgerblue4","purple3"),
  size = 3,
  font.main = list(size = 34),          # Main title font size
  font.submain = list(size = 33),       # Subtitle font size
  font.x = list(size = 30),             # x-axis label font size
  font.y = list(size = 30),             # y-axis label font size
  font.tickslab = list(size = 27),      # Tick label font size
  font.legend = list(size = 29)         # Legend font size
)
p








