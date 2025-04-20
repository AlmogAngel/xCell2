library(tidyverse)


source("/bigdata/almogangel/xCell2_dev/paper_figures/load_benchmark_data.R")


# Human

xCell2results <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig1D/fig1d_xcell2_res_alpha_0.5.rds")


benchmark_correlations_spr_raw <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = FALSE, cMethod = "spearman")
colnames(benchmark_correlations_spr_raw)[3:4] <- c("cor_spearman", "p_spearman")

benchmark_correlations_prs_raw <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = FALSE, cMethod = "pearson")
colnames(benchmark_correlations_prs_raw)[3:4] <- c("cor_pearson", "p_pearson")

benchmark_correlations_spr_raw |> 
  left_join(benchmark_correlations_prs_raw) |> 
  select(method, celltype, cor_spearman, p_spearman, cor_pearson, p_pearson, n, ref, val) |> 
  writexl::write_xlsx("/bigdata/almogangel/xCell2_dev/paper_figures/supp/benchmark_correlations_human.xlsx")


benchmark_correlations_spr_weighted <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "spearman")
benchmark_correlations_spr_weighted <- benchmark_correlations_spr_weighted |> 
  group_by(method, ref) |> 
  summarise(ref_cor_mean = mean(ref_cor)) |> 
  rename(mean_spearman = ref_cor_mean)


benchmark_correlations_prs_weighted  <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "pearson")
benchmark_correlations_prs_weighted <- benchmark_correlations_prs_weighted |> 
  group_by(method, ref) |> 
  summarise(ref_cor_mean = mean(ref_cor)) |> 
  rename(mean_pearson = ref_cor_mean)


benchmark_correlations_spr_weighted |> 
  left_join(benchmark_correlations_prs_weighted) |> 
  writexl::write_xlsx("/bigdata/almogangel/xCell2_dev/paper_figures/supp/benchmark_mean_correlations_human.xlsx")


# Mouse

xCell2results <- readRDS("/bigdata/almogangel/xCell2_dev/paper_figures/data/fig2C/fig2c_xcell2_mouse_res.rds")

benchmark_correlations_mouse_spr <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "spearman")
benchmark_correlations_mouse_prs <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "pearson")

# Save for supplementary
benchmark_correlations_mouse_spr_raw <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = FALSE, cMethod = "spearman")
colnames(benchmark_correlations_mouse_spr_raw)[3:4] <- c("cor_spearman", "p_spearman")
benchmark_correlations_mouse_prs_raw <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = FALSE, cMethod = "pearson")
colnames(benchmark_correlations_mouse_prs_raw)[3:4] <- c("cor_pearson", "p_pearson")

benchmark_correlations_mouse_spr_raw |> 
  left_join(benchmark_correlations_mouse_prs_raw) |> 
  select(method, celltype, cor_spearman, p_spearman, cor_pearson, p_pearson, n, ref, val) |> 
  writexl::write_xlsx("/bigdata/almogangel/xCell2_dev/paper_figures/benchmark_correlations_mouse.xlsx")


benchmark_correlations_spr_weighted <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "spearman")
benchmark_correlations_spr_weighted <- benchmark_correlations_spr_weighted |> 
  group_by(method, ref) |> 
  summarise(ref_cor_mean = mean(ref_cor)) |> 
  rename(mean_spearman = ref_cor_mean)


benchmark_correlations_prs_weighted  <- get_xcell2_correlations(xCell2results = xCell2results, round_results = 3, weight_cors = TRUE, cMethod = "pearson")
benchmark_correlations_prs_weighted <- benchmark_correlations_prs_weighted |> 
  group_by(method, ref) |> 
  summarise(ref_cor_mean = mean(ref_cor)) |> 
  rename(mean_pearson = ref_cor_mean)


benchmark_correlations_spr_weighted |> 
  left_join(benchmark_correlations_prs_weighted) |> 
  writexl::write_xlsx("/bigdata/almogangel/xCell2_dev/paper_figures/supp/benchmark_mean_correlations_mouse.xlsx")