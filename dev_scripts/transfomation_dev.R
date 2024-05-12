
signatures_ctoi <- signatures[startsWith(names(signatures), paste0(ctoi, "#"))]


# Mixture ------------

# Score mixture - all sigs
mix_ranked <- singscore::rankGenes(mix)
scores_ctoi_mix <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
})
scores_ctoi_mix <- cbind(scores_ctoi_mix, frac = truth_mat[ctoi, colnames(mix)], sample = colnames(mix))
scores_ctoi_mix <- as_tibble(scores_ctoi_mix) %>%
  pivot_longer(-c(frac, sample), names_to = "sig", values_to = "score") %>%
  mutate(frac = as.numeric(frac),
         score = as.numeric(score),
         type = "Mixture")


mix_data <- scores_ctoi_mix %>%
  group_by(type, sample, frac) %>%
  summarise(score = mean(score)) %>%
  group_by(type) %>%
  mutate(score = score-min(score)) %>%
  ungroup()


plot(mix_data$score, mix_data$frac)
cor(mix_data$score, mix_data$frac, method = "s")

# Learn linear transformation parameters
scores.in <- apply(scores_ctoi_mix, 2, function(x){
  x-min(x)
})

pred <- round(predict(model, scores.in, type = "response"), 8)
cor(pred, fracs, method = "s")
plot(pred, fracs)

# Filtering data ------------

ds2use <- unlist(lapply(filtering_data$truth, function(x){
  ctoi %in% rownames(x)
}))

filt_ds_sigs <- parallel::mclapply(names(filtering_data$mixture[ds2use]), function(ds){

  ctoi_filt_data <- filtering_data$mixture[[ds]]
  fracs <- filtering_data$truth[[ds]][ctoi,]
  fracs <- fracs[!is.na(fracs)]
  samples <- intersect(names(fracs), colnames(ctoi_filt_data))
  fracs <- fracs[samples]
  fracs <- as.numeric(fracs)

  filt_ranked <- singscore::rankGenes(ctoi_filt_data[,samples])
  scores_ctoi_filt <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    suppressWarnings(singscore::simpleScore(filt_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
  })
  scores_ctoi_filt <- cbind(scores_ctoi_filt, frac = fracs, sample = samples)
  scores_ctoi_filt <- as_tibble(scores_ctoi_filt) %>%
    pivot_longer(-c(frac, sample), names_to = "sig", values_to = "score") %>%
    mutate(frac = as.numeric(frac),
           score = as.numeric(score),
           type = paste0(ds, "-Filtering"))

  scores_ctoi_filt

}, mc.cores = 25)
filt_ds_sigs <- bind_rows(filt_ds_sigs)


# Learn linear transformation parameters
filt2parm <- filt_ds_sigs %>%
  group_by(type) %>%
  summarise(tp = list(nls(score ~ a * frac^b, start = list(a=1, b=1), control = list(maxiter = 500)))) %>%
  rowwise() %>%
  mutate(a = coef(tp)[[1]][1],
         b = coef(tp)[[2]][1])

# # Perform linear transformation parameters
# scores_ctoi_filt_transfomed <- scores_ctoi_filt %>%
#   left_join(sig2parm, by = "sig") %>%
#   rowwise() %>%
#   mutate(score = (score^(1/b)) / a)
#
# # Find filtering data shift value and slope
# filt2shift2slope <- scores_ctoi_filt_transfomed %>%
#   group_by(sig) %>%
#   summarise(lm = list(lm(score~frac))) %>%
#   rowwise() %>%
#   mutate(filt_shift = coef(lm)[[1]],
#          filt_slope = coef(lm)[[2]]) %>%
#   ungroup()


# plot(scores_ctoi_filt_transfomed[scores_ctoi_filt_transfomed$sig == "B cell#_0.05_2.322_12",]$score,
#      scores_ctoi_filt_transfomed[scores_ctoi_filt_transfomed$sig == "B cell#_0.05_2.322_12",]$frac)
# plot(scores_ctoi_filt[scores_ctoi_filt$sig == "B cell#_0.05_2.322_12",]$score,
#      scores_ctoi_filt[scores_ctoi_filt$sig == "B cell#_0.05_2.322_12",]$frac)

# Simulations ----

# Learn linear transformation parameters


sim2params <- simulations_scored[[ctoi]] %>%
  group_by(sim, frac) %>%
  summarise(score =  mean(score)) %>%
  mutate(score = score-min(score)) %>%
  filter(frac > 0.01) %>%
  summarise(tp = list(stats::nls(score ~ a * frac^b, start = list(a=1, b=1), control = list(maxiter = 500)))) %>%
  rowwise() %>%
  mutate(a = coef(tp)[[1]][1],
         b = coef(tp)[[2]][1]) %>%
  select(sim, a, b) %>%
  ungroup()

# Perform linear transformation on simulation
scores_ctoi_sim_transfomed <- simulations_scored[[ctoi]]  %>%
  left_join(sig2parm.sim, by = c("sim", "sig")) %>%
  rowwise() %>%
  mutate(score = (score^(1/b)) / a) %>%
  ungroup()

# Find simulation data shift value and slope
sim2shift2slope <- scores_ctoi_sim_transfomed %>%
  group_by(sim, sig) %>%
  summarise(lm = list(lm(score~frac))) %>%
  rowwise() %>%
  mutate(shift = coef(lm)[[1]],
         slope = coef(lm)[[2]]) %>%
  ungroup()


# plot(sim.lm.tbl[sim.lm.tbl$sig == "B cell#_0.05_2.322_12" & sim.lm.tbl$sim == "sim-1",]$score,
#      sim.lm.tbl[sim.lm.tbl$sig == "B cell#_0.05_2.322_12" & sim.lm.tbl$sim == "sim-1",]$frac)
# plot(scores_ctoi_sim_transfomed[scores_ctoi_sim_transfomed$sig == "B cell#_0.05_2.322_12" & scores_ctoi_sim_transfomed$sim == "sim-1",]$score,
#      scores_ctoi_sim_transfomed[scores_ctoi_sim_transfomed$sig == "B cell#_0.05_2.322_12" & scores_ctoi_sim_transfomed$sim == "sim-1",]$frac)

# Perform linear transformation on filtering data and calculate shift and slope
filt2sim2sig2shift2slope <- lapply(unique(sig2parm.sim$sim), function(sim_id){

  sim_param <- scores_ctoi_sim_transfomed %>%
    filter(sim == sim_id) %>%
    select(sig, a, b) %>%
    unique()

  scores_ctoi_filt %>%
    left_join(sim_param, by = "sig") %>%
    rowwise() %>%
    mutate(score_filt = (score_filt^(1/b)) / a) %>%
    group_by(sig) %>%
    summarise(lm = list(lm(score_filt~frac_filt))) %>%
    rowwise() %>%
    mutate(filt_shift = coef(lm)[[1]],
           filt_slope = coef(lm)[[2]]) %>%
    ungroup() %>%
    mutate(sim = sim_id) %>%
    select(sim, sig, filt_shift, filt_slope)

})
filt2sim2sig2shift2slope <- bind_rows(filt2sim2sig2shift2slope)


# Choose sim who have the closest shift value and slope to the filtering data
best_sim <- sim2shift2slope %>%
  left_join(filt2sim2sig2shift2slope, by = c("sim", "sig")) %>%
  mutate(shift_dist = abs(shift - filt_shift),
         slope_dist = abs(slope - filt_slope)) %>%
  ungroup() %>%
  mutate(shift_ranked = rank(-shift_dist),
         slope_ranked = rank(-slope_dist)) %>%
  mutate(sim_sig_score = shift_ranked + slope_ranked) %>%
  group_by(sim) %>%
  summarise(sim_score = mean(sim_sig_score)) %>%
  arrange(-sim_score) %>%
  slice(1) %>%
  pull(sim)

sim2shift2slope %>%
  left_join(filt2sim2sig2shift2slope, by = c("sim", "sig")) %>%
  mutate(shift_dist = abs(shift - filt_shift),
         slope_dist = abs(slope - filt_slope)) %>%
  filter(sim == best_sim)

# Get shift value
sig2shift <- sim2shift2slope %>%
  filter(sim == best_sim) %>%
  select(sig, shift)

# Get transformation parameters
sig2param <- scores_ctoi_sim_transfomed %>%
  filter(sim == best_sim) %>%
  select(sig, a, b) %>%
  unique() %>%
  left_join(sig2shift, by = "sig")



# Prepare all data -----

# Mixture
mix_data <- scores_ctoi_mix %>%
  group_by(type, sample, frac) %>%
  summarise(score = mean(score)) %>%
  group_by(type) %>%
  mutate(score = score-min(score))

#  left_join(sig2param, by = "sig") %>%
#  rowwise() %>%
#  mutate(score = (score^(1/b)) / a) %>%
#  mutate(score = ifelse(score < 0, 0, score)) %>%
#  select(-c(a, b, shift)) %>%
#  ungroup()

# Filtering
filt_data <- filt_ds_sigs %>%
  group_by(type, sample, frac) %>%
  summarise(score = mean(score)) %>%
  group_by(type) %>%
  mutate(score = score-min(score))

# Simulation
sim_data <- scores_ctoi_sim_transfomed %>%
  filter(sim == best_sim) %>%
  filter(frac >= 0.01) %>%
  select(frac, sample, sig, score, type) %>%
  left_join(sig2param, by = "sig") %>%
  mutate(score = score-shift) %>%
  mutate(score = ifelse(score < 0, 0, score)) %>%
  select(frac, sample, sig, score, type)


# Plot ----

data_combined <- rbind(mix_data, filt_data)



data_combined %>%
  filter(type == "Mixture") %>%
  summarise(cor=cor(score, frac, method = "s"))
  ungroup() %>%
  ggplot(., aes(x=frac, y=score, col=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)




data_combined %>%
  filter(sig == "B cell#_0.05_1.322_12") %>%
  #group_by(type, sample, frac) %>%
  #summarise(score = mean(score)) %>%
  #ungroup() %>%
  ggplot(., aes(x=frac, y=score, col=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# New - use xCell2's outputs -----

mix_data <- scores_ctoi_mix %>%
  left_join(select(params[[ctoi]], sig, a, b, shift), by = "sig") %>%
  rowwise() %>%
  mutate(score = (score^(1/b)) / a) %>%
  mutate(score = score-shift) %>%
  mutate(score = ifelse(score < 0, 0, score)) %>%
  select(-c(a, b, shift)) %>%
  ungroup()


filt_data <- scores_ctoi_filt %>%
  left_join(select(params[[ctoi]], sig, a, b, shift), by = "sig") %>%
  rowwise() %>%
  mutate(score = (score_filt^(1/b)) / a) %>%
  mutate(score = score-shift) %>%
  mutate(score = ifelse(score < 0, 0, score)) %>%
  mutate(frac = frac_filt) %>%
  select(frac, sample, sig, score, type) %>%
  ungroup()

sim_data <- simulations_scored[[ctoi]] %>%
  filter(sim == params[[ctoi]]$sim) %>%
  left_join(select(params[[ctoi]], sig, a, b, shift), by = "sig") %>%
  rowwise() %>%
  mutate(score = (score^(1/b)) / a) %>%
  mutate(score = score-shift) %>%
  mutate(score = ifelse(score < 0, 0, score)) %>%
  mutate(type = "Simulation",
         sample = paste0(sim, "-", frac)) %>%
  select(frac, sample, sig, score, type) %>%
  ungroup()


data_combined <- rbind(mix_data, filt_data, sim_data)


data_combined %>%
  filter(sig == "monocyte#_0.05_1.585_7") %>%
  #group_by(type, sample, frac) %>%
  #summarise(score = mean(score)) %>%
  #ungroup() %>%
  ggplot(., aes(x=frac, y=score, col=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)





















scores_ctoi_sim <- scores_ctoi_sim %>%
  rowwise() %>%
  mutate(score = score - shift_value) %>%
  mutate(score = ifelse(score < 0 , 0, score)) %>%
  mutate(score = (score^(1/coef(xx)[[2]])) / coef(xx)[[1]] )

scores_ctoi_mix <- scores_ctoi_mix %>%
  rowwise() %>%
  mutate(score = score - shift_value)




# scores_combined <- rbind(scores_ctoi_mix, scores_ctoi_sim, scores_ctoi_filt)
# scores_combined <- rbind(scores_ctoi_mix, scores_ctoi_filt)
scores_combined <- rbind(scores_ctoi_mix, scores_ctoi_sim)
# scores_combined <- scores_ctoi_mix

scores_combined %>%
  group_by(type, sample, frac) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  #mutate(score = score - shift_value) %>%
  #mutate(score = ifelse(score < 0 , 0, score)) %>%
  #mutate(score = (score^(1/coef(xx)[[2]])) / coef(xx)[[1]] ) %>%
  #filter(frac <= 0.15) %>%
  ggplot(., aes(x=frac, y=score, col=type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

