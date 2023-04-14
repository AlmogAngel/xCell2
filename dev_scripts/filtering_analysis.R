library(tidyverse)





makeMixture2 <- function(ctoi, ct_mat = pure_ct_mat, dep = dep_list, fractions = mixture_fractions){


  ctoi_fracs <- sapply(fractions, function(frac){
    ct_mat[,ctoi] * frac
  })

  no_dep_controls <- Rfast::rowMedians(ct_mat[,!colnames(ct_mat) %in% dep[[ctoi]]])

  controls_fracs <- sapply(fractions, function(frac){
    no_dep_controls * (1-frac)
  })

  mix <- ctoi_fracs + controls_fracs
  colnames(mix) <- fractions

  return(mix)
}

ctoi <- "CD8+ T-cells"

# Signature CT ----
all_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct == sample_ct) %>%
  filter(signature %in% grubbs) %>%
  drop_na() %>%
  pull(score)

top_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct == sample_ct & signature %in% top_sigs) %>%
  drop_na() %>%
  pull(score)

df <- data.frame(sigs = c(rep("All", length(all_sigs_scores)), rep("Top", length(top_sigs_scores))),
                          score = c(all_sigs_scores, top_sigs_scores))
ggplot(df, aes(sigs, score, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sigs, score, col=sigs), alpha = 0.5)


# Not signature CT ----

all_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct != sample_ct) %>%
  filter(signature %in% grubbs) %>%
  drop_na() %>%
  group_by(signature) %>%
  summarise(avg_score = mean(score)) %>%
  mutate(sigs = "All")

top_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct != sample_ct & signature %in% top_sigs) %>%
  drop_na() %>%
  group_by(signature) %>%
  summarise(avg_score = mean(score)) %>%
  mutate(sigs = "Top")

rbind(all_sigs_scores, top_sigs_scores) %>%
  ggplot(., aes(sigs, avg_score, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sigs, avg_score, col=sigs), alpha = 0.5)

# Correlation of scores with similar cell types ----
all_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct != sample_ct) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(similarity = cor_mat[signature_ct, sample_ct]) %>%
  select(signature, score, similarity) %>%
  group_by(signature) %>%
  summarise(cor = cor(score, similarity, method = "pearson")) %>%
  mutate(sigs = "All")


top_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct != sample_ct & signature %in% top_sigs[1:20]) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(similarity = cor_mat[signature_ct, sample_ct]) %>%
  select(signature, score, similarity) %>%
  group_by(signature) %>%
  summarise(cor = cor(score, similarity, method = "pearson")) %>%
  mutate(sigs = "Top")


rbind(all_sigs_scores, top_sigs_scores) %>%
  ggplot(., aes(sigs, cor, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2)

# Correlation with CTOI reference mixtures ----
makeMixture2 <- function(ctoi = "B-cells", ct_mat = pure_ct_mat, dep = dep_list, fractions = mixture_fractions){


  ctoi_fracs <- sapply(fractions, function(frac){
    ct_mat[,ctoi] * frac
  })

  no_dep_controls <- Rfast::rowMedians(ct_mat[,!colnames(ct_mat) %in% dep[[ctoi]]])

  controls_fracs <- sapply(fractions, function(frac){
    no_dep_controls * (1-frac)
  })

  mix <- ctoi_fracs + controls_fracs
  colnames(mix) <- fractions

  return(mix)
}


b_mix <- makeMixture2(ctoi = "B-cells")
signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), ctoi)]


mix_ranked <- singscore::rankGenes(b_mix)
scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
  singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
})
colnames(scores) <- names(signatures_ctoi)
rownames(scores) <- colnames(mix_ranked)

cors <- apply(scores, 2, function(x){
  cor(as.numeric(rownames(scores)), x, method = "pearson")
})

cors_top <- cors[names(cors) %in% top_sigs[1:20]]

df <- data.frame(sigs = c(rep("All", length(cors)), rep("Top", length(cors_top))),
                 score = c(cors, cors_top))
ggplot(df, aes(sigs, score, fill=sigs)) +
  geom_boxplot() +
  geom_jitter()

# Correlation with not CTOI reference mixtures ----

signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), ctoi)]

all_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == "B-cells" & signature_ct != sample_ct) %>%
  drop_na() %>%
  select(sample_ct) %>%
  unique() %>%
  rowwise() %>%
  mutate(mix = list(makeMixture2(ctoi = sample_ct))) %>%
  mutate(scores = list(sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(singscore::rankGenes(mix), upSet = sig, centerScore = FALSE)$TotalScore
  })))

all_sigs_scores <- all_sigs_scores %>%
  mutate(cors = list(apply(scores, 2, function(x){
    cor(mixture_fractions, x, method = "pearson")
  }))) %>%
  select(sample_ct, cors) %>%
  unnest(cols = c(cors)) %>%
  mutate(signature = rep(names(signatures_ctoi), length(unique(all_sigs_scores$sample_ct))))

all_sigs_scores %>%
  filter(signature %in% grubbs) %>%
  rowwise() %>%
  mutate(sigs = ifelse(signature %in% top_sigs, "Top", "All")) %>%
  ggplot(., aes(sample_ct, cors, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sample_ct, cors, col=sigs), alpha = 0.5)

all_sigs_scores %>%
  filter(signature %in% grubbs) %>%
  group_by(signature) %>%
  summarise(avg_cor = mean(cors)) %>%
  rowwise() %>%
  mutate(sigs = ifelse(signature %in% top_sigs, "Top", "All")) %>%
  ggplot(., aes(sigs, avg_cor, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sigs, avg_cor, col=sigs), alpha = 0.5)


all_sigs_scores2 <- all_sigs_scores %>%
  filter(sample_ct != "B-cells" & signature %in% grubbs) %>%
  group_by(signature) %>%
  summarise(avg_cor = mean(cors)) %>%
  mutate(sigs = "All")

#all_sigs_scores2 <- all_sigs_scores %>%
#  filter(sample_ct != "B-cells") %>%
#  drop_na() %>%
#  select(sample_ct, cors) %>%
#  mutate(sigs = "All")

#top_sigs_scores <- all_sigs_scores %>%
#  filter(sample_ct != "B-cells" & signature %in% top_sigs[1:20]) %>%
#  drop_na() %>%
#  select(sample_ct, cors) %>%
#  mutate(sigs = "Top")

top_sigs_scores <- all_sigs_scores %>%
  filter(sample_ct != "B-cells" & signature %in% top_sigs[1:20]) %>%
  group_by(signature) %>%
  summarise(avg_cor = mean(cors)) %>%
  mutate(sigs = "Top")

rbind(all_sigs_scores2, top_sigs_scores) %>%
  ggplot(., aes(sigs, avg_cor, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sigs, avg_cor, col=sigs), alpha = 0.5)

#  Grubbs' test  ----

all_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == ctoi) %>%
  drop_na() %>%
  #rowwise() %>%
  #mutate(similarity = cor_mat[signature_ct, sample_ct]) %>%
  #filter(similarity > 0.9) %>%
  group_by(signature) %>%
  summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 20, opposite = FALSE, two.sided = FALSE)$p.value) %>%
  mutate(sigs = "All")

top_sigs_scores <- scores_mat_tidy_blood %>%
  filter(signature_ct == ctoi & signature %in% top_sigs) %>%
  drop_na() %>%
  #rowwise() %>%
  #mutate(similarity = cor_mat[signature_ct, sample_ct]) %>%
  #filter(similarity > 0.9) %>%
  group_by(signature) %>%
  summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 20, opposite = FALSE, two.sided = FALSE)$p.value) %>%
  mutate(sigs = "Top")

rbind(all_sigs_scores, top_sigs_scores) %>%
  filter(grubbs_pvalue < 0.05) %>%
  ggplot(., aes(sigs, grubbs_pvalue, fill=sigs)) +
  geom_boxplot() +
  geom_jitter(aes(sigs, grubbs_pvalue, col=sigs), alpha = 0.5)



# Filtering! -----

grubbs <- scores_mat_tidy_blood %>%
  filter(signature_ct == ctoi) %>%
  drop_na() %>%
  group_by(signature) %>%
  summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 20, opposite = FALSE, two.sided = FALSE)$p.value) %>%
  filter(if (n() >= 5) grubbs_pvalue <= 0.02 else grubbs_pvalue < 1) %>%
  pull(signature)

all(top_sigs %in% grubbs)
which(top_sigs %in% grubbs)


signatures_ctoi.grubbs <- signatures_ctoi[names(signatures_ctoi) %in% grubbs]

all_genes <- unlist(lapply(names(signatures_ctoi.grubbs), function(s){
  GSEABase::geneIds(signatures_ctoi.grubbs[[s]])
}))






top_markers <- names(sort(table(all_genes), decreasing = T)[1:50])

sigs2check <- unlist(lapply(names(signatures_ctoi.grubbs), function(s){
  sig_genes <- GSEABase::geneIds(signatures_ctoi.grubbs[[s]])
  if (!all(top_markers %in% sig_genes)) {
    s
  }
}))

which(top_sigs %in% sigs2check)


