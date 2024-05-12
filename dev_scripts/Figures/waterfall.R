library(tidyverse)

data.in = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_sigs_10apr/withEssen/xcell2.filteredSigsCors.10apr.rds")
data.in = bind_rows(data.in)


unique(data.in$is_filtered)

# data.in2 <- data.in %>%
#   filter(is_filtered != "Filtered")

data <- data.in %>%
  group_by(ref, val, celltype, is_filtered) %>%
  summarise(cor = mean(cor)) %>%
  # group_by(ref, val, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  # group_by(ref, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  spread(key = is_filtered, value = cor) %>%
  # mutate(delta_cor = `Filtered + Essential Genes` - `Unfiltered`) %>%
  mutate(delta_cor = `yes` - `no`) %>%
  #select(ref, val, celltype, delta_cor) %>%
  #select(ref, val, delta_cor) %>%
  select(ref, delta_cor) %>%
  drop_na()


data_sorted <- data %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  mutate(x_axis = row_number())

# Create a waterfall plot
ggplot(data_sorted, aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "darkred", high = "darkblue", mid = "yellow", midpoint = 0) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Delta Cor", title = "Waterfall Plot of Delta Cor")



# Ontology


data.no = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_sigs_10apr/noOnto/xcell2.sigs.noOnto.10apr.rds")
data.no = bind_rows(data.no)
data.no$with_onto <- "no"

data.yes = readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/xcell2_sigs_10apr/withOnto/xcell2.sigs.withOnto.10apr.rds")
data.yes = bind_rows(data.yes)

data.in = rbind(data.yes, data.no)


data <- data.in %>%
  group_by(ref, val, celltype, with_onto) %>%
  summarise(cor = mean(cor)) %>%
  # group_by(ref, val, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  # group_by(ref, is_filtered) %>%
  # summarise(cor = mean(cor)) %>%
  spread(key = with_onto, value = cor) %>%
  # mutate(delta_cor = `Filtered + Essential Genes` - `Unfiltered`) %>%
  mutate(delta_cor = `yes` - `no`) %>%
  #select(ref, val, celltype, delta_cor) %>%
  #select(ref, val, delta_cor) %>%
  select(ref, delta_cor) %>%
  drop_na()


data_sorted <- data %>%
  ungroup() %>%
  arrange(delta_cor) %>%
  mutate(x_axis = row_number())

red_color_palette <- colorRampPalette(c("darkred", "red"))(11)
blue_color_palette <- colorRampPalette(c("lightblue1", "darkblue"))(36)
color_palette <- c(red_color_palette, blue_color_palette)

# Create a waterfall plot
ggplot(data_sorted, aes(x = x_axis, y = delta_cor, fill = delta_cor)) +
  geom_bar(stat = "identity") +
  #scale_fill_gradientn(low = "darkred", high = "darkblue", midpoint = 0, guide = FALSE,space = "Lab",  limits = c(min(data_sorted$delta_cor), max(data_sorted$delta_cor))) +
  scale_fill_gradientn(colors = color_palette, guide = FALSE, ) +
  geom_hline(yintercept = mean(data_sorted$delta_cor), linetype = "dashed") +
  annotate("text", x = 100, y = mean(data_sorted$delta_cor)+0.005, label = paste0("Mean delta: ", round(mean(data_sorted$delta_cor), 2)), vjust = -0.5) +
  theme(axis.title.x = element_blank()) +
  labs(x = "", y = "Delta Correlation", title = "") +
  scale_y_continuous(limits = c(-0.08, 0.23), breaks = seq(-0.08, 0.23, 0.04)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank())












data_sorted %>%
  ggplot(., aes(x=reorder(celltype, delta_cor, median) , y=delta_cor)) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

data_sorted %>%
  filter(celltype == "myeloid cell")

data.in %>%
  group_by(ref, val, celltype, is_filtered) %>%
  summarise(cor = median(cor)) %>%
  group_by(ref, val, is_filtered) %>%
  summarise(cor = median(cor)) %>%
  group_by(ref, is_filtered) %>%
  summarise(cor = median(cor)) %>%
  ggplot(., aes(x=is_filtered, y=cor, fill=is_filtered)) +
  geom_boxplot() +
  geom_jitter()


y <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.SigsCorsPca.8mar.rds")
y=bind_rows(y)

library(RColorBrewer)
display.brewer.all()
palette <- brewer.pal(n = 10, name = "Set3")

jitter_position <- position_jitter(width = 0.1, height = 0)

y %>%
  filter(val == "SDY420" & ref == "lm22") %>%
  filter(pc %in% unique(y$pc)[1:10]) %>%
  ggplot(., aes(x=pc, y=truePcsCors)) +
  geom_violin(fill="lightgray", alpha=0.8) +
  geom_point(aes(fill=celltype), position=position_jitter(width=0.1), size=3, shape=21, color="black", stroke=0.5) +
  stat_summary(fun=median, geom="errorbar", width=0.5,
               aes(ymax=after_stat(y), ymin=after_stat(y)), color="red") +
  coord_cartesian(ylim=c(-1, 1)) +
  scale_color_brewer(palette="Set3") +
  scale_y_continuous(breaks=seq(-1, 1, by=0.1)) +
  labs(y = "", x="", fill="")+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

y %>%
  filter(pc %in% unique(y$pc)[1:10]) %>%
  mutate(truePcsCors = abs(truePcsCors)) %>%
  ggplot(., aes(x=pc, y=truePcsCors)) +
  geom_violin(fill="lightgray", alpha=0.8) +
  stat_summary(fun=median, geom="errorbar", width=0.5,
               aes(ymax=after_stat(y), ymin=after_stat(y)), color="red") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(values = palette) +
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  labs(y = "", x="")+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))


y %>%
  filter(pc %in% unique(y$pc)[1:10]) %>%
  mutate(truePcsCors = abs(truePcsCors)) %>%
  group_by(ref, val, celltype) %>%
  summarise(truePcsCors = max(truePcsCors)) %>%
  group_by(ref, val) %>%
  summarise(truePcsCors = mean(truePcsCors)) %>%
  ggplot(., aes(x=ref, y=truePcsCors)) +
  geom_violin(fill="lightgray", alpha=0.8) +
  stat_summary(fun=median, geom="errorbar", width=0.5,
               aes(ymax=after_stat(y), ymin=after_stat(y)), color="red") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_color_manual(values = palette) +
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  labs(y = "", x="")+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

y %>%
  filter(pc %in% unique(y$pc)[1:20]) %>%
  mutate(truePcsCors = abs(truePcsCors)) %>%
  ggplot(., aes(x=pc, y=truePcsCors)) +
  geom_boxplot()

y %>%
  mutate(truePcsCors = abs(truePcsCors)) %>%
  group_by(ref, val, celltype) %>%
  summarise(truePcsCors = max(truePcsCors)) %>%
  group_by(ref, val) %>%
  summarise(truePcsCors = mean(truePcsCors)) %>%
  ggplot(., aes(x=ref, y=truePcsCors)) +
  geom_boxplot()

y


y <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.sigs.10mar.rds")

sigsResults <- parallel::mclapply(1:nrow(vals.refs.res), function(i){


  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
  valType <- vals.refs.res[i,]$val_type
  valDataset <- vals.refs.res[i,]$val_dataset
  refName <- vals.refs.res[i,]$ref_name
  valDataType <- vals.refs.res[i,]$val_data_type[[1]]

  sigs <- y[[i]]

  # Return signatures correlations
  if (TRUE) {

    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
    samples <- intersect(colnames(truth_mat), colnames(mix.in))
    truth_mat <- truth_mat[,samples]
    genes2use <-intersect(rownames(ref.in), rownames(mix.in))

    mix <- mix.in[genes2use,samples]

    mix_ranked <- singscore::rankGenes(mix)
    celltypes <- intersect(rownames(truth_mat), labels$label)

    all_ctoi_res <- lapply(celltypes, function(ctoi){

      # All sigs
      signatures_ctoi <- sigs[startsWith(names(sigs), paste0(ctoi, "#"))]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })



      fracs <- truth_mat[ctoi,colnames(mix_ranked)]
      all_sigs_cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })




      tibble(ref = refName[[1]], val = valDataset, celltype = ctoi,
             sig_name = names(all_sigs_cors),
             cor = all_sigs_cors) %>%
      rowwise() %>%
      mutate(genes = list(sigs[[sig_name]])) %>%
        return(.)
    })
    names(all_ctoi_res) <- celltypes

    return(bind_rows(all_ctoi_res))
  }



}, mc.cores = 40)
sigsResults <- bind_rows(sigsResults)


sigsResults %>%
  ungroup() %>%
  filter(celltype == "CD8-positive, alpha-beta T cell") %>%
  group_by(ref, val, celltype) %>%
  mutate(mean_cor = mean(cor)) %>%
  ungroup() %>%
  unnest_longer(genes) %>%
  mutate(delta_cor = cor - mean_cor) %>%
  group_by(ref, val, genes) %>%
  summarise(gene_delta_cor = sum(delta_cor)) %>%
  ungroup() %>%
  mutate(gene_delta_cor = ifelse(gene_delta_cor > 0 , 1, -1)) %>%
  group_by(genes) %>%
  summarise(gene_delta_cor = mean(gene_delta_cor)) %>% View()


yy %>%
  unnest_longer(genes) %>%
  ggplot(., aes(x=reorder(genes, cor, median) , y=cor)) +
  geom_boxplot() +
  geom_hline(yintercept = median(yy$cor), linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))



load("/bigdata/almogangel/xCell2/data/celltype.data.rda")

z <- celltype.data %>%
  filter(xCell2_labels == "eosinophil") %>%
  #mutate(essential_genes = list(c("CSTA", "LYZ"))) %>%
  slice(1) %>%
  pull(essential_genes)

x <- celltype.data %>%
  filter(xCell2_labels != "myeloid cell")

celltype.data <- rbind(x,z)

save(celltype.data, file = "/bigdata/almogangel/xCell2/data/celltype.data.rda")



y <- readRDS("/bigdata/almogangel/xCell2_data/benchmarking_data/results/correlations/xcell2.sigs.noEssential.10mar.rds")
sigsResults <- parallel::mclapply(1:nrow(vals.refs.res), function(i){


  # Load data
  val_ref <- paste0(vals.refs.res[i,]$val_dataset, "_", vals.refs.res[i,]$ref_name[[1]])
  print(val_ref)
  mix.in <- cyto.vals$mixtures[[vals.refs.res[i,]$val_type]][[vals.refs.res[i,]$val_dataset]]
  ref.in <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$ref
  labels <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$labels
  lineage_file <- refsRDSList[[vals.refs.res[i,]$ref_type]][[vals.refs.res[i,]$ref_name[[1]]]]$lineage_file
  refType <- ifelse(vals.refs.res[i,]$ref_type == "rna_seq", "rnaseq", vals.refs.res[i,]$ref_type)
  valType <- vals.refs.res[i,]$val_type
  valDataset <- vals.refs.res[i,]$val_dataset
  refName <- vals.refs.res[i,]$ref_name
  valDataType <- vals.refs.res[i,]$val_data_type[[1]]

  sigs <- y[[i]]

  # Return signatures correlations
  if (TRUE) {

    truth_mat <- cyto.vals$truth[[valType]][[valDataset]]
    samples <- intersect(colnames(truth_mat), colnames(mix.in))
    truth_mat <- truth_mat[,samples]
    genes2use <-intersect(rownames(ref.in), rownames(mix.in))

    mix <- mix.in[genes2use,samples]

    mix_ranked <- singscore::rankGenes(mix)
    celltypes <- intersect(rownames(truth_mat), labels$label)

    all_ctoi_res <- lapply(celltypes, function(ctoi){

      # All sigs
      signatures_ctoi <- sigs[startsWith(names(sigs), paste0(ctoi, "#"))]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        suppressWarnings(singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore)
      })



      fracs <- truth_mat[ctoi,colnames(mix_ranked)]
      all_sigs_cors <- apply(scores, 2, function(x){
        cor(x, fracs, method = "spearman", use = "pairwise.complete.obs")
      })




      tibble(ref = refName[[1]], val = valDataset, celltype = ctoi,
             sig_name = names(all_sigs_cors),
             cor = all_sigs_cors) %>%
        rowwise() %>%
        mutate(genes = list(sigs[[sig_name]])) %>%
        return(.)
    })
    names(all_ctoi_res) <- celltypes

    return(bind_rows(all_ctoi_res))
  }



}, mc.cores = 40)
sigsResults <- bind_rows(sigsResults)

data.in %>%
  filter(is_filtered == "no")
