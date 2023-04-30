
# allRes <- readRDS("/bigdata/almogangel/twelve_years_decon_paper/analysis/results/accuracy/allRes.rds")


allRes <- allRes[,c(1,2,4)]
cbrx.df <- data.frame("dataset" = names(cbrx_res), "method" = "CIBERSORTx", "CCorr" = cbrx_res)
allRes <- rbind(allRes, cbrx.df)
allRes <- allRes[rownames(allRes) != "Liver",]

df <- pivot_wider(allRes, names_from = method, values_from = CCorr)
df <- data.frame(df[,-1], row.names = df$dataset, check.names = FALSE)
#df[df < 0] <- NA
methods_sorted <- names(sort(apply(df, 2, function(x){median(x, na.rm = T)}), decreasing = F))

allRes$method <- factor(allRes$method, levels = methods_sorted)



allRes %>%
  as_tibble() %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2"), "yes", "no")) %>%
  ggplot(., aes(x=method, y=CCorr)) +
  geom_boxplot(aes(fill=is_xcell2), position = position_dodge(1), outlier.shape = NA,  alpha = 0.5) +
  geom_jitter(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1), labels = as.character(seq(-1,1,0.1))) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96", "#8B4513", "#2F4F4F", "#54FF9F")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  coord_flip() +
  theme_linedraw() +
  theme(plot.title = element_text(size=22, hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(colour = "#1A1A1A", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "white"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, angle = 10, hjust=1, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  labs(y = "Spearman r", title = "Tabula Sapiens Validation", x = NULL, colour = NULL, fill = NULL)
