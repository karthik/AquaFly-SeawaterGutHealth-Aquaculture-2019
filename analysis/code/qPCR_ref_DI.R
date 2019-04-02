library(tidyverse) # tidy, manipulate and plot data
library(broom) # tidy statistical outputs
library(cowplot) #  a ggplot2 add-on to provide a publication-ready theme  

# Import and tidy data ################################################################################################
df <- read.csv("./data/clean_data/AqFl2_qPCR_ref_DI.csv", header = T, na.strings = c(""), stringsAsFactors = FALSE)
head(df, n = 20L)
str(df)

# Tidy data and calculate new variables
df <- df %>% 
  arrange(Ref_gene, desc(Diet)) %>% 
  mutate(Quantity = `^`(PE, -Cq)) %>% # calculate mRNA quantity. PE: primer efficiency
  group_by(Ref_gene, Diet)  %>%
  mutate(Mean = mean(Quantity), SD = sd(Quantity)) %>% # calculate average mRNA quantity per diet and its sd
  ungroup() %>%
  group_by(Ref_gene) %>%
  mutate(Quantity_std = scale(Quantity), # standardize mRNA quantity of each reference gene
         Quantity_nm = Quantity / Quantity[1]) # normalize mRNA quantity of each reference gene to the first sample

# Exploratory analysis ################################################################################################ 
# Convert character to factor, specifying the desired orders to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Sample_ID <- factor(df$Sample_ID, levels = unique(df$Sample_ID)) # use the exact order of Sample_ID column

# Point diagram -----------------------------------------------------------------------------------
# make a dataframe containing the grand mean of mRNA quantity per reference gene
gm <- df %>% 
  group_by(Ref_gene) %>%
  summarize(gmean = mean(Quantity))

# make point diagram 
# Show individual expression profile
p1 <- df %>% 
  ggplot(aes(x = Sample_ID, y = Quantity, group = 1)) +
  geom_point(aes(color = Diet)) +
  geom_line() +
  facet_wrap(~ Ref_gene, ncol = 1, scales = "free_y") +
  labs(y = "mRNA quantity", tag = "A") +
  scale_y_continuous(labels = scales::scientific) + 
  theme_bw() +
  geom_hline(aes(yintercept = gmean), gm, color = "blue") # The blue line shows the grand mean.

# Show normalized expression levels of the 4 reference genes in one plot
p2 <- df %>% 
  ggplot(aes(x = Sample_ID, y = Quantity_nm, color = Ref_gene)) +
  geom_point() +
  geom_line(aes(group = Ref_gene)) +
  labs(y = "normalized mRNA quantity", tag = "B") +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") 

# Combine the plots
plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(3, 2))

ggsave('./analysis/exploratory_analysis/qPCR_ref_DI_pointDiagram.pdf', 
       units = "in", 
       dpi = 300,
       height = 10,
       width = 8)

# Boxplot -----------------------------------------------------------------------------------------
# convert Sample_ID as character, otherwise ggrepel won't work properly
df$Sample_ID <- as.character(df$Sample_ID)

# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Make boxplot using standardized mRNA quantity
df %>% 
  group_by(Ref_gene) %>%
  mutate(outlier = is_outlier(Quantity_std)) %>%
  ggplot(aes(x = Ref_gene, y = Quantity_std, label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Diet), size = 3, shape = 21, position = position_jitterdodge(0.2)) +
  ggrepel::geom_label_repel() + # label outliers with Sample_ID
  ylab("Standardized mRNA quantity") +
  theme_cowplot()

ggsave('./analysis/exploratory_analysis/qPCR_ref_DI_boxPlot.pdf', units = "in", dpi = 300)

# Barplot -----------------------------------------------------------------------------------------
df %>%
  ggplot(aes(x = Diet, y =  Mean, fill = Diet)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), size = 0.3, width = 0.2) + # add error bar (sd)
  facet_wrap(~ Ref_gene, nrow = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) + 
  ylab("mRNA quantity") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave('./analysis/exploratory_analysis/qPCR_ref_DI_barPlot.pdf', units = "in", dpi = 300)

# Rank reference genes by coefficient of variation and F statistic ###################################################
# Calculate coefficient of variance (cv), measuring overall stability of reference genes across all the samples
cv <- df %>% 
  group_by(Ref_gene) %>%
  summarize(Cq_max = max(Cq), 
            Cq_min = min(Cq), 
            Cq_range = max(Cq) - min(Cq), 
            Quantity_mean = mean(Quantity), 
            Quantity_SD = sd(Quantity),
            Quantity_CV = 100 * Quantity_SD / Quantity_mean) %>%
  mutate(CV_rank = dense_rank(Quantity_CV)) %>% # rank reference genes by CV
  arrange(Ref_gene)

# Calculate F-statistic, measuring reference gene stability across experimental treatments
F_stat <- df %>% 
  group_by(Ref_gene) %>% 
  do(tidy(aov(Quantity ~ Diet, data = .))) %>% # run ANOVA for each reference gene and tidy statistical outputs
  na.omit() %>%
  select(Ref_gene, statistic) %>%
  rename(Quantity_F = statistic) %>%
  as.data.frame() %>%
  mutate(F_rank = dense_rank(Quantity_F)) %>% # rank reference genes by F statistic
  arrange(Ref_gene) 

# Summary 
smr <- full_join(cv, F_stat, by = "Ref_gene")

# Format digits
smr$Quantity_mean <- formatC(smr$Quantity_mean, format = "e", digits = 1)
smr$Quantity_SD <- formatC(smr$Quantity_SD, format = "e", digits = 2)
smr$Quantity_CV <- formatC(smr$Quantity_CV, format = "f", digits = 1)
smr$Quantity_F <- formatC(smr$Quantity_F, format = "f", digits = 4)

# Export summary table
pdf('./results/reference_gene_ranks/ref_gene_rank_DI.pdf', width = 12)
gridExtra::grid.table(smr)
dev.off()

# Get session info
writeLines(capture.output(sessionInfo()), "./analysis/code/qPCR_ref_DI_sessionInfo.txt")
