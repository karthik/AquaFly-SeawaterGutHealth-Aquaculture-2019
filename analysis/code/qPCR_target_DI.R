library(tidyverse) # tidy, manipulate and plot data
library(cowplot) #  a ggplot2 add-on to provide a publication-ready theme  
library(nlme)  # run linear mixed model
library(ComplexHeatmap) # make heatmap

# Import and tidy data ################################################################################################
# Set system locale to "Greek" so that greek letters in the gene names can be recognized
Sys.setlocale(category = "LC_ALL", locale = "Greek")

df <- read.csv("./data/clean_data/AqFl2_qPCR_target_DI.csv", 
               header = T, 
               na.strings = c(""),
               encoding = "UTF-8", 
               stringsAsFactors = FALSE)

# The encoding of the csv file need to be UTF-8, otherwise greek letters won't be recognized. To do that:
# 1. Save the excel sheet as a Unicode text.
# 2. Open your saved file in Notepad and replace all tab characters with commas (","). 
#   2.1 Select a tab character (select and copy the space between two column headers.
#   2.2 Open the "Find and Replace" window and replace all tab characters with comma. 
# 3. Click Save As, name the file, and change the Encoding: to UTF-8.
# 4. Change the file extension from "*.txt" to "*.csv".
# NB: You can view the UTF-8 file in excel, but don't save it when close it. Otherwise, you need to do it again.

head(df, n = 20L)
str(df)

# Tidy data and calculate new variables 
df <- df %>% 
  rename(Sample_ID = X.U.FEFF.Sample_ID) %>% # remove "X.U.FEFF.", which was reserved to tell that the file is in UTF-8 format
  select(-Comments) %>% # remove the "Comments" column as it contains many NA values
  filter(Target != "il6") %>% # remove il6 from the analysis as the expression level was too low 
  na.omit() %>% # remove rows with NA value
  mutate(NE_gapdh = `^`(Target_PE,  -Target_Cq) / `^`(gapdh_PE, -gapdh_Cq), # normalzie expression to gapdh
         NE_ranpo2 = `^`(Target_PE,  -Target_Cq) / `^`(rnapo2_PE, -rnapo2_Cq), # normalzie expression to rnapo2
         NE_hprt1 = `^`(Target_PE,  -Target_Cq) / `^`(hprt1_PE, -hprt1_Cq), # normalzie expression to hprt1
         NE_actb = `^`(Target_PE,  -Target_Cq) / `^`(actb_PE, -actb_Cq), # enormalzie expression to actb
         MNE = (NE_gapdh * NE_ranpo2 * NE_hprt1 * NE_actb)^(1/4)) # calculate geometric mean of normalzied expressions

head(df, n = 20L) # For sanity, check the dataframe again

# Exploratory analysis ################################################################################################
# Convert numeric variable to character
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Convert character to factor, specifying the desired orders to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Target <- factor(df$Target, levels = c("myd88", "tnfα", "il1β", "il8", "cd3γδ", "cd8β", "mhcI", "ifnγ", 
                                          "il17a", "foxp3", "il10", "tgfβ1", "il4", "plin2", "mtp", "apoa1",  
                                          "apoa4", "apob", "chk", "pcyt1a", "fabp2b", "aqp8ab", "cldn15",  
                                          "cldn25b", "zo1", "cdh1", "muc2",  "cyp1a1", "mta", "hsp70", "sod1", 
                                          "cat", "casp6", "pcna", "mmp13")) 

# Violin plot to check data distribution ----------------------------------------------------------
df %>% 
  ggplot(aes(x = Diet, y = MNE)) + 
  geom_violin(aes(fill = Diet), trim = FALSE) +
  facet_wrap(~ Target, nrow = 9, scales = "free_y") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange") + # add mean and SD 
  scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Gene expression profile in DI", y = "MNE") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('./results/exploratory_analysis/qPCR_target_DI_violin.pdf', 
       units = "in", 
       dpi = 300, 
       encoding = "CP1253", # print greek letters in the correct format
       height = 12, 
       width = 8)

# Box plot to check outliers and explore for cluster effect, i.e., net pen ------------------------
# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

df %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group for each gene
  ggplot(aes(x = Diet, y = MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen), size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~Target, nrow = 9, scales = "free_y") +
  ggrepel::geom_label_repel() +  # label outliers with Sample_ID
  scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Gene expression profile in DI", y = "MNE") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('./results/exploratory_analysis/qPCR_target_DI_boxPlot.pdf', 
       units = "in", 
       dpi = 300, 
       encoding = "CP1253", 
       height = 24, 
       width = 16)

# Some samples were identified as outliers in several genes. 
# Check sampling records and raw data for any warning notes or mistakes

# Heatmap -----------------------------------------------------------------------------------------
# Convert dataframe from "long" to "wide" format
df_w <- reshape2::dcast(df, Sample_ID + Diet + Net_pen ~ Target, value.var = "MNE") %>%
  arrange(Diet, Net_pen) # arrange samples by diet and net_pen makes it easier to annotate heatmap

# Make a numeric matrix for plotting heatmap
nm  <- select(df_w, -Sample_ID, -Diet, -Net_pen) %>% 
  data.matrix() %>% 
  scale() %>% # due to the huge differences in the normalized expression levels, data need to be scaled
  t() # transpose the matrix so that samples can be shown on columns and genes on rows

colnames(nm) <- df_w[ , 1] # assign Sample_ID as column names

# Make column annoation for heatmap
annCol = HeatmapAnnotation(df = df_w[ , 2:3], 
                           col = list(Diet = c("REF" = "red", "IM" = "blue"), # customize color
                                      Net_pen = c("101" = "#1B9E77", "106" = "#D95F02","111" = "#7570B3",
                                                  "104" = "#E7298A","105" = "#66A61E","110" = "#E6AB02")
                           ), # The color panel was generated by RColorBrewer::brewer.pal(n = 6, name = "Dark2")
                           show_annotation_name = TRUE,
                           annotation_name_offset = unit(2, "mm"))

# Make row annoation for heatmap
x <- c("Immune modulation", "lipid metabolism", "barrier function", "xenobiotic metabolism")
times <- c(13, 8, 6, 8)
annRow <- data.frame(func = rep(x, times = times))
annRow$func <- factor(annRow$func, levels = c("Immune modulation",
                                              "lipid metabolism",
                                              "barrier function",
                                              "xenobiotic metabolism"))

pdf('./results/exploratory_analysis/qPCR_target_DI__heatmap.pdf', 
    encoding = "CP1253", 
    height = 10,
    width = 8.5)

Heatmap(nm, 
        name = "SNE", 
        split = annRow$func, # split heatmap rows and perform clustering with each functional category
        gap = unit(1, "mm"), 
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        clustering_distance_rows = "spearman", # cluster genes based on spearman correlation.
        column_names_gp = gpar(fontsize = 10),
        column_title = "Sample ID",
        column_title_side = "bottom",
        row_names_gp = gpar(fontsize = 12, fontface = "italic"),
        row_names_side = "right",
        row_dend_side = "right",
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 12),
        row_dend_width = unit(2.5, "cm"),
        rect_gp = gpar(col = "gray50", lty = 1, lwd = 0.2), # color, line type and size of cell boarder
        top_annotation = annCol, # column annotation
        col = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

dev.off()

# Statistics ##########################################################################################################
# Split the dataframe by "Target"
df_spl <- split(df, f = df$Target)

# Run linear mixed model for each data frame using lapply()
lme <- lapply(df_spl, function(x) lme(MNE ~ Diet, 
                                      random = ~ 1|Net_pen, 
                                      varPower(form = ~ fitted(.)), 
                                      data = x))

# Get model summary and extract p values of the Diet effect.
smr <- lapply(lme, function(x) summary(x))
anv <- lapply(lme, function(x) anova(x))

p_val <- data.frame(row.names = NULL,
                    "Gene_of_interest" = names(lme),
                    "p_raw" = unlist(lapply(anv, function(x) x[[4]][[2]]))
)

p_val$p_adj <- p.adjust(p_val$p_raw, method = "fdr")

# Format digits of p values
p_val$p_raw  <- formatC(p_val$p_raw, format = "f", digits = 3)
p_val$p_adj  <- formatC(p_val$p_adj, format = "f", digits = 3)

# Model diagnostics -------------------------------------------------------------------------------
# Extract residuals
res <- lapply(lme, function(x) resid(x))

# 1.Linearity and homoskedasticity of residuals
pdf(file = "./results/statistics/qPCR_target_DI_resid_plot.pdf", encoding = "CP1253") 

lapply(
  seq_along(res), # add index to the elements in the list "res" so that element names and contents can be extracted
  function(x) 
  {
    plot(res[[x]], main = names(res)[x], ylab = "residuals") # make residual plot for each model
    abline(0,0)
  }
)

dev.off() 

# 2.Absence of collinearity. When more than one fixed effects are included, the collinearity shoulb be checked.
# Not applicable to the present dataset, which has only one fixed effect

# 3.Normality of residuals

pdf(file = "./results/statistics/qPCR_target_DI_qqplot.pdf", encoding = "CP1253")

lapply(
  seq_along(res), 
  function(x) 
  {
    qqnorm(res[[x]], main = paste("Normal Q-Q Plot: ", names(res)[x])) # make qq plot for each model
    qqline(res[[x]])
  }
)

dev.off()

# 4.Absence of influential data points

# Visual inspection of residual plot and qq plot revealed obvious deviations from homoscedasticity and normality for 
# some genes, including ... 
# The linear mixed effect model need to be run by modelling the variance using the varPower function.

# Make Figure 4 #######################################################################################################
# add raw and adjusted p values to row annotation
annRow1 = rowAnnotation(p_raw = row_anno_text(p_val$p_raw, gp = gpar(fontsize = 10)),
                        p_adj = row_anno_text(p_val$p_adj, gp = gpar(fontsize = 10)),
                        show_annotation_name = TRUE,
                        annotation_name_gp = gpar(fontsize = 10),
                        annotation_name_rot = c(0),
                        annotation_width = unit(c(1, 1), "cm"))   

tiff('./results/figures/Figure 4.tiff', 
     compression = "lzw", 
     units = "in", 
     res = 300,
     height = 10,
     width = 8.5)

annRow1 + Heatmap(nm, # the order of row annotation and heatmap matters determines how they're arranged
                  name = "SNE", 
                  split = annRow$func, 
                  gap = unit(1, "mm"), 
                  cluster_columns = FALSE,
                  cluster_rows = TRUE,
                  clustering_distance_rows = "spearman", 
                  column_names_gp = gpar(fontsize = 10),
                  column_title = "Sample ID",
                  column_title_side = "bottom",
                  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
                  row_names_side = "right",
                  row_dend_side = "right",
                  row_title_side = "right",
                  row_title_gp = gpar(fontsize = 12),
                  row_dend_width = unit(2.5, "cm"),
                  rect_gp = gpar(col = "gray50", lty = 1, lwd = 0.2), 
                  top_annotation = annCol, 
                  col = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

dev.off()

# Get session info 
writeLines(capture.output(sessionInfo()), "./code/qPCR_target_DI_sessionInfo.txt")
