library(tidyverse) # tidy, manipulate and plot data
library(ggpubr) # ggplot2-based publication-ready plots
library(ComplexHeatmap) # make heatmap
library(lme4) # linear mixed effect modelling (LMMs)
library(influence.ME) # detect influential observations in the LMMs
library(pbkrtest) # parametric bootstrap for testing fixed effects in LMMs

# Import and tidy data #########################################################
# Set system locale to "Greek" so that greek letters in the gene names can be recognized
Sys.setlocale(category = "LC_ALL", locale = "Greek")

# Import data 
df <- read_csv("data/clean_data/AqFl2_qPCR_Target_DI.csv", 
               col_names = T, na = "") 

# The encoding of the csv file need to be UTF-8, otherwise greek letters won't be recognized. To do that:
# 1. Save the excel sheet as a Unicode text.
# 2. Open your saved file in Notepad and replace all tab characters with commas (","). 
# 3. Click Save As, name the file and change the Encoding to UTF-8.
# 4. Change the file extension from ".txt" to ".csv".
# NB: You can view the UTF-8 file in excel, but don't save it when close it. Otherwise, you have to do it again.

head(df, n = 20L)
str(df)

# Tidy data and calculate new variables 
df <- df %>% 
  select(-Comments) %>% # remove the "Comments" column as it contains many NA values
  filter(Target != "il6") %>% # remove il6 from the analysis as the expression level was too low 
  na.omit() %>% # remove rows with NA value
  mutate(NE_gapdh  = `^`(Target_PE,  -Target_Cq) / `^`(gapdh_PE, -gapdh_Cq), # normalzie expression to gapdh
         NE_ranpo2 = `^`(Target_PE,  -Target_Cq) / `^`(rnapo2_PE, -rnapo2_Cq), # normalzie expression to rnapo2
         NE_hprt1  = `^`(Target_PE,  -Target_Cq) / `^`(hprt1_PE, -hprt1_Cq), # normalzie expression to hprt1
         NE_actb   = `^`(Target_PE,  -Target_Cq) / `^`(actb_PE, -actb_Cq), # enormalzie expression to actb
         MNE       = (NE_gapdh * NE_ranpo2 * NE_hprt1 * NE_actb)^(1/4)) # calculate geometric mean of normalzied expressions

head(df, n = 20L) 

# Exploratory analysis #########################################################
# Convert numeric variable to character
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Convert character to factor, specifying the desired orders to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Target <- factor(df$Target, 
                    levels = c("myd88", "tnfα", "il1β", "il8", "cd3γδ", "cd8β", 
                               "mhcI", "ifnγ","il17a", "foxp3", "il10", "tgfβ1", 
                               "il4", "plin2", "mtp", "apoa1", "apoa4", "apob", 
                               "chk", "pcyt1a", "fabp2b", "aqp8ab", "cldn15", 
                               "cldn25b", "zo1", "cdh1", "muc2",  "cyp1a1", "mta", 
                               "hsp70", "sod1", "cat", "casp6", "pcna", "mmp13")) 

# Violin plot to check data distribution ---------------------------------------
df %>% 
  ggplot(aes(x = Diet, y = MNE)) + 
  geom_violin(aes(fill = Diet), trim = T) +
  facet_wrap(~ Target, nrow = 9, scales = "free_y") +
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1), 
               geom     = "pointrange") + # add mean and SD 
  scale_y_continuous(limits = c(0, NA), 
                     expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Gene expression profile in DI", y = "MNE") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('analysis/exploratory_analysis/qPCR_target_DI_violin.pdf', 
       units = "in", 
       dpi = 300, 
       encoding = "CP1253", # print greek letters in the correct format
       height = 12, 
       width = 8)

# Boxplot to check outliers ----------------------------------------------------
# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Set seed to make jitter points reproducible
set.seed(1910)

# Make boxplots
df %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group for each gene
  ggplot(aes(x = Diet, y = MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen), size = 2, shape = 21, 
             position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, nrow = 9, scales = "free_y") +
  ggrepel::geom_label_repel() +  # label outliers with Sample_ID
  scale_y_continuous(limits = c(0, NA), 
                     expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Gene expression profile in DI", y = "MNE") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('analysis/exploratory_analysis/qPCR_target_DI_boxPlot.pdf', 
       units = "in", 
       dpi = 300, 
       encoding = "CP1253", 
       height = 24, 
       width = 16)

# Heatmap ----------------------------------------------------------------------
# Convert dataframe from "long" to "wide" format
df_w <- reshape2::dcast(df, Sample_ID + Diet + Net_pen ~ Target, 
                        value.var = "MNE") %>%
  arrange(Diet, Net_pen) # arrange samples by diet and net_pen makes it easier to annotate heatmap

# Make a numeric matrix for plotting heatmap
mat <- select(df_w, -Sample_ID, -Diet, -Net_pen) %>% 
  data.matrix() %>% 
  scale() %>% # due to the huge differences in the normalized expression levels, data need to be scaled
  t() %>% # transpose the matrix so that samples can be shown on columns and genes on rows
  `colnames<-`(df_w[ , 1]) # assign Sample_ID as column names

# Make column annoations for heatmap
annCol <- HeatmapAnnotation(df = df_w[ , 2:3], 
                            col = list(Diet = c("REF" = "red", "IM" = "blue"), # customize color
                                       Net_pen = c("101" = "#1B9E77", 
                                                   "106" = "#D95F02",
                                                   "111" = "#7570B3", 
                                                   "104" = "#E7298A",
                                                   "105" = "#66A61E",
                                                   "110" = "#E6AB02")), # RColorBrewer::brewer.pal(n = 6, name = "Dark2")
                            show_annotation_name = T,
                            annotation_name_offset = unit(2, "mm"))
# Make heatmap
pdf('analysis/exploratory_analysis/qPCR_target_DI_heatmap.pdf', 
    encoding = "CP1253", 
    height = 10,
    width = 8.5)

Heatmap(mat, 
        name = "SNE", 
        cluster_rows = T,
        clustering_distance_rows = "spearman", # cluster genes based on spearman correlation
        clustering_method_rows = "ward.D2",
        row_dend_side = "right",
        row_dend_width = unit(2.5, "cm"),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 12, fontface = "italic"),
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 12),
        cluster_columns = T,
        clustering_distance_columns = "euclidean", # cluster samples based on euclidean distance
        clustering_method_columns = "ward.D2",
        column_dend_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_title = "Sample ID",
        column_title_side = "bottom",
        show_parent_dend_line = F,
        rect_gp = gpar(col = "gray50", lty = 1, lwd = 0.2), # color, line type and size of cell boarder
        top_annotation = annCol, # column annotation
        col = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

dev.off()

# Statistics ###################################################################
# Split the dataframe by target genes
df_spl <- split(df, f = df$Target)

# As the number of random-effect levels is small (6), some fitted models may be 
# singular. Hence, we'll first fit models using maximal likelihood estimation 
# and filter variables showing singular fits

# Based on the exploratory anaylysis, we'll also log transform the response 
# variables to stabilize the mean-variance relationship

lmerML <- lapply(df_spl, function(x) lmer(log10(MNE) ~ Diet + (1|Net_pen), 
                                          REML = F, data = x)) 

# 21 fitted models are singular. We'll not fit LMMs to these variables, instead 
# we'll run Welch' t-test to compare group means

# Welch's t-test ---------------------------------------------------------------
# Find out target genes with singular fits
target_sf <- map_df(lmerML, ~VarCorr(.x)$Net_pen[1]) %>% 
  gather(key = Target) %>%
  filter(value == 0)

# Make a list of dataframes containing target genes with singular fits
df_sf <- filter(df, Target %in% target_sf$Target) %>%
  mutate(Target = as.character(Target), # Target was previously converted to facotr with 34 levels
         MNE = log10(MNE)) %>% # log transform to stabalize mean-variance relationship
  split(f = .$Target)

# Welch's t-test
welch_t <- lapply(df_sf, function(x) t.test(MNE ~ Diet, x)) 

# Check the normality assumption visually via Q-Q plots
qq_welchT <- lapply(
  seq_along(df_sf), 
  function(x) 
  {
    ggqqplot(df_sf[[x]], "MNE", color = "Diet", 
             facet.by = "Diet", title = names(df_sf)[x]) 
    
  }
)

ggexport(filename = "analysis/model_diagnostics/qPCR_DI_welch-t_qqplot.pdf",
         plotlist = qq_welchT, nrow = 2, ncol = 2, width = 8)

# Fit LMMs with REML -----------------------------------------------------------
# Make a list of dataframes containing target genes not showing singular fits
df_nsf <- filter(df, !Target %in% target_sf$Target) %>%
  mutate(Target = as.character(Target), 
         MNE = log10(MNE)) %>% 
  split(f = .$Target)

# Fit LMMs
lmerREML <- lapply(df_nsf, 
                   function(x) 
                   lmer(MNE ~ Diet + (1|Net_pen), REML = T, data = x)) 

### Model diagnostics ###

## 1.Residual analysis ##
# Extract residuals
res <- lapply(lmerREML, function(x) resid(x))

# Check homogeneity and normality of residuals
pdf(file = "analysis/model_diagnostics/qPCR_DI_LMMs_residual.pdf", 
    encoding = "CP1253",
    width = 8) 

mapply(function(x, y){
  par(mfrow = c(2, 2), mar = c(5, 5, 2, 2), cex.lab = 1.5)
  
  plot(res[[x]], xlab = "Fitted values", ylab = "Residuals") # plot residuals against fitted values
  abline(0,0)
  
  boxplot(res[[x]] ~ Diet, data = y, xlab = "Diet", ylab = "Residuals") # plot residuals against Diet
  abline(h = 0, lty = 2)
  
  boxplot(res[[x]] ~ Net_pen, data = y, xlab = "Net pen", ylab = "Residuals") # plot residuals against Net pen
  abline(h = 0, lty = 2)
  
  qqnorm(res[[x]], main = "Normal Q-Q Plot") # make a Q-Q plot 
  qqline(res[[x]])
  
  par(oma = c(0,0,2,0))
  title(main = names(res)[x], outer = T)}, 
  
  x = seq_along(res), 
  y = df_nsf
)

dev.off() 

## 2.Influence analysis ##

# Based on the residual analysis, we'll run influence analysis for il8, apob, 
# mtp and fabp2b.

# Due to a bug in the influence.ME package, we can't apply influence() function 
# to elements in the list "lmerREML" directly. We'll iteratively make a dataframe
# subset, fit LMMs and run the influence analyses.

for (var in c("il8", "apob", "mtp", "fabp2b")) { 
  # Subset dataframe and fit LMMs again
  df_sub <- filter(df, Target == var) 
  lmer_sub <- lmer(MNE ~ Diet + (1|Net_pen), REML = T, data = df_sub) 
  
  # Determine the extent to which the parameter estimates in the model change, 
  # when each of the observations is deleted from the data iteratively 
  estex <- influence(lmer_sub, obs = T)
  
  # Detect influence observations via Cook's distance
  plot(estex, which = "cook", cutoff = 4/nrow(df_sub), sort = T, 
       main = var, xlab = "Cook's Distance", ylab = "Row Number")
  
  # Detect influence observations via sigtest
  print(paste0("sigtest for ", var))
  
  if(fixef(lmer_sub)[2] > 0){
    print(sigtest(estex, test = 1.96)$Diet)
  } else if (fixef(lmer_sub)[2] < 0) {
    print(sigtest(estex, test = -1.96)$Diet)
  } else {
    print("NB: the paramter of the fixed effect was estimated to be zero!")
  }
} 

# Get p values of the Diet effect ----------------------------------------------
## Get p values from welch' t test ##
p_welch <- map_df(welch_t, ~.x$p.value) %>%
  gather(key = Gene, value = p_raw)

## Get p values from LMMs using parametric bootstrap ##
# Remove singular fits from LMMs estimated by maximal likelihood
lmerML_nsf <- lmerML[c(names(df_nsf))]

# Construct reduced models
lmerML_nsf_nofix <- lapply(lmerML_nsf, function(x) update(x, . ~ . - Diet))

# Create clusters for parallel computing to speed up parametric bootstrap 
nc <- parallel::detectCores() 
clus <- parallel::makeCluster(rep("localhost", nc)) 

# Set seed to make parametric bootstrap reproducible
set.seed(1910)

# Parametric bootstrap comparisons 
pb <- map2(lmerML_nsf, lmerML_nsf_nofix, ~PBmodcomp(.x, .y, cl = clus))

# Stop the clusters
parallel::stopCluster(clus)

# Get the p values of parametric bootstrap comparisons
p_pb <- map_df(pb, ~.x$test["PBtest", "p.value"]) %>%
  gather(key = Gene, value = p_raw)

# Gather p values 
p_val <- bind_rows(p_welch, p_pb) %>%
  mutate(p_adj = p.adjust(p_raw, method = "fdr"), # FDR based p value correction
         p_raw = formatC(p_raw, format = "f", digits = 3), # format digits
         p_adj = formatC(p_adj, format = "f", digits = 3)) %>% 
  arrange(match(Gene, row.names(mat))) # match the order of genes in the matrix used for plotting heatmap

# Make Figure 4 ################################################################
# Make a dataframe for splitting heatmap rows by gene functions
splitRow <- select(df, Target, Target_func) %>%
  unique() %>%
  mutate(Target_func = gsub("\\_", " ", Target_func),
         Target_func = factor(Target_func, 
                              levels = c("immune modulation",
                                         "lipid metabolism",
                                         "barrier function", 
                                         "xenobiotic metabolism")))

# Add raw and adjusted p values to row annotations
annRow <- rowAnnotation(p_raw = anno_text(p_val$p_raw, gp = gpar(fontsize = 10)),
                        p_adj = anno_text(p_val$p_adj, gp = gpar(fontsize = 10)),
                        annotation_width = unit(c(1, 1), "cm"))   

tiff('results/figures/Figure 4.tiff', 
     compression = "lzw", 
     units = "in", 
     res = 300,
     height = 10,
     width = 8.5)

Heatmap(mat, 
        name = "SNE", 
        row_split = splitRow$Target_func, # split heatmap rows and perform clustering with each functional category
        cluster_row_slices = F,
        row_gap = unit(1, "mm"), 
        cluster_columns = T,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        column_dend_side = "top",
        column_split = df_w$Diet, 
        cluster_column_slices = F,
        column_gap = unit(1, "mm"),
        cluster_rows = T,
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
        show_parent_dend_line = F,
        left_annotation = annRow, # row annotation
        top_annotation = annCol, # column annotation
        rect_gp = gpar(col = "gray50", lty = 1, lwd = 0.2), # color, line type and size of cell boarder
        col = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

# The show_annotation_name() for anno_text() is turned off in this version of 
# ComplexHeatmap. we'll add annotation names "p_raw" and "p_adj" via 
# decorate_annotation()
decorate_annotation("p_raw", slice = 4, {
  grid.text("p_raw", y = unit(-2, "mm"), hjust = "left", 
            gp = gpar(fontsize = 10, fontface = "bold"))
})

decorate_annotation("p_adj", slice = 4, {
  grid.text("p_adj", y = unit(-2, "mm"), hjust = "left", 
            gp = gpar(fontsize = 10, fontface = "bold"))
})

dev.off()

# Get session info 
writeLines(capture.output(sessionInfo()), "analysis/code/qPCR_target_DI_sessionInfo.txt")