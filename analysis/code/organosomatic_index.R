library(tidyverse) # tidy, manipulate and plot data
library(ggpubr) # make ggplot2 based publication-ready plots
library(ggsignif) # add p values to plots
library(lme4) # linear mixed effect modelling (LMMs)
library(influence.ME) # detect influential observations in the LMMs
library(pbkrtest) # parametric bootstrap for testing fixed effects in LMMs

# Import and tidy data #########################################################
df <- read_csv("data/raw_data/AqFl2_organosomatic_index.csv", 
               col_names = T, 
               na = "")

head(df)
str(df)

# Calculate organosomatic index (OSI)
df <- df %>%
  mutate(OSI = 100 * Organ_weight / Body_weight)

# Exploratory analysis ######################################################### 
# Define a function for setting the number of decimals shown on axis
fmt_dcimals <- function(decimals = 0){
  function(x) format(x, nsmall = decimals, scientific = F)
}

# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Convert numeric variables to character
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Convert character to factor, specifying the desired orders to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))

# Violin plot to check data distribution ---------------------------------------
ggplot(df, aes(x = Diet, y = OSI)) + 
  geom_violin(aes(fill = Diet), trim = F) +
  facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1), 
               geom = "pointrange") + # add mean and SD 
  scale_y_continuous(limits = c(0, NA), 
                     labels = fmt_dcimals(1), 
                     expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Organosomatic indices", y = "%") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('analysis/exploratory_analysis/OSI_violin.pdf', units = "in", dpi = 300)

# Box plot to check outliers and explore random effect, i.e., net pen ----------
df %>%
  group_by(Gut_segment, Diet) %>%
  mutate(outlier = is_outlier(OSI)) %>% 
  ggplot(aes(x = Diet, y = OSI, label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(0.2)) +
  facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
  ggrepel::geom_label_repel() + # label outliers using Sample_ID
  scale_y_continuous(limits = c(0, NA), 
                     labels = fmt_dcimals(1), 
                     expand = expand_scale(mult = c(0, 0.1))) +
  labs(title = "Organosomatic indices", y = "%") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('analysis/exploratory_analysis/OSI_boxPlot.pdf', units = "in", dpi = 300)

# Statistics ###################################################################
# Split the data frame by Gut_segment
df_spl <- split(df, f = df$Gut_segment)

# As the number of random-effect levels is small (6), some fitted models may be 
# singular. Hence, we'll first fit models using maximal likelihood estimation and 
# filter variables producing singular fits

lmerML <- lapply(df_spl, function(x) lmer(OSI ~ Diet + (1|Net_pen), 
                                          REML = F, data = x)) 

# 2 fitted models are singular. We'll not fit LMMs to these variables, instead 
# we'll run Welch' t-test to compare group means

# Welch's t-test ---------------------------------------------------------------
# Find out the variables with singular fits
OSI_sf <- map_df(lmerML, ~VarCorr(.x)$Net_pen[1]) %>% 
  gather(key = Gut_segment) %>%
  filter(value == 0)

# Make a list of dataframes containing variables with singular fits
df_sf <- filter(df, Gut_segment %in% OSI_sf$Gut_segment) %>%
  mutate(Gut_segment = factor(Gut_segment, levels = c("PI", "MI"))) %>% 
  split(f = .$Gut_segment)

# Welch's t-test
welch_t <- lapply(df_sf, function(x) t.test(OSI ~ Diet, x)) 

# Check the normality assumption visually via Q-Q plots
qq_welchT <- lapply(
  seq_along(df_sf), 
  function(x) 
  {
    ggqqplot(df_sf[[x]], "OSI", color = "Diet", 
             facet.by = "Diet", title = names(df_sf)[x]) 
    
  }
)

ggexport(filename = "analysis/model_diagnostics/OSI_welch-t_qqplot.pdf",
         plotlist = qq_welchT, nrow = 2, width = 8)
    
# Fit LMMs with REML -----------------------------------------------------------  
# Get the variable not showing a singular fit
df_nsf <- filter(df, !Gut_segment %in% OSI_sf$Gut_segment) %>%
  mutate(Gut_segment = as.character(Gut_segment)) 

# Fit LMMs to "DI"
lmerREML <- lmer(OSI ~ Diet + (1|Net_pen), REML = T, data = df_nsf) 

### Model diagnostics ###

## 1.Residual analysis ##
# Extract residuals
res <- resid(lmerREML)

# Check homogeneity and normality of residuals
pdf(file = "analysis/model_diagnostics/OSI_LMMs_residual.pdf", width = 8) 

par(mfrow = c(2, 2), mar = c(5, 5, 2, 2), cex.lab = 1.5)
    
plot(res, xlab = "Fitted values", ylab = "Residuals") # plot residuals against fitted values
abline(0,0)
    
boxplot(res ~ Diet, data = df_nsf, xlab = "Diet", ylab = "Residuals") # plot residuals against Diet
abline(h = 0, lty = 2)
  
boxplot(res ~ Net_pen, data = df_nsf, xlab = "Net pen", ylab = "Residuals") # plot residuals against Net pen
abline(h = 0, lty = 2)
  
qqnorm(res, main = "Normal Q-Q Plot") # make a Q-Q plot 
qqline(res)
  
par(oma = c(0,0,2,0))
title(main = "DI", outer = T) 
  
dev.off() 

# 2.Influence analysis 
# Based on the residual analysis, we don't need to run influence analysis 

# Get p values of the Diet effect ----------------------------------------------
## Get p values from LMMs using parametric bootstrap 
# Construct a reduced model
lmerML_nofix <- update(lmerML[["DI"]], . ~ . - Diet)

# Create clusters for parallel computing to speed up parametric bootstrap 
nc <- parallel::detectCores() 
clus <- parallel::makeCluster(rep("localhost", nc)) 

# Parametric bootstrap comparisons 
pb <- PBmodcomp(lmerML[["DI"]], lmerML_nofix, seed = 1910, cl = clus)

# Stop the clusters
parallel::stopCluster(clus)

# Gather p values
p_val <- data.frame("Gut_segment" = names(lmerML),
                    "p_values"    = c(welch_t$PI$p.value,
                                      welch_t$MI$p.value,
                                      pb$test["PBtest", "p.value"]))  
                    
# Make Figure 1 ################################################################
# Initial plot
fig <- df %>% 
  ggplot(aes(x = Diet, y = OSI)) +
  geom_boxplot(aes(fill = Diet), outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), 
                     labels = fmt_dcimals(1), 
                     expand = expand_scale(mult = c(0, 0.1))) +  
  labs(y = "Organosomatic index (%)") +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

# Add significant p values to the plot -----------------------------------------
# Make a datafraome for the p-value annotation
ann <- filter(p_val, p_values <= 0.05) %>% # retain rows with p vales smaller than 0.05 for the annotation
       mutate(p_values = formatC(p_values, format = "f", digits = 3), # format digits of p values
              p_values = paste("p = ", p_values), # add label "p =" to p values
              start    = rep("REF", nrow(.)), # start position of p value labeling on x axis
              end      = rep("IM", nrow(.)), # end position of p value labeling on x axis
              y        = rep(0.47, nrow(.)) # position of p value label on y axis
              ) 
  

# Add the significant p values. The warning about the missing aesthetics can be ignored.
fig + geom_signif(data = ann, 
                  aes(xmin = start, 
                      xmax = end, 
                      annotations = p_values, 
                      y_position = y),
                  textsize = 4, 
                  tip_length = 0,
                  manual = T)

ggsave('results/figures/Figure 1.tiff', units = "in", dpi = 300, compression = "lzw")

# Get session info
writeLines(capture.output(sessionInfo()), "analysis/code/organosomatic_index_sessionInfo.txt")
