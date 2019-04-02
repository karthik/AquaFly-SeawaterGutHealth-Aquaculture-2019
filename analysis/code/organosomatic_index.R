library(tidyverse) # tidy, manipulate and plot data
library(cowplot) #  a ggplot2 add-on to provide a publication-ready theme
library(ggsignif) # add p values to plots
library(lme4) # linear mixed effect modelling
library(nlme) # linear mixed effect modelling, allowing heterogeneity
library(influence.ME) # detect influential observations in the linear mixed effect model

# Import and tidy data #########################################################
df <- read.csv("./data/raw_data/AqFl2_organosomatic_index.csv", 
               header = T, 
               na.strings = c(""), 
               stringsAsFactors = FALSE)

head(df, n = 20L)
str(df)

# Calculate organosomatic index (OSI)
df <- df %>%
  mutate(OSI = 100 * Organ_weight / Body_weight)

# Exploratory analysis ######################################################### 
# Define a function for setting the number of decimals shown on axis
fmt_dcimals <- function(decimals = 0){
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}

# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Convert numeric variable to character
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Convert character to factor, specifying the desired order to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))

# Violin plot to check data distribution ---------------------------------------
ggplot(df, aes(x = Diet, y = OSI)) + 
  geom_violin(aes(fill = Diet), trim = FALSE) +
  facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
  stat_summary(fun.data = "mean_sdl", 
               fun.args = list(mult = 1), 
               geom = "pointrange") + # add mean and SD 
  scale_y_continuous(limits = c(0, NA), 
                     labels = fmt_dcimals(1), 
                     expand = expand_scale(mult = c(0, 0.1))) + 
  labs(title = "Organosomatic indices", y = "%") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('./analysis/exploratory_analysis/OSI_violin.pdf', units = "in", dpi = 300)

# Box plot to check outliers and explore for cluster effect, i.e., net pen -----
df %>%
  group_by(Gut_segment, Diet) %>%
  mutate(outlier = is_outlier(OSI)) %>% # mark outliers within each diet for each gut segment
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
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('./analysis/exploratory_analysis/OSI_boxPlot.pdf', units = "in", dpi = 300)

# Statistics ###################################################################
# Split the data frame by "Gut_segment"
df_spl <- split(df, f = df$Gut_segment)

# Run linear mixed model for each gut_segment
lmer <- lapply(df_spl, 
               function(x) 
               lmer(OSI ~ Diet + (1|Net_pen), data = x))                                         
              
# Model diagnostics ------------------------------------------------------------
# 1.Residual analysis 
# Extract residuals
res <- lapply(lmer, function(x) resid(x))

# Check homogeneity and normality of residuals
pdf(file = "./analysis/model_diagnostics/OSI_residual.pdf", width = 8) 

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
  
  par(oma=c(0,0,2,0))
  title(main = names(res)[x], outer = TRUE)}, 
  x = seq_along(res), 
  y = df_spl
)

dev.off() 

# MI showed heterogeneity with one potential outlier

# 2.Influence analysis 
# Based on the residual analysis, we'll run influence analysis for MI.
# As the influence() function cann't be applied to an element in a list, we'll
# have to fit the model again
df_mi <- filter(df, Gut_segment == "MI") 
lmer_mi <- lmer(OSI ~ Diet + (1|Net_pen), data = df_mi) 

# Detect influence observations using Cook's distance and sigtest
estex.obs <- influence(lmer_mi, obs = TRUE)

plot(estex.obs, which = "cook", 
     cutoff = 4/36, sort = TRUE, 
     xlab = "Cook´s Distance",
     ylab = "Row Number")

sigtest(estex.obs, test = -1.96)$Diet

# Sample 56 has influence on the parameter estimation but deosn't change the 
# significant level of Diet effect.

# Update models ----------------------------------------------------------------
# To deal with heterogeneity, we'll use varIdent() function from the nlme package 
lme_mi <- lme(OSI ~ Diet,
              random  = ~ 1|Net_pen, 
              weights = varIdent(form = ~ 1|Net_pen), # allow different variance for different net pen
              data    = df_mi)

# Compare to the original model
AIC(lmer[["MI"]], lme_mi)

# Get p values of the Diet effect ----------------------------------------------
# Refit the models using maximum likelihood estimation
lmer_ML <- lapply(lmer[-2], function(x) update(x, . ~ ., REML = FALSE))
lme_mi_ML <- update(lme_mi, . ~ ., method = "ML")

# Likelihood ratio test of fixed effect using the drop1() function
lmer_fixef <- lapply(lmer_ML, function(x) drop1(x, test = "Chi"))
lme_mi_fixef <- drop1(lme_mi_ML, test = "Chi")

# Get p values of the Diet effect
p_val <- data.frame(row.names = NULL,
                    "Gut_segment" = names(lmer),
                    "p_values" = c(lmer_fixef[["PI"]]["Diet", "Pr(Chi)"],
                                   lme_mi_fixef["Diet", "Pr(>Chi)"],
                                   lmer_fixef[["DI"]]["Diet", "Pr(Chi)"])
                    )  
                    
# Make Figure 1 ################################################################
# Initial plot
fig1 <- df %>% 
  ggplot(aes(x = Diet, y = OSI)) +
  geom_boxplot(aes(fill = Diet), outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), 
                     labels = fmt_dcimals(1), 
                     expand = expand_scale(mult = c(0, 0.1))) +  
  labs(y = "Organosomatic index (%)") +
  theme_cowplot() +
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
fig1 + geom_signif(data = ann, 
                   aes(xmin = start, 
                       xmax = end, 
                       annotations = p_values, 
                       y_position = y),
                   textsize = 4, 
                   tip_length = 0,
                   manual = TRUE)

ggsave('./results/figures/Figure 1.tiff', units = "in", dpi = 300, compression = "lzw")

# Get session info
writeLines(capture.output(sessionInfo()), "analysis/code/organosomatic_index_sessionInfo.txt")
