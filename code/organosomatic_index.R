library(tidyverse) # tidy, manipulate and plot data
library(cowplot) #  a ggplot2 add-on to provide a publication-ready theme  
library(nlme) # run linear mixed model
library(ggsignif) # add p values to plot

# Import and tidy data ################################################################################################
df <- read.csv("./data/raw_data/AqFl2_organosomatic_index.csv", 
               header = T, 
               na.strings = c(""), 
               stringsAsFactors = FALSE)

head(df, n = 20L)
str(df)

# Calculate organosomatic index (OSI)
df <- df %>%
  mutate(OSI = 100 * Organ_weight / Body_weight)
  
# Exploratory analysis ################################################################################################ 
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
df$Diet <- factor(df$Diet, levels = c("REF", "IM33", "IM66", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))

# Violin plot to check data distribution ----------------------------------------------------------
df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% # select data for "REF" and "IM100"
  ggplot(aes(x = Diet, y = OSI)) + 
    geom_violin(aes(fill = Diet), trim = FALSE) +
    facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange") + # add mean and SD 
    scale_y_continuous(limits = c(0, NA), labels = fmt_dcimals(1), expand = expand_scale(mult = c(0, 0.1))) + 
    labs(title = "Organosomatic indices", y = "%") +
    theme_cowplot() +
    scale_fill_brewer(palette = "Dark2")

ggsave('./results/exploratory_analysis/OSI_violin.pdf', units = "in", dpi = 300)

# Box plot to check outliers and explore for cluster effect, i.e., net pen ------------------------

df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>%
  group_by(Gut_segment, Diet) %>%
  mutate(outlier = is_outlier(OSI)) %>% # mark outliers within each diet for each gut segment
  ggplot(aes(x = Diet, y = OSI, label = ifelse(outlier, Sample_ID, NA))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = Net_pen), size = 2, shape = 21, position = position_jitterdodge(0.2)) +
    facet_wrap(~ Gut_segment, nrow = 1, scales="free_y") +
    ggrepel::geom_label_repel() + # label outliers using Sample_ID
    scale_y_continuous(limits = c(0, NA), labels = fmt_dcimals(1), expand = expand_scale(mult = c(0, 0.1))) +
    labs(title = "Organosomatic indices", y = "%") +
    theme_cowplot() +
    scale_fill_brewer(palette = "Dark2")

ggsave('./results/exploratory_analysis/OSI_boxPlot.pdf', units = "in", dpi = 300)
# For sanity, check sampling records and raw data of outliers for any warning notes or mistakes

# Statistics ##########################################################################################################
# Split the data frame by "Gut_segment"
df_spl <- df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  split(f = .$Gut_segment)

# Run linear mixed model for each data frame using lapply()
lme <- lapply(df_spl, function(x) lme(OSI ~ Diet, random = ~1|Net_pen, varPower(form = ~ fitted(.)), data = x))

# Get model summary and extract p values of the Diet effect
smr <- lapply(lme, function(x) summary(x))
anv <- lapply(lme, function(x) anova(x))

p_val <- data.frame(row.names = NULL,
                    "Gut_segment" = names(lme),
                    "p_values" = unlist(lapply(anv, function(x) x[[4]][[2]]))
                    ) %>% 
  arrange(desc(Gut_segment))    
  
# Model diagnostics -------------------------------------------------------------------------------
# Extract residuals
res <- lapply(lme, function(x) resid(x))

# 1.Linearity and homoskedasticity of residuals
pdf(file = "./results/statistics/OSI_resid_plot.pdf") 

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
pdf(file = "./results/statistics/OSI_qqplot.pdf")

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

# Visual inspection of residual plot and qq plot did not reveal any obvious deviations from homoscedasticity or 
# normality. There's, however, one leverage point in MI. From the boxplot we can tell that leaving out this data 
# point won't affect the statistical conclusion (significant or not). Hence, there's no need to run a new test.

# Make Figure 1 #######################################################################################################
# Initial plot
fig1 <- df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>%
  ggplot(aes(x = Diet, y = OSI)) +
    geom_boxplot(aes(fill = Diet), outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    facet_wrap(~ Gut_segment, nrow = 1, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA), labels = fmt_dcimals(1), expand = expand_scale(mult = c(0, 0.1))) +  
    labs(y = "Organosomatic index (%)") +
    theme_cowplot() +
    scale_fill_brewer(palette = "Dark2")

# Add significant p values to the plot ------------------------------------------------------------
# Make a datafraome for the p-value annotation
ant <- p_val %>%
  arrange(Gut_segment) %>%
  filter(p_values <= 0.05) %>% # retain rows with p vales smaller than 0.05 for the annotation
  mutate(start = rep("REF", nrow(.))) %>% # start position of p value labeling on x axis
  mutate(end = rep("IM100", nrow(.))) %>% # end position of p value labeling on x axis
  mutate(y = c(0.47)) # check fig1 and manually decide the y intercerpt for showing p value
  
# Add the significant p values. The warning about the missing aesthetics can be ignored.
fig1 + geom_signif(data = ant, 
                   aes(xmin = start, 
                       xmax = end, 
                       annotations = formatC(as.numeric(p_values), format = "f", digits = 4), # fromat p value digits
                       y_position = y),
                   textsize = 4, 
                   tip_length = 0,
                   manual = TRUE)

ggsave('./results/figures/Figure 1.tiff', units = "in", dpi = 300, compression = "lzw")

# Get session info
writeLines(capture.output(sessionInfo()), "./code/organosomatic_index_sessionInfo.txt")