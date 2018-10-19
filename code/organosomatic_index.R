library(tidyverse) # tidy, manipulate and plot data
library(ggrepel) # labelling

df <- read.csv("./data/AqFl2_organ_somatic_index.csv", header = T, na.strings = c(""))
head(df, n = 20L)
str(df)

# convert numeric variable to string
df <- within(df, {
    Sample_ID <- as.character(Sample_ID)
    Net_pen <- as.character(Net_pen)
})

# Calculate organosomatic index
df <- df %>%
  mutate(Organ_somatic_index = 100*Organ_weight/Body_weight)
  
# Reorder varialbes in the desired order 
df$Diet <- factor(df$Diet, levels = c("REF", "IM33", "IM66", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))
                
# Exploratory analysis ################################################################################################ 
# Violin plot to check data distribution
df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  ggplot(aes(x = Diet, y = Organ_somatic_index)) + 
    geom_violin(aes(fill = Diet), trim = FALSE) +
    facet_wrap(~Gut_segment, nrow = 1, scales = "free_y") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange") + # add mean and SD 
    expand_limits(y = 0) + 
    labs(title = "Organosomatic indices", y = "Organosomatic index") +
    theme(axis.text.x=element_text(size = 18), 
          axis.text.y=element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          strip.text = element_text(face="bold", size=18),
          legend.text=element_text(face="bold",size=18),
          legend.title=element_text(face="bold",size=20),
          plot.title = element_text(size = 20, face = "bold",hjust = 0.5)
          ) 

# Box plot to check outliers
# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Make box plot with dots
df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>%
  group_by(Gut_segment, Diet) %>%
  mutate(outlier = is_outlier(Organ_somatic_index)) %>% # mark outliers within each diet group
  ggplot(aes(x = Diet, y = Organ_somatic_index,label = ifelse(outlier, Sample_ID, NA))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = Net_pen), size = 2, shape = 21, position = position_jitterdodge(0.2)) +
    facet_wrap(~ Gut_segment, nrow = 1, scales="free_y") +
    geom_label_repel() +
    labs(title = "Organosomatic indices", y = "Organosomatic index") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold",hjust = 0.5),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.direction="horizontal", 
          legend.position="bottom") 

# raincloud plot


# Statistics ##########################################################################################################
library(nlme) # running linear mixed model

# Convert Diet and Target back to character
df$Diet <- as.character(df$Diet)
df$Gut_segment <- as.character(df$Gut_segment)

# Split the data frames by "Gut_segment"
df_spl <- df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  split(f = .$Gut_segment)

# lapply() linear mixed model to each data frame in the list 
lme <- lapply(df_spl, function(x) lme(Organ_somatic_index ~ Diet, random = ~1|Net_pen, data = x))

# Get model summary and extract p values of the Diet effect.
summary <- lapply(lme, function(x) summary(x))
anova <- lapply(lme, function(x) anova(x))

p_values <- data.frame(row.names = NULL,
                       "Gut_segment" = names(lme),
                       "p_values" = unlist(lapply(anova, function(x) x[[4]][[2]]))
                       ) %>% 
  arrange(desc(Gut_segment)) 
  
  
# Model diagnostics ---------------------------------------------------------------------------------------------------
# Extract residuals
res <- lapply(lme, function(x) resid(x))

# 1.Linearity and homoskedasticity of residuals
pdf(file = "./results/statistics/resid_plot_organosomatic_indices.pdf") # turn on pdf graphics device

lapply(
  seq_along(res), # add index to the elements in the list
  function(x) 
  {
    plot(res[[x]], main = names(res)[x], ylab = "residuals")
    abline(0,0)
  }
)

dev.off() # turn off the device

# 2.Absence of collinearity. When more than one fixed effects are included, the collinearity shoulb be checked.
# Not applicable to the present dataset, which has only one fixed effect

# 3.Normality of residuals

pdf(file = "./results/statistics/qqplot_organosomatic_indices.pdf")

lapply(
  seq_along(res), 
  function(x) 
  {
    qqnorm(res[[x]], main = paste("Normal Q-Q Plot: ", names(res)[x]))
    qqline(res[[x]])
  }
)

dev.off()

# 4.Absence of influential data points

# Make Figure 1 #######################################################################################################
# Reorder varialbes in the desired order 
df$Diet <- factor(df$Diet, levels = c("REF", "IM33", "IM66", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))

# Make a function for setting the number of decimals of y axis
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

# Initial plot
fig1 <- df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>%
  ggplot(aes(x = Diet, y = Organ_somatic_index)) +
    geom_boxplot(aes(fill = Diet), outlier.shape = NA, show.legend = FALSE) +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    facet_wrap(~ Gut_segment, nrow = 1, scales="free_y") +
    scale_y_continuous(limits = c(0, NA), labels = fmt_dcimals(1)) +
    labs(y = "Organosomatic index") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 16)
          ) 

# Add significant p values to the plot
# First, make a datafraome for the annotation
anno <- p_values %>%
  arrange(Gut_segment) %>%
  filter(p_values <= 0.05) %>%
  mutate(start = rep("REF", nrow(.))) %>%
  mutate(end = rep("IM100", nrow(.))) %>%
  mutate(y = c(0.47)) 
  
# Add the p values. The warning about the missing aesthetics can be ignored.
library(ggsignif) # adding significant levels to plot

fig1 + geom_signif(data = anno, 
                   aes(xmin = start, 
                       xmax = end, 
                       annotations = formatC(as.numeric(p_values), digits = 3), 
                       y_position = y),
                   textsize = 4, 
                   tip_length = 0.01,
                   manual = TRUE)
# Save the figure
ggsave("Figure 1. Organosomatic indices.tiff", 
       units = "in", 
       dpi = 300, 
       compression = "lzw", 
       path = "./results/figures")