library(dplyr) # data manipulation
library(ggplot2) # plotting
library(nlme) # running linear mixed model

df <- read.csv("./data/AqFl2_organ_somatic_index.csv", header = T, na.strings = c(""))
df

# convert numeric variable to string
df <- within(df, {
    Sample_ID <- as.character(Sample_ID)
    Net_pen <- as.character(Net_pen)
})

# Reorder varialbes in the desired order
df$Diet <- factor(df$Diet, levels = c("REF", "IM33", "IM66", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))
                
# Exploratory analysis ####################################################################################### 
# Check outliers and data distribution
# tiff(filename = "./results/figures/organosomatic_indices.eps", units = "in", res = 300, compression = "lzw") #File too large. Won't work. 
df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  ggplot(., aes(x = Diet, y = Organ_somatic_index)) + 
  geom_boxplot(aes(fill = Diet)) + 
  facet_wrap(~Gut_segment, ncol = 3, scales = "free_y") + 
  expand_limits(y = 0) + 
  labs(title = "Organosomatic indices", y = "Organosomatic index") +
  theme(axis.text.x=element_text(size = 18),
        axis.text.y=element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        #strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face="bold", size=18),
        legend.text=element_text(face="bold",size=18),
        legend.title=element_text(face="bold",size=20),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5)) 
# dev.off()
ggsave("organosomatic_indices.tiff", units = "in", dpi = 300, compression = "lzw", path = "./results/figures")

# Statistics ##########################################################################################################
# Convert Diet and Target back to character, otherwise the split function won't work correctly
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
  "Gut segment" = names(lme),
  "p values" = unlist(lapply(anova, function(x) x[[4]][[2]]))
)
  
  
# Model diagnostics ---------------------------------------------------------------------------------------------------
# Extract residuals
res <- lapply(lme, function(x) resid(x))

# 1.Linearity and homoskedasticity of residuals
pdf(file = "./results/statistical outputs/resid_plot_organosomatic_indices.pdf", encoding = "CP1253") # turn on pdf graphics device
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

pdf(file = "./results/statistical outputs/qqplot_organosomatic_indices.pdf", encoding = "CP1253")
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
