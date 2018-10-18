library(dplyr) # data manipulation
library(ggplot2) # plotting
library(reshape2) # convert data format (long/wide)
library(ordinal) # for running ordered logistic regression

# import data
df <- read.csv("./data/AqFl2_histology.csv", header = T, na.strings = c(""))
df

# convert numeric Variable to string
df <- within(df, {
  Net_pen <- as.character(Net_pen)
  Sample_ID <- as.character(Sample_ID)
})

# code for making contingency table from raw data. Not used in downstream analysis#####################################
df_contingency <- df %>% 
  tidyr::gather("Variable", "Class", 4:ncol(df)) %>% # convert data to the long format
  na.omit() %>% 
  group_by(Diet, Variable, Class) %>% 
  mutate(observation = n()) %>% # # count the number of fish assigned with different Classes of histilogical characteristics evaluated
  dcast(., Variable + Diet ~ Class, value.var="observation") %>% # convert data to the wide format
  arrange(desc(Variable), desc(Diet))

# Plotting ############################################################################################################
# edit and summarize data for plotting
df_plot <- df %>% 
  tidyr::gather("Variable", "Class", 4:ncol(df)) %>% # convert dataframe to the "long" format
  na.omit() %>% 
  group_by(Diet, Variable, Class) %>% 
  mutate(observation = n()) %>% # count the number of fish assigned with different Classes of histilogical characteristics evaluated 
  ungroup() %>% 
  group_by(Diet, Variable) %>% 
  mutate(percent = observation/n()) %>% # calculate the percentage of different clesses within the same diet 
  ungroup() %>% 
  group_by(Diet, Variable, Class, percent) %>% 
  summarize()

# reoder the 'Class' and 'Variables' as desired
df_plot$Class <- factor(df_plot$Class, levels = c("Severe", "Marked", "Moderate", "Mild", "Normal"))
df_plot$Variable <- factor(df_plot$Variable, levels = c("PI_hpv", "PI_lpc", "PI_smc", "MI_hpv", "MI_lpc", "MI_smc", 
                                                        "DI_snv", "DI_lpc", "DI_smc", "DI_mfh"))
# make stacked bar plot
# define the color scheme
my_col = c(Normal = "royalblue2", Mild = "tan1", Moderate = "tomato", Marked = "orangered", Severe = "red3")

p_histo <- ggplot(df_plot, aes(Diet, percent, fill = Class)) + 
  geom_bar(stat = "identity") + 
  facet_grid(~Variable) + 
  ggtitle("Histology") + 
  scale_x_discrete(limits = c("REF", "IM100")) + 
  scale_fill_manual(values = my_col) + 
  scale_y_continuous(limits = c(0, 1.02), breaks = 0:5*0.2, expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(size = 15, 
                                  face = "bold", hjust = 0.5))

# statistics ##########################################################################################################
# convert dataframe to the "long" format
df_long <- df %>% 
  tidyr::gather("Variable", "Class", 4:ncol(df)) %>% 
  na.omit() 
 
# convert "Class" to ordered factor
df_long$Class <- ordered(df_long$Class, levels = c("Normal", "Mild", "Moderate", "Marked", "Severe"))

# Split the data frames by "Variable" 
df_spl <- split(df_long, f = df_long$Variable)

# run cumulative link mixed model for each dataframe, treating "Diet" as fixed effect and "Net_pen" as random effect
mod_ful <- lapply(df_spl, function(x) clmm(Class ~ Diet + (1|Net_pen), data = x))

mod_sub <- lapply(df_spl, function(x) clm(Class ~ Diet, data = x))

summary_ful <- lapply(mod_ful, function(x) summary(x))
summary_sub <- lapply(mod_sub, function(x) summary(x))
anova_sub <- lapply(mod_sub, function(x) anova(x))

mod_select <- data.frame(row.names = NULL, 
                         "Variable" = names(summary_ful),
                         "AIC_ful" = unlist(lapply(summary_ful, function(x) x[["info"]][["AIC"]])),
                         "AIC_sub" = unlist(lapply(summary_sub, function(x) x[["info"]][["AIC"]])),
                         "Hessian_ful" = unlist(lapply(summary_ful, function(x) x[["info"]][["cond.H"]])),
                         "Hessian_sub" = unlist(lapply(summary_sub, function(x) x[["info"]][["cond.H"]]))
)

p_values_sub <- data.frame(row.names = NULL, 
                           "Variable" = names(anova_sub),
                           "p_raw" = unlist(lapply(anova_sub, function(x) x[[3]]))
)

sig <- p_values_sub[8:10, ] %>%
  mutate(p_adjusted =  p.adjust(p_raw, method = "bonferroni"))
  
# Model assumption checking
nominal_test_sub <- lapply(mod_sub, function(x) nominal_test(x))



