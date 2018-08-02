library(dplyr) # data manipulation
library(ggplot2) # plotting
library(reshape2) # convert data format (long/wide)
library(stats) # for running fisher's exact test

# import data
df <- read.csv("z:/PhD/Lab/Histology/AquaFl2/R/AquaFly2-Histology-R.csv", header = T, na.strings = c(""))
head(df)

# convert numeric Variable to string
df <- within(df, {
  Net_pen <- as.character(Net_pen)
  Sample_ID <- as.character(Sample_ID)
})

# code for making contingency table from raw data, not used in downstream analysis######################################
df_contingency <- df %>% 
  tidyr::gather("Variable", "Class", 4:ncol(df)) %>% # convert data to the long format
  na.omit() %>% 
  group_by(Diet, Variable, Class) %>% 
  mutate(observation = n()) %>% # # count the number of fish assigned with different Classes of histilogical characteristics evaluated
  dcast(., Variable + Diet ~ Class, value.var="observation") %>% # convert data to the wide format
  arrange(desc(Variable), desc(Diet))

# Plotting #############################################################################################################
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

####################################################################################################################### 
# convert dataframe to the "long" format
df_long <- df %>% 
  tidyr::gather("Variable", "Class", 4:ncol(df)) %>% 
  na.omit()

# Split the data frames by "Variable" 
df_split <- split(df_long, f = df_long$Variable)

# `lapply()` fisher's exaxt test to each data frame in the list and store results as a list  
fisher_test_list <- lapply(df_split, function(x) fisher.test(x$Diet, x$Class))  

# Extract p values, store them in a data frame and apply multiple comparison correction
fisher_test <- data.frame(row.names = NULL,
                          "Histological_characteristic" = names(fisher_test_list),
                          "p_raw" = unlist(lapply(fisher_test_list, function(x) x[[1]]))) # Extract p values of fisher's test from the list

fisher_test$p_adjusted <- p.adjust(fisher_test$p_raw, method = "fdr")

