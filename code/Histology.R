library(tidyverse) # tidy, manipulate and plot data
library(cowplot) # a ggplot2 add-on to provide a publication-ready theme 
library(gridExtra) # combine figures
library(ggsignif) # # add p values to plot
library(ordinal) # run ordered logistic regression

# Import and tidy data ################################################################################################
df <- read.csv("./data/raw_data/AqFl2_histology.csv", header = T, na.strings = c(""), stringsAsFactors = FALSE)
head(df, n = 20L) 
str(df)

# Tidy data
df <- df %>% 
  gather("Variable", "Rank", 4:ncol(df)) %>% # convert data to from "wide" to "long" format
  na.omit() %>% # remove rows with "NA"
  mutate(Split = Variable, # use as an index for splitting the dataframe for running ordered logistic regression
         Gut_segment = gsub("\\_.*", "", Variable), # create a new column specifying the gut segment
         Variable = gsub("*.I\\_", "", Variable)) # edit the "Variable" column to remove the "PI_/MI_/DI_" 
                                                  # *. repalce anything before; replace anything after .*  
                                                  # to replace punctuation characters, add double space (\\) prefix

# Exploratory analysis ################################################################################################
# Convert numeric variables to character
df <- within(df, {
  Net_pen <- as.character(Net_pen)
  Sample_ID <- as.character(Sample_ID)
})

# Convert character to factor, specifying the desired orders to be shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))
df$Rank <- factor(df$Rank, levels = c("Normal", "Mild", "Moderate", "Marked", "Severe"))
df$Variable <- factor(df$Variable, levels = c("hpv", "snv", "smc", "lpc", "mfh"))

# Make a contingency table to get an idea of data distribution
df_ctg <- df %>% 
  group_by(Gut_segment, Diet, Variable, Rank) %>% 
  mutate(observation = n()) %>% # count the number of fish assigned with different ranks of histilogical characteristics 
  reshape2::dcast(Gut_segment + Variable + Diet ~ Rank, value.var = "observation") %>% # convert data to the wide format
  arrange(Gut_segment, Variable, Diet)

# Make histogram
ggplot(df, aes(Rank, fill = Net_pen)) +
  geom_histogram(#width = 1, # set width to 1 if removing the gap between bars is desired
                 stat="count") +
  facet_grid(Gut_segment + Diet ~ Variable) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Dark2")

ggsave('./results/exploratory_analysis/histology_histogram.pdf', units = "in", dpi = 300)

# Statistics ##########################################################################################################
# Aggregate the levels of ranks, reduced to 3 levels  
df <- df %>%
  mutate(Rank_agr = gsub("Normal", "1", Rank), 
         Rank_agr = gsub("Mild|Moderate", "2", Rank_agr),
         Rank_agr = gsub("Marked|Severe", "3", Rank_agr))

# Convert ranks to ordered factor
df$Rank <- ordered(df$Rank, levels = c("Normal", "Mild", "Moderate", "Marked", "Severe"))
df$Rank_agr <- ordered(df$Rank_agr, levels = c("1", "2", "3"))

# Split the dataframe by "Variable" 
df_spl <- split(df, f = df$Split)

# Run ordered logistic regression for each dataframe and compare models
mod_ful <- lapply(df_spl, function(x) clmm(Rank ~ Diet + (1|Net_pen), data = x, Hess = TRUE, nAGQ = 10))
mod_sub <- lapply(df_spl, function(x) clm(Rank ~ Diet, data = x))

smr_ful <- lapply(mod_ful, function(x) summary(x))
smr_sub <- lapply(mod_sub, function(x) summary(x))

mod_select <- data.frame(row.names = NULL, 
                         "Variable" = names(smr_ful),
                         "AIC_ful" = unlist(lapply(smr_ful, function(x) x[["info"]][["AIC"]])),
                         "AIC_sub" = unlist(lapply(smr_sub, function(x) x[["info"]][["AIC"]])),
                         "Hessian_ful" = unlist(lapply(smr_ful, function(x) x[["info"]][["cond.H"]])),
                         "Hessian_sub" = unlist(lapply(smr_sub, function(x) x[["info"]][["cond.H"]]))
)

# Run anova to get p values of diet effect
anv_sub <- lapply(mod_sub, function(x) anova(x))
p_val <- data.frame(row.names = NULL, 
                    "Variable" = names(anv_sub),
                    "p_raw" = unlist(lapply(anv_sub, function(x) x[[3]]))
                    )
p_val <- p_val %>% 
  mutate(Gut_segment = gsub("\\_.*", "", Variable), 
         Variable = gsub("*.I\\_", "", Variable)) %>%
  group_by(Gut_segment) %>%
  nest() %>%
  mutate(p_adj = map(data, ~p.adjust(.x$p_raw, method = "holm"))) %>%
  unnest()
         
# Model assumption checking via nominal test
nominal_test <- lapply(mod_sub, function(x) nominal_test(x))

# Make Figure 2 #######################################################################################################
# Calculate new variable "Percent"for plotting stacked bar plot
df_bar <- df %>% 
  group_by(Gut_segment, Diet, Variable, Rank) %>% 
  mutate(observation = n()) %>% # number of fish assigned with different ranks
  ungroup() %>% 
  group_by(Gut_segment, Diet, Variable) %>% 
  mutate(Percent = 100 * observation/n()) %>% # calculate the percentage of different ranks within the same diet 
  ungroup() %>% 
  group_by(Gut_segment, Diet, Variable, Rank, Percent) %>% 
  summarize()

# Define color schemes and choose one you prefer
my_col1 = c(Normal = "royalblue2", Mild = "tan1", Moderate = "tomato", Marked = "orangered", Severe = "red3")
my_col2 = c(Normal = "royalblue2", Mild = "peachpuff1", Moderate = "tan1", Marked = "tomato", Severe = "red3")

# Use different fecet labels by changing the underlying factor level names
levels(df_bar$Variable) <- c("Hypervacuolization", 
                             "Supranuclear vacuolization", 
                             "Submucosal cellularity",
                             "Lamina propria cellularity",
                             "Mucosal fold height")

# Make stacked barplot for each gut segment
p1 <- df_bar %>%
  filter(Gut_segment == "PI") %>%
  ggplot(aes(Diet, Percent, fill = forcats::fct_rev(Rank))) + # forcats::fct_rev() reverses stacked bars
  geom_bar(stat = "identity") +
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col2) + 
  scale_y_continuous(limits = c(0, 105), breaks = 0:5*20, expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "PI", y = "%") +
  guides(fill = guide_legend(title = "Rank")) + # use "Rank" as legend title
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0), 
        legend.position = "none")

p2 <- df_bar %>%
  filter(Gut_segment == "MI") %>%
  ggplot(aes(Diet, Percent, fill = forcats::fct_rev(Rank))) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col2, drop = FALSE) + # drop = FALSE forces legend to show all categories
  scale_y_continuous(limits = c(0, 105), breaks = 0:5*20, expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "MI", y = "%") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "none")

p3 <- df_bar %>%
  filter(Gut_segment == "DI") %>%
  ggplot(aes(Diet, Percent, fill = forcats::fct_rev(Rank))) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col2) + 
  scale_y_continuous(limits = c(0, 105), breaks = 0:5*20, expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "DI", y = "%") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0), 
        legend.position = "none")

# Add significant p values to the plots -----------------------------------------------------------
# Make a datafraome for the p-value annotation
ant <- p_val %>%
  filter(Gut_segment == "PI", Variable == "hpv") %>%
  mutate(Variable = gsub("hpv", "Hypervacuolization", Variable)) %>%
  select(Variable, p_adj) %>%
  mutate(start = "REF",
         end = "IM100",
         y = 100)

# Add the p values. The warning about the missing aesthetics can be ignored
p1 + geom_signif(data = ant,
                 aes(xmin = start, 
                     xmax = end, 
                     annotations = formatC(as.numeric(p_adj), format = "f", digits = 3), 
                     y_position = y),
                 textsize = 4, 
                 tip_length = 0,
                 manual = TRUE)

# P values cann't be added. Error message: Error in check_factor(f) : object 'Rank' not found.
# P values can be added if the barplots are not stacked by Rank, using ggsignif or geom_text + geom_segment

# Combine figures ---------------------------------------------------------------------------------
# Extract legend from one of the figures and use it as the shared legend
legend <- get_legend(p1 + theme(legend.position = "right"))

# Make a list of grobs 
gl <- list("1" = p1, "2" = p2, "3" = p3, "4" = legend)  

# Make a layout matrix to gruide the layout of figures
lom <- cbind(c(1, 2, 3),
             c(1, 2, 3), 
             c(1, 2, 3), 
             c("NA", 4, 3)
            )

# Use gridExtra::grid.arrange to combine the figures and export the final figure as Figure 2
tiff('./results/figures/Figure 2.tiff', 
     compression = "lzw",
     units = "in", 
     res = 300, 
     height = 12,
     width = 10)

grid.arrange(grobs = gl, layout_matrix = lom)

dev.off()

# Get session info
writeLines(capture.output(sessionInfo()), "./code/histology_sessionInfo.txt")