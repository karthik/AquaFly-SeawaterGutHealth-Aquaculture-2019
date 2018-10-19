library(dplyr)
library(ggplot2)
library(ggrepel)
library(broom)

df <- read.csv("C:/Users/ljt89/OneDrive - Norwegian University of Life Sciences/Code demo/AqFl2_qPCR_PI_Ref.csv", 
               header =T, na.strings=c(""))

df <- read.csv("Z:/PhD/Lab/qPCR/AqFl2/R/AqFl2_qPCR_PI_Ref.csv", header =T, na.strings=c(""))

head(df)

# Tidy data and calculate new variables
df1 <- df %>% 
  select(-Comments) %>% # remove the "Comments" column
  na.omit() %>% # remove rows with NA value
  arrange(Ref_gene, desc(Diet)) %>% # order rows by "Ref_gene" and "Diet"  
  mutate(Quantity = `^`(PE, -Cq)) %>% # calculate mRNA quantity
  group_by(Ref_gene, Diet)  %>%
  mutate(Mean_diet = mean(Quantity), SD = sd(Quantity)) # calculate average mRNA quantity per diet

# convert "Diet" and "Target" to factor variables with the desired order
df1$Diet <- factor(df1$Diet, levels = c("REF", "IM100"))

## Plotting ###########################################################################################################
## Point diagram ------------------------------------------------------------------------------------------------------
# make Sample_ID as factor variable
df1$Sample_ID <- factor(df1$Sample_ID, levels = unique(df1$Sample_ID))

# make a dataframe from which grand mean of mRNA quantity can be plotted
gm <- df1 %>% 
  group_by(Ref_gene) %>%
  summarize(mean = mean(Quantity))

# make point diagram 
p_dot <- df1 %>% 
  ggplot(aes(x = Sample_ID, y = Quantity, group = 1)) +
  geom_point(aes(color = Diet)) +
  geom_line() +
  facet_wrap(~ Ref_gene, ncol = 1,scales = "free_y") +
  theme_bw() +
  geom_hline(aes(yintercept = mean), gm, color = "blue") # The blue line shows the grand mean. 

## Boxplot ------------------------------------------------------------------------------------------------------------
# convert Sample_ID as character, otherwise the labelling by Sample_ID won't work correctly
df1$Sample_ID <- as.character(df1$Sample_ID)

# Define a function for identifying data outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Make boxplot using standardized mRNA quantity
p_box <- df1 %>% 
  group_by(Ref_gene) %>%
  mutate(Quantity_std = scale(Quantity))  %>% # the values are subtracted by the mean (centered) and divided by the standard deviation (scaled)
  mutate(outlier = is_outlier(Quantity_std)) %>%
  ggplot(aes(x=Ref_gene, y=Quantity_std, label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Diet),size = 3, shape = 21, position = position_jitterdodge(0.2)) +
  geom_label_repel() +
  theme_bw()

## Barplot ------------------------------------------------------------------------------------------------------------
p_bar <- df1 %>%
  ggplot(aes(x = Diet, y =  Mean_diet, fill = Diet)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin = Mean_diet, ymax = Mean_diet + SD), 
                size = 0.3, 
                width= 0.2) + 
  facet_wrap(~ Ref_gene, nrow = 1,scales = "free_y") +
  theme_bw() +
  expand_limits(y = 0) +
  ylab("mRNA quantity") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

## Rank reference genes by coefficient of variation and F statistic ###################################################
cv <- df1 %>% 
  group_by(Ref_gene) %>%
  summarize(Cq_Max = max(Cq), 
            Cq_Min = min(Cq), 
            Cq_range = max(Cq) - min(Cq), 
            Quantity_mean = mean(Quantity), 
            Quantity_SD = sd(Quantity), 
            Quantity_cv = 100 * Quantity_SD / Quantity_mean) %>%
  mutate(cv_rank = dense_rank(Quantity_cv)) %>%
  arrange(Ref_gene)


## Rank reference genes by F statistics
F_statistic <- df1 %>% 
  group_by(Ref_gene) %>% 
  do(tidy(aov(Quantity ~ Diet, data = .))) %>%
  na.omit() %>%
  select(Ref_gene, statistic) %>%
  rename(Quantity_F = statistic) %>%
  as.data.frame() %>%
  mutate(F_rank = dense_rank(Quantity_F)) %>%
  arrange(Ref_gene) %>%
  tbl_df()

## Summary
Summary <- full_join(cv, F_statistic, by = "Ref_gene")

