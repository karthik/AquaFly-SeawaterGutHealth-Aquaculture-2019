library(dplyr) # data manipulation
library(ggplot2) # plotting
library(ggrepel) # labelling
library(reshape2) # convert data format (long/wide)
library(NMF) # makeing heatmap
library(lme4) ## running linear mixed model
library(broom) # tidy statistical outputs

# start the session in greek so that greek letters can be recognised
Sys.setlocale(category = "LC_ALL", locale = "Greek")

df <- read.csv("C:/Users/ljt89/OneDrive - Norwegian University of Life Sciences/Code demo/AqFl2_qPCR_DI_Target.csv", 
               header = T, na.strings=c(""), encoding = "UTF-8")

df <- read.csv("Z:/PhD/Lab/qPCR/AqFl2/R/AqFl2_qPCR_PI_Target.csv", header =T, na.strings=c(""), encoding = "UTF-8")

head(df)

# Tidying data and calculate new variables ###########################################################################
# data filtering and calculating new variables
df1 <- df %>% 
  rename(Sample_ID = X.U.FEFF.Sample_ID) %>% # remove "X.U.FEFF." from the column name, which are the bytes reserved to tell that the file is in UTF-8 format.
  select(-Comments) %>% # remove the "Comments" column
  filter(!Target %in% c("il6", "tnfα")) %>% # remove il6 and tnfα from analysis as their expression levels are too low to be reliably quantified
  na.omit() %>% # remove rows with NA value
  mutate(NE_gapdh = `^`(Target_PE,  -Target_Cq) / `^`(gapdh_PE, -gapdh_Cq), # expression normalzied to gapdh
         NE_ranpo2 = `^`(Target_PE,  -Target_Cq) / `^`(rnapo2_PE, -rnapo2_Cq), # expression normalzied to rnapo2
         NE_hprt1 = `^`(Target_PE,  -Target_Cq) / `^`(hprt1_PE, -hprt1_Cq), # expression normalzied to hprt1
         NE_actb = `^`(Target_PE,  -Target_Cq) / `^`(actb_PE, -actb_Cq), # expression normalzied to actb
         MNE = (NE_gapdh * NE_ranpo2 * NE_hprt1 * NE_actb)^(1/4)) # geometric mean of normalzied expression

head(df1)

# Plotting ###########################################################################################################
# Before plotting, first, convert numeric variable to string
df1 <- within(df1, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Then, convert "Diet" and "Target" to factor variables with the desired order
df1$Diet <- factor(df1$Diet, levels = c("REF", "IM100"))

df1$Target <- factor(df1$Target, levels = c("myd88","il1β","il6", "il8", "tnfα", "cd3γδ", "cd8β", "mhcI", "ifnγ",
                                            "il17a", "foxp3", "il10", "tgfβ1", "il4", "plin2", "mtp", "apoa1", "apoa4", 
                                            "apob", "chk", "pcyt1a", "fabp2b", "aqp8ab", "cldn15", "cldn25b", "zo1",
                                            "cdh1", "muc2",  "cyp1a1", "mta", "hsp70",  "sod1", "cat", "casp6", "pcna", 
                                            "mmp13")) 

df1$Target_func <- factor(df1$Target_func, levels = c("Immune_modulation", "lipid_droplet_turnover", "barrier_function", "xenobiotic_matebolism"))                                             

# Make boxplot --------------------------------------------------------------------------------------------------------
# First, define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}
# Plotting
p_box <- df1 %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group
  ggplot(aes(x=Diet, y=MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen),size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, ncol = 4,scales="free_y") +
  geom_label_repel() +
  theme_bw() +
  ggtitle("Gene expression profile in PI") +
  theme(axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        legend.direction="horizontal", 
        legend.position="bottom") 

# Make heatmap -------------------------------------------------------------------------------------------------------
# convert the data format from "long" to "wide"
df_wide <- dcast(df1, Sample_ID + Diet + Net_pen ~ Target, value.var="MNE") %>%
  arrange(Diet, Net_pen) 

# make a numeric matrix for heatmap plotting
dm  <- select(df_wide, -Sample_ID, -Diet, -Net_pen) %>%
  data.matrix()
# assign Sample_ID as row name and transpose the matrix
rownames(dm) <- df_wide[,1]
dm_t <- t(dm)
# make column and row annoation 
annCol <- select(df_wide, Diet, Net_pen)

x <- c("Immune_modulation", "lipid_droplet_turnover", "barrier_function", "xenobiotic_matebolism")
times <- c(12, 8, 6, 8)
func <- rep(x, times = times)
annRow <- data.frame(function_category = func)

# plotting
p_heatmap <- aheatmap(dm_t, 
                      border_color = "grey60",
                      Rowv = NA,
                      Colv = NA,
                      scale = "row",  
                      annCol = annCol,
                      annRow = annRow,
                      fontsize = 14,
                      cellwidth = 22,
                      cellheight = 16)


# Make Barplot showing fold change of gene expression -----------------------------------------------------------------
df4 <- df1 %>%
  filter(Diet == "REF") %>%
  group_by(Target) %>%
  mutate(MNE_group_mean_REF = mean(MNE)) %>%
  ungroup() %>%
  select(Target_func, Target, MNE_group_mean_REF) %>%
  mutate(rowName = row.names(.)) 

df5 <- df1 %>%
  filter(Diet == "IM100") %>%
  rename(MNE_IM100 = MNE) %>%
  mutate(rowName = row.names(.)) %>%
  select(rowName, MNE_IM100) %>%
  full_join(df4, ., by = "rowName") %>%
  mutate(Fold_change = MNE_IM100 / MNE_group_mean_REF) %>%
  group_by(Target_func, Target) %>%
  summarize(N = n(), Mean_fold_change = mean(Fold_change), SD = sd(Fold_change), SE = SD / N^(1/2))

df2 <- df1 %>%
  group_by(Target_func, Target, Diet) %>%
  summarize(MNE_group_mean = mean(MNE)) %>%
  filter(Diet == "REF") %>%
  select(Target_func, Target, MNE_group_mean) %>%
  rename(MNE_group_mean_REF = MNE_group_mean)

df3 <- df1 %>%
  group_by(Target, Diet) %>%
  summarize(MNE_group_mean = mean(MNE)) %>%
  filter(Diet == "IM100") %>%
  select(Target, MNE_group_mean) %>%
  rename(MNE_group_mean_IM100 = MNE_group_mean) %>%
  full_join(df2, ., by = "Target") %>%
  mutate(Fold_change = MNE_group_mean_IM100 / MNE_group_mean_REF)

p_bar <- ggplot(df5, aes(x = Target, y = Mean_fold_change, fill = Target_func)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  coord_flip() +
  geom_errorbar(aes(ymin = Mean_fold_change, ymax = Mean_fold_change + SE), 
                size = 0.3, 
                width= 0.2) +
  scale_y_continuous(limits = c(0, 2), breaks = 0:10*0.2, expand = c(0, 0)) +
  geom_hline(yintercept = 1, color = "blue", linetype = 2) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_blank())

# linear mixed model assumption checking ############################################################################


# run linear mixed model and tidy model output ######################################################################
lmer <- df1 %>% 
  group_by(Target) %>%
  do(tidy(aov(Target_Cq ~ Diet , data = .))) %>%
  rename(F_statistic = statistic) %>%
  select(Gene_Name, F_statistic) %>%
  mutate(F_rank = dense_rank(F_statistic)) %>%
  arrange(Gene_Name)


