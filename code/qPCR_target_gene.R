library(dplyr) # data manipulation
library(ggplot2) # plotting
library(ggrepel) # labelling
library(reshape2) # convert data format (long/wide)
library(NMF) # makeing heatmap
library(nlme) ## running linear mixed model
library(car) # tidy statistical outputs

#######################################
#                                     #
# Analysis of gene expression in PI   #
#                                     #
#######################################

# Start the session in greek so that greek letters can be recognised
Sys.setlocale(category = "LC_ALL", locale = "Greek")

## Import data. The encoding of the csv file need to be UTF-8, otherwise greek letters won't be recognized. To do that:
# 1. Save the excel sheet as a Unicode text.
# 2. Open your saved file in Notepad and replace all tab characters with commas (",").
#Select a tab character (select and copy the space between two column headers)
#Open the "Find and Replace" window and replace all tab characters with comma .
# 3. Click Save As, name the file, and change the Encoding: to UTF-8.
# 4. Change the file extension from "*.txt" to "*.csv".

df_pi <- read.csv("C:/Users/ljt89/OneDrive - Norwegian University of Life Sciences/Code demo/AqFl2_qPCR_PI_Target.csv", 
                  header = T, na.strings=c(""), encoding = "UTF-8")

df_pi <- read.csv("Z:/PhD/Lab/qPCR/AqFl2/R/AqFl2_qPCR_PI_Target.csv", header =T, na.strings=c(""), encoding = "UTF-8")

head(df_pi)

# Tidying data and calculate new variables ##########################################################################################################
df_pi1 <- df_pi %>% 
  rename(Sample_ID = X.U.FEFF.Sample_ID) %>% # remove "X.U.FEFF." from the column name, which was reserved to tell that the file is in UTF-8 format.
  select(-Comments) %>% # remove the "Comments" column
  filter(!Target %in% c("il6", "tnfα")) %>% # remove il6 and tnfα from analysis as their expression levels are too low to be reliably quantified
  na.omit() %>% # remove rows with NA value
  mutate(NE_gapdh = `^`(Target_PE,  -Target_Cq) / `^`(gapdh_PE, -gapdh_Cq), # expression normalzied to gapdh
         NE_ranpo2 = `^`(Target_PE,  -Target_Cq) / `^`(rnapo2_PE, -rnapo2_Cq), # expression normalzied to rnapo2
         NE_hprt1 = `^`(Target_PE,  -Target_Cq) / `^`(hprt1_PE, -hprt1_Cq), # expression normalzied to hprt1
         NE_actb = `^`(Target_PE,  -Target_Cq) / `^`(actb_PE, -actb_Cq), # expression normalzied to actb
         MNE = (NE_gapdh * NE_ranpo2 * NE_hprt1 * NE_actb)^(1/4)) # geometric mean of normalzied expression

head(df_pi1)

# Plotting ##########################################################################################################################################
# Before plotting, first, convert numeric variable to string
df_pi1 <- within(df_pi1, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Then, convert "Diet" and "Target" to factor variables with the desired order
df_pi1$Diet <- factor(df_pi1$Diet, levels = c("REF", "IM100"))

df_pi1$Target <- factor(df_pi1$Target, levels = c("myd88", "il1β", "tnfα", "il6", "il8", "cd3γδ", "cd8β", "mhcI", "ifnγ",
                                                  "il17a", "foxp3", "il10", "tgfβ1", "il4", "plin2", "mtp", "apoa1", "apoa4", 
                                                  "apob", "chk", "pcyt1a", "fabp2b", "aqp8ab", "cldn15", "cldn25b", "zo1",
                                                  "cdh1", "muc2",  "cyp1a1", "mta", "hsp70",  "sod1", "cat", "casp6", "pcna", 
                                                  "mmp13")) 

# Make boxplot --------------------------------------------------------------------------------------------------------------------------------------
# First, define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}
# Make boxplot for genes under different functional categories
box_imu_pi <- df_pi1 %>%
  filter(Target_func == "Immune_modulation") %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group
  ggplot(aes(x=Diet, y=MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen),size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, ncol = 4,scales="free_y") +
  geom_label_repel() +
  theme_bw() +
  ggtitle("Gene expression profile in PI: Immune modulation") +
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

box_lpi_pi <- df_pi1 %>%
  filter(Target_func == "lipid_metabolism") %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group
  ggplot(aes(x=Diet, y=MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen),size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, ncol = 4,scales="free_y") +
  geom_label_repel() +
  theme_bw() +
  ggtitle("Gene expression profile in PI: lipid metabolism") +
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

box_brr_pi <- df_pi1 %>%
  filter(Target_func == "barrier_function") %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group
  ggplot(aes(x=Diet, y=MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen),size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, ncol = 3,scales="free_y") +
  geom_label_repel() +
  theme_bw() +
  ggtitle("Gene expression profile in PI: barrier function") +
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

box_xno_pi <- df_pi1 %>%
  filter(Target_func == "xenobiotic_matebolism") %>%
  group_by(Target, Diet) %>%
  mutate(outlier = is_outlier(MNE)) %>% # mark outliers within each diet group
  ggplot(aes(x=Diet, y=MNE,label = ifelse(outlier, Sample_ID, NA))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Net_pen),size = 2, shape = 21, position = position_jitterdodge(0.2)) +
  facet_wrap(~ Target, ncol = 4,scales="free_y") +
  geom_label_repel() +
  theme_bw() +
  ggtitle("Gene expression profile in PI: xenobiotic matebolism") +
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

# Make heatmap --------------------------------------------------------------------------------------------------------------------------------------
# Convert the data format from "long" to "wide"
df_pi_wide <- dcast(df_pi1, Sample_ID + Diet + Net_pen ~ Target, value.var="MNE") %>%
  arrange(Diet, Net_pen) 

# Make a numeric matrix for heatmap plotting
nmp  <- select(df_pi_wide, -Sample_ID, -Diet, -Net_pen) %>%
  data.matrix()
# Assign Sample_ID as row name and transpose the matrix
rownames(nmp) <- df_pi_wide[,1]
nmp_t <- t(nmp)
# Make column and row annoation 
annCol <- select(df_pi_wide, Diet, Net_pen)

x <- c("Immune_modulation", "lipid_metabolism", "barrier_function", "xenobiotic_matebolism")
times <- c(12, 8, 6, 8)
func <- rep(x, times = times)
annRow <- data.frame(function_category = func)
annRow$function_category <- factor(annRow$function_category, levels = c("Immune_modulation",
                                                                        "lipid_metabolism",
                                                                        "barrier_function",
                                                                        "xenobiotic_matebolism"))
# Plotting
tiff('heatmap_PI.tiff', units="in", width=16, height=9, res=300)
aheatmap(nmp_t, 
         border_color = "grey60",
         #filename = "p_heatmap.tiff",
         Rowv = NA,
         Colv = NA,
         fontsize = 14,
         cellwidth = 22,
         cellheight = 16,
         scale = "row",  
         annCol = annCol,
         annRow = annRow)             
dev.off()

# Make Barplot showing fold change of gene expression -----------------------------------------------------------------------------------------------
df_pi2 <- df_pi1 %>%
  filter(Diet == "REF") %>%
  group_by(Target) %>%
  mutate(MNE_group_mean_REF = mean(MNE)) %>%
  ungroup() %>%
  select(Target_func, Target, MNE_group_mean_REF) %>%
  mutate(rowName = row.names(.)) 

df_pi3 <- df_pi1 %>%
  filter(Diet == "IM100") %>%
  rename(MNE_IM100 = MNE) %>%
  mutate(rowName = row.names(.)) %>%
  select(rowName, MNE_IM100) %>%
  full_join(df_pi2, ., by = "rowName") %>%
  mutate(Fold_change = MNE_IM100 / MNE_group_mean_REF) %>%
  group_by(Target_func, Target) %>%
  summarize(N = n(), Mean_fold_change = mean(Fold_change), SD = sd(Fold_change), SE = SD / N^(1/2))

df_pi3$Target_func <- factor(df_pi3$Target_func, levels = c("xenobiotic_matebolism",
                                                            "barrier_function",
                                                            "lipid_metabolism",
                                                            "Immune_modulation"))

bar_pi <- ggplot(df_pi3, aes(x = Target, y = Mean_fold_change, fill = Target_func)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  coord_flip() +
  geom_errorbar(aes(ymin = Mean_fold_change, ymax = Mean_fold_change + SE), 
                size = 0.3, 
                width= 0.3) +
  scale_y_continuous(limits = c(0, 2.25), breaks = 0:10*0.25, expand = c(0, 0)) +
  labs(title = "Gene expression profile in PI", y = "fold change", tag = "A") +
  geom_hline(yintercept = 1, color = "blue", linetype = 2) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, face = "bold"),
        axis.text.y=element_text(size = 12, face = "bold.italic"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank())


# Run linear mixed model and tidy model output #####################################################################################################
# Convert Diet and Target back to character, otherwise the split function won't work correctly
df_pi1$Diet <- as.character(df_pi1$Diet)
df_pi1$Target <- as.character(df_pi1$Target)

# Split the data frames by "Target"
df_pi1_spl <- split(df_pi1, f = df_pi1$Target)

# `lapply()` linear mixed model to each data frame in the list and store results as a list
lme_pi <- lapply(df_pi1_spl, function(x) lme(MNE ~ Diet, random = ~1|Net_pen, data = x))

# Model summary and extract p values of the Diet effect.
summary_pi <- lapply(lme_pi, function(x) summary(x))

p_values_pi <-
  
  
  # Model diagnostics ############################################################################################################
# Extract residuals
res_pi <- lapply(lme_pi, function(x) resid(x))

# 1.Linearity and homoskedasticity of residuals
pdf("resid_plot_PI.pdf", encoding = "CP1253") # turn on pdf graphics device
lapply(
  seq_along(res_pi), # add index to the elements in the list
  function(x) 
  {
    plot(res_pi[[x]], main = names(res_pi)[x], ylab = "residuals")
    abline(0,0)
  }
)
dev.off() # turn off the device

pdf("resid_plot_std_PI.pdf")
lapply(lme_pi, function(x) plot(x))
dev.off()

# 2.Absence of collinearity. When more than one fixed effects are included, the collinearity shoulb be checked.
# Not applicable to the present data, which have only one fixed effect

# 3.Normality of residuals

pdf("qqplot_PI.pdf", encoding = "CP1253")
lapply(
  seq_along(res_pi), 
  function(x) 
  {
    qqnorm(res_pi[[x]], main = paste("Normal Q-Q Plot: ", names(res_pi)[x]))
    qqline(res_pi[[x]])
  }
)
dev.off()

# 4.Absence of influential data points


#######################################
#                                     #
# Analysis of gene expression in DI   #
#                                     #
#######################################


# combine barplot to make Fig.2

