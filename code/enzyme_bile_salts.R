library(dplyr) # data manipulation
library(ggplot2) # plotting

df <- read.csv("./data/AqFl2_bile_salts.csv", header = T, na.strings = c(""))
df

# convert numeric variable to string
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Reorder varialbes in the desired order
df$Diet <- factor(df$Diet, levels = c("REF", "IM33", "IM66", "IM100"))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI1","PI2", "MI", "DI1", "DI2"))

# Exploratory analysis ####################################################################################### 
# Inpsect raw data (distribution and outliers)

df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  group_by(Gut_segment, Diet) %>%
  summarize(N = n(), Mean = mean(Bile_salts_concentration), SD = sd(Bile_salts_concentration), SE = SD / N^(1/2)) %>%
  ggplot(aes(x = Gut_segment, y = Mean, fill = Diet)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SE), 
                size = 0.3, 
                width= 0.3,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0, 250), breaks = 0:5*50, expand = c(0, 0)) +
  theme_bw() +
  theme_classic()+
  labs(title = "Bile salts", x = "Gut segment", y = "µmol/g DM") +
  theme(axis.text.x=element_text(size = 16, face = "bold"),
        axis.text.y=element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.text=element_text(face="bold",size=18),
        legend.title=element_text(face="bold",size=20),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5)) 

ggsave("Bile salts.tiff", units = "in", dpi = 300, compression = "lzw", path = "./results/figures")

# Statistics ##########################################################################################################
# Convert Diet and Target back to character, otherwise the split function won't work correctly
df$Diet <- as.character(df$Diet)
df$Gut_segment <- as.character(df$Gut_segment)

# Split the data frames by "Gut_segment"
df_spl <- df %>% 
  filter(Diet %in% c('REF', 'IM100')) %>% 
  split(f = .$Gut_segment)

# lapply() linear mixed model to each data frame in the list 
t_test <- lapply(df_spl, function(x) t.test(Bile_salts_concentration ~ Diet, data = x))

# Extract p values
p_values <- data.frame(row.names = NULL,
                       "Gut segment" = names(t_test),
                       "p values" = unlist(lapply(t_test, function(x) x[[3]]))
)


