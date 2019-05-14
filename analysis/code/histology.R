library(tidyverse) # tidy, manipulate and plot data
library(cowplot) # a ggplot2 add-on to provide a publication-ready theme 
library(ggsignif) # # add p values to plot
library(gridExtra) # combine figures
library(ordinal) # run ordered logistic regression
library(broom) # tidy model outputs

# Import and tidy data #########################################################
df <- read_csv("data/raw_data/AqFl2_histology.csv", col_names = T, na = "") 

head(df, n = 20L) 
str(df)

# Tidy data
df <- df %>% 
  gather("Variable", "Rank", 4:ncol(df)) %>% # convert data to the "long" format
  na.omit() %>% # remove rows with "NA"
  mutate(Gut_segment = gsub("\\_.*", "", Variable), # create a Gut_segment column 
         Variable = gsub("*.I\\_", "", Variable)) # remove "PI_/MI_/DI_" from "Variable" column 

# *. repalce anything before; replace anything after .*  
# to replace punctuation characters, add double space (\\) prefix
         
# Exploratory analysis #########################################################
# Convert variables to factors, specifying the desired orders shown on plots
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Net_pen <- factor(df$Net_pen, levels = unique(df$Net_pen))
df$Gut_segment <- factor(df$Gut_segment, levels = c("PI", "MI", "DI"))
df$Variable <- factor(df$Variable, levels = c("hpv", "snv", "smc", "lpc", "mfh"))
df$Rank <- factor(df$Rank, levels = c("Normal", "Mild", "Moderate", "Marked", "Severe"))
                                      
# Make histogram
ggplot(df, aes(Rank, fill = Net_pen)) +
  geom_histogram(#width = 1, # remove the gap between bars 
                 stat="count") +
  facet_grid(Gut_segment + Diet ~ Variable) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Dark2")

ggsave('analysis/exploratory_analysis/histology_histogram.pdf', 
       units = "in", 
       dpi = 300,
       height = 8,
       width = 8)

# Statistics ###################################################################
# Aggregate the levels of ranks to 4 
df <- df %>%
  mutate(Rank_agr = gsub("Normal", "1", Rank), 
         Rank_agr = gsub("Mild", "2", Rank_agr),
         Rank_agr = gsub("Moderate", "3", Rank_agr),
         Rank_agr = gsub("Marked|Severe", "4", Rank_agr))

# Convert aggregated ranks to ordered factors
df$Rank_agr <- ordered(df$Rank_agr, levels = c("1", "2", "3", "4"))

# Run ordinal logistic regression and compare models ---------------------------
olr <- df %>%
  group_by(Gut_segment, Variable) %>% 
  nest() %>% 
  mutate(mod_clmm   = map(data, ~clmm(Rank_agr ~ Diet + (1|Net_pen),
                                     data = .x, Hess = T, nAGQ = 10)), # full model with random effect. 
         mod_clm    = map(data, ~clm(Rank_agr ~ Diet, data = .x)), # sub model without random effect
         smr_clmm   = map(mod_clmm, ~summary(.x)), # model summary
         smr_clm    = map(mod_clm, ~summary(.x)),
         Hess_clmm  = map_dbl(smr_clmm, ~pluck(.x, "condHess")), # extract condition number of Hessian
         Hess_clm   = map_dbl(smr_clm, ~pluck(.x, "cond.H")),
         AIC_clmm   = map_dbl(mod_clmm, ~pluck(glance(.x), "AIC")), # extract AIC
         AIC_clm    = map_dbl(mod_clm, ~pluck(glance(.x), "AIC")),
         ranef_var  = map_dbl(mod_clmm, ~(VarCorr(.x))[[1]][1,1]), # extract variance of random effect
         ICC        = map_dbl(ranef_var, ~(.x)/(.x + pi^2/3)), # calcualte intra-class correlation
         lrt_ranef  = map2(mod_clm, mod_clmm, anova), # likelihood ratio test of random effect. Somehow, mod_clm has to be placed before mod_clmm  
         p_ranef    = map_dbl(lrt_ranef, ~.x[2, "Pr(>Chisq)"]/2) # p value of random effect. Since the likelihood ratio test 
        )                                                        # is on the boundary of the parameter space, a more correct
                                                                 # p-value was obtained by halving the p-value.

# Overview of model comparison
olr %>%
  select(1:2, 8:11, 13, 15)
  
# Fit full model to MI_hpv and MI_smc ------------------------------------------
mod_clmm <- olr %>%
  select(1:3) %>%
  filter(Gut_segment == "MI" & Variable %in% c("hpv", "smc")) %>%
  mutate(mod_clmm      = map(data, ~clmm2(Rank_agr ~ Diet, random = Net_pen,
                                          data = .x, Hess = T, nAGQ = 10)),
         mod_clmm_nom  = map(data, ~clmm2(Rank_agr ~ 1, random = Net_pen, 
                                          data = .x, Hess = T, nAGQ = 10,
                                          nominal = ~Diet)), # nominal effect is not available in clmm 
         mod_clmm_sca  = map(data, ~clmm2(Rank_agr ~ Diet, random = Net_pen, 
                                          data = .x, Hess = T, nAGQ = 10,
                                          scale = ~Diet)), # scale effect is not available in clmm
         lrt_nom       = map2(mod_clmm, mod_clmm_nom, anova), # likelihood ratio test of proportional odds assumption via nominal effect
         lrt_sca       = map2(mod_clmm, mod_clmm_sca, anova), # likelihood ratio test of proportional odds assumption via scale effect
         p_nom         = map_dbl(lrt_nom, ~.x[2, "Pr(Chi)"]),
         p_sca         = map_dbl(lrt_sca, ~.x[2, "Pr(Chi)"]),
         mod_clmm_drop = map(data, ~clmm2(Rank_agr ~ 1, random = Net_pen,
                                          data = .x, Hess = T, nAGQ = 10)), # drop1 works for clmm objects but not clmm2 objects
         lrt_fixef     = map2(mod_clmm, mod_clmm_drop, anova), # likelihood ratio test of the diet effect
         p_diet        = map_dbl(lrt_fixef, ~.x[2, "Pr(Chi)"])
  )

# Check random effect 
# Define a function to plot random effect (see clmm2 tutorial from ordinal pakcage)
plot.ranef <- function(mod) { # mod should be a clmm2 object
  rn <- names(mod[["stDev"]]) 
  nc <- length(mod[["ranef"]]) 
  ci <- mod$ranef + qnorm(0.975) * sqrt(mod$condVar) %o% c(-1, 1)
  ord.re <- order(mod$ranef)
  ci <- ci[order(mod$ranef), ]
  plot(1:nc, mod$ranef[ord.re], axes = F, ylim = range(ci),
       xlab = rn, ylab = paste0(rn, " effect"))
  axis(1, at = 1:nc, labels = ord.re)
  axis(2)
  for(i in 1:nc) segments(i, ci[i,1], i, ci[i, 2])
  abline(h = 0, lty = 2)
} 

plot.ranef(mod_clmm[[1, "mod_clmm"]]) # random effect of MI_hpv
plot.ranef(mod_clmm[[2, "mod_clmm"]]) # random effect of MI_smc

# Fit sub model to other response variables ------------------------------------ 
# Nominal_test/scale_test doesn't work with purrr::map(), lapply() or foor loop. 
# To test the mode assumption, we need to run clm() with or without nominal/scale 
# effect and perform likelihood ratio tests
mod_clm <- olr %>%
  select(1:3, mod_clm) %>%
  filter(Gut_segment %in% c("PI", "DI") |
           Gut_segment == "MI" & !Variable %in% c("hpv", "smc")) %>%
  mutate(mod_clm_nom = map(data, ~clm(Rank_agr ~ 1, nominal = ~Diet, data = .x)), # add nominal effect 
         mod_clm_sca = map(data, ~clm(Rank_agr ~ Diet, scale = ~Diet, data = .x, # add scale effect 
                                      control = list(maxIter = 300L))), # increase maximum iterations as some models failed to converge
         lrt_nom     = map2(mod_clm, mod_clm_nom, anova), # likelihood ratio test of proportional odds assumption via nominal effect
         lrt_sca     = map2(mod_clm, mod_clm_sca, anova), # likelihood ratio test of proportional odds assumption via scale effect
         p_nom       = map_dbl(lrt_nom, ~.x[2, "Pr(>Chisq)"]),
         p_sca       = map_dbl(lrt_sca, ~.x[2, "Pr(>Chisq)"])
  ) 

# The warning message that "(1)Hessian is numerically singular: parameters are 
# not uniquely determined. In addition: Absolute convergence criterion was met, 
# but relative criterion was not met" is an indication that the the ML estimates
# of some of the threshold parameters are at infinity. Nontheless, the value of 
# the log-likelihood is accurately determined. Hence, we can trust results from
# the likelihood ratio tests 

# Note that p value of nominal test for PI_lpc is missing. Let's look at the 
# correspondent likelihood ratio test
mod_clm[[2, "lrt_nom"]]

# Turns out the two models are identical. The p value should be 1. Let's replace
# NA with 1
mod_clm <- mod_clm %>%
  mutate(p_nom = replace_na(p_nom, 1)) 

# Next, let's get clm modles that meet the proportional odds assumption 
mod_clm_po <- mod_clm %>%
  filter(p_nom >= 0.05 & p_sca >= 0.05) %>%
  mutate(lrt_diet = map(mod_clm, ~drop1(.x, test = "Chi")), # likelihood ratio test of diet effect
         p_diet   = map_dbl(lrt_diet, ~.x["Diet", "Pr(>Chi)"]))  
         
# And, clm models that breach the proportional odds assumption 
mod_clm_npo <- mod_clm %>%
  filter(p_nom < 0.05 | p_sca < 0.05) %>%
  mutate(lrt_diet = map(mod_clm_sca, ~drop1(.x, test = "Chi")), # likelihood ratio test of diet effect
         p_diet   = map_dbl(lrt_diet, ~.x["Diet", "Pr(>Chi)"])) 

# Finally, get p values of diet effect ----------------------------------------- 
p_diet <- rbind(mod_clmm[ , c("Gut_segment", "Variable", "p_diet")],  
                mod_clm_po[ , c("Gut_segment", "Variable", "p_diet")],
                mod_clm_npo[ , c("Gut_segment", "Variable", "p_diet")])

# Apply multiple comparison correction for each gut segment
p_diet <- p_diet %>%
  group_by(Gut_segment) %>% # correct p values for each gut segment
  nest() %>%
  mutate(p_diet_adj = map(data, ~p.adjust(.x$p_diet, method = "holm"))) %>% 
  unnest() %>% # somehow, unnest change the column position of p_diet_adj
  mutate(Gut_segment = as.character(Gut_segment)) %>% # Convert to character so that it can be sorted
  arrange(desc(Gut_segment)) %>%
  select(Gut_segment, Variable, p_diet, p_diet_adj) # reorder columns

# Make Figure 2 ################################################################
# Calculate new variable "Percent"for plotting stacked bar plot
df_bar <- df %>% 
  group_by(Gut_segment, Diet, Variable, Rank) %>% 
  mutate(observation = n()) %>% # number of fish assigned with different ranks
  ungroup() %>% 
  group_by(Gut_segment, Diet, Variable) %>% 
  mutate(Percent = 100 * observation/n()) %>% 
  ungroup() %>% 
  group_by(Gut_segment, Diet, Variable, Rank, Percent) %>% 
  summarize()

# Define the color scheme 
my_col = c(Normal = "royalblue2", Mild = "peachpuff1", Moderate = "tan1", 
           Marked = "tomato", Severe = "red3")

# Change fecet labels 
levels(df_bar$Variable) <- c("Hypervacuolization", 
                             "Supranuclear vacuolization", 
                             "Submucosal cellularity",
                             "Lamina propria cellularity",
                             "Mucosal fold height")

# Make stacked barplot for each gut segment
p1 <- df_bar %>%
  filter(Gut_segment == "PI") %>%
  ggplot(aes(Diet, Percent)) + 
  geom_bar(aes(fill = forcats::fct_rev(Rank)), # fct_rev() reverses stacked bars
           stat = "identity") + 
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col) + 
  scale_y_continuous(limits = c(0, 105), 
                     breaks = 0:5*20, 
                     expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "PI", y = "%") +
  guides(fill = guide_legend(title = "Rank")) + 
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0), 
        legend.position = "none")

p2 <- df_bar %>%
  filter(Gut_segment == "MI") %>%
  ggplot(aes(Diet, Percent)) + 
  geom_bar(aes(fill = forcats::fct_rev(Rank)), stat = "identity") +
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col, 
                    drop = F) + # forces legend to show all categories
  scale_y_continuous(limits = c(0, 105), 
                     breaks = 0:5*20, 
                     expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "MI", y = "%") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "none")

p3 <- df_bar %>%
  filter(Gut_segment == "DI") %>%
  ggplot(aes(Diet, Percent)) + 
  geom_bar(aes(fill = forcats::fct_rev(Rank)), stat = "identity") +
  facet_wrap(~ Variable, nrow = 1) +
  scale_fill_manual(values = my_col) + 
  scale_y_continuous(limits = c(0, 105), 
                     breaks = 0:5*20, 
                     expand = expand_scale(mult = c(0, 0.05))) +
  labs(title = "DI", y = "%") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0), 
        legend.position = "none")

# Add p values to the plots ----------------------------------------------------
# Make p vlaue annotations
ann <- p_diet %>%
  mutate(Variable = gsub("hpv", "Hypervacuolization", Variable), 
         Variable = gsub("lpc", "Lamina propria cellularity", Variable),
         Variable = gsub("mfh", "Mucosal fold height", Variable),
         Variable = gsub("smc", "Submucosal cellularity", Variable),
         Variable = gsub("snv", "Supranuclear vacuolization", Variable),
         p_diet_adj = formatC(p_diet_adj, format = "f", digits = 3), # format digits of p values
         p_diet_adj = paste("p = ", p_diet_adj), # add label "p =" to p values
         start = rep("REF", nrow(.)), # start position of p value label on x axis
         end = rep("IM", nrow(.)), # end position of p value label on x axis
         y = rep(103, nrow(.)) # position of p value label on y axis
  ) 

# Add p values to the plots. The warning about the missing aesthetics can be ignored
p1 <- p1 + geom_signif(data = filter(ann, Gut_segment == "PI"),
                       aes(xmin = start, 
                           xmax = end, 
                           annotations = p_diet_adj, 
                           y_position = y),
                       textsize = 4, 
                       tip_length = 0,
                       manual = T)

p2 <- p2 + geom_signif(data = filter(ann, Gut_segment == "MI"),
                       aes(xmin = start, 
                           xmax = end, 
                           annotations = p_diet_adj, 
                           y_position = y),
                       textsize = 4, 
                       tip_length = 0,
                       manual = T)

p3 <- p3 + geom_signif(data = filter(ann, Gut_segment == "DI"),
                       aes(xmin = start, 
                           xmax = end, 
                           annotations = p_diet_adj, 
                           y_position = y),
                       textsize = 4, 
                       tip_length = 0,
                       manual = T)
# Combine figures --------------------------------------------------------------
# Extract legend from one of the figures and use it as the shared legend
legend <- get_legend(p1 + theme(legend.position = "right"))

# Make a list of grobs 
gl <- list("1" = p1, "2" = p2, "3" = p3, "4" = legend)  

# Make a layout matrix to guide the layout of figures
layout <- cbind(c(1, 2, 3),
                c(1, 2, 3), 
                c(1, 2, 3), 
                c("NA", 4, 3))

# Use gridExtra::grid.arrange to combine the figures
tiff('results/figures/Figure 2.tiff', 
     compression = "lzw",
     units = "in", 
     res = 300, 
     height = 14,
     width = 12)

grid.arrange(grobs = gl, layout_matrix = layout)

dev.off()

# Get session info
writeLines(capture.output(sessionInfo()), "analysis/code/histology_sessionInfo.txt")
