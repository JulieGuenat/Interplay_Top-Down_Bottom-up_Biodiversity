################################################################################
############### CHAPTER I : EFFECT OF TREATMENTS ON DIVERSITY ##################
##################       Julie Morgane Guenat       ############################
#################          December 2024            ############################             
################################################################################

setwd("C:/Users/jguenat/OneDrive - Université de Lausanne/PhD_Analyses/Chapter_1_Effect_of_treat_on_div")


# 1. Packages ####
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(vegan)
library(lme4)
library(car)
library(effectsize)
library(broom.mixed)
library(cowplot)
library(broom.mixed)
library(flextable)
library(officer)


###############################################################################
##################### ALPHA DIVERSITY ANALYSES DNA ############################
###############################################################################

########## Initial Datasets DNA ###########################

data <- read.csv("../Cleaning_data/final_datasets/1.Euka_DNA_Taxo_verified_pres_abs.csv", header = T, sep = ",")
samp <- read.csv("../Cleaning_data/final_datasets/sample_mapping_names.csv", header = T, sep = ",")
samples_pla <- read.csv("../Cleaning_data/sample_names.csv", sep= ";", header=T)
metadata <- read.csv("../Cleaning_data/Samples_info_Planaqua.csv", header= T, sep=";", dec= ",")

########## MOTUs richness ###########################
# 1. Creating the MOTUs richness dataset ####
# first we consider each line as a different species. 
# Meaning that here, we study the MOTUs richness. 

## 1.1 Terrestrial + aquatic MOTUs ####
Motus_rich_all <- data[,c(10:194)] #select the columns corresponding to pres/abs
Motus_rich_all <- as.data.frame(colSums(Motus_rich_all)) # computing the colsum == # of MOTUs 
colnames(Motus_rich_all) <- "Motus_rich_all"
Motus_rich_all$sample_id <- row.names(Motus_rich_all) 

## 1.2 Aquatic MOTUS ####
aquatic <- data %>% filter(aquatic==1)
Motus_rich_aqua <- aquatic[,c(10:194)] 
Motus_rich_aqua <- as.data.frame(colSums(Motus_rich_aqua))
colnames(Motus_rich_aqua) <- "Motus_rich_aqua"
Motus_rich_aqua$sample_id <- row.names(Motus_rich_aqua) 

## 1.3 join MOtus richness ####
a1.Motus_rich<-merge(Motus_rich_all, Motus_rich_aqua, by="sample_id")

rm(Motus_rich_all, Motus_rich_aqua)

########## Phylum richness ###########################
# 2. Creating the Phylum richness dataset ####
phylum<- data %>%
  dplyr::select(., c(2, 8:194)) %>% 
  group_by(phylum) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

## 2.1 Terrestrial + aquatic Phylum ####
Phylum_rich_all <- phylum[,c(4:188)] #select the columns corresponding to pres/abs
Phylum_rich_all <- as.data.frame(colSums(Phylum_rich_all)) # computing the colsum == # of MOTUs 
colnames(Phylum_rich_all) <- "Phylum_rich_all"
Phylum_rich_all$sample_id <- row.names(Phylum_rich_all) 

## 2.2 Aquatic MOTUS ####
aquatic <- phylum %>% filter(aquatic==1)
Phylum_rich_aqua <- aquatic[,c(4:188)] 
Phylum_rich_aqua <- as.data.frame(colSums(Phylum_rich_aqua))
colnames(Phylum_rich_aqua) <- "Phylum_rich_aqua"
Phylum_rich_aqua$sample_id <- row.names(Phylum_rich_aqua) 

## 2.3 join MOtus richness ####
a2.Phylum_rich<-merge(Phylum_rich_all, Phylum_rich_aqua, by="sample_id")

rm(phylum, Phylum_rich_all, Phylum_rich_aqua, aquatic)

########## Class richness ###########################
# 3. Creating the Class richness dataset ####
Class<- data %>%
  dplyr::select(., c(3, 8:194)) %>% 
  group_by(class) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

## 3.1 Terrestrial + aquatic Class ####
Class_rich_all <- Class[,c(4:188)] #select the columns corresponding to pres/abs
Class_rich_all <- as.data.frame(colSums(Class_rich_all)) # computing the colsum == # of MOTUs 
colnames(Class_rich_all) <- "Class_rich_all"
Class_rich_all$sample_id <- row.names(Class_rich_all) 

## 3.2 Aquatic MOTUS ####
aquatic <- Class %>% filter(aquatic==1)
Class_rich_aqua <- aquatic[,c(4:188)] 
Class_rich_aqua <- as.data.frame(colSums(Class_rich_aqua))
colnames(Class_rich_aqua) <- "Class_rich_aqua"
Class_rich_aqua$sample_id <- row.names(Class_rich_aqua) 

## 3.3 join MOtus richness ####
a3.Class_rich<-merge(Class_rich_all, Class_rich_aqua, by="sample_id")

rm(Class, Class_rich_all, Class_rich_aqua, aquatic)

########## Order richness ###########################
# 4. Creating the Order richness dataset ####
Order<- data %>%
  dplyr::select(., c(4, 8:194)) %>% 
  group_by(order) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

## 4.1 Terrestrial + aquatic Order ####
Order_rich_all <- Order[,c(4:188)] #select the columns corresponding to pres/abs
Order_rich_all <- as.data.frame(colSums(Order_rich_all)) # computing the colsum == # of MOTUs 
colnames(Order_rich_all) <- "Order_rich_all"
Order_rich_all$sample_id <- row.names(Order_rich_all) 

## 4.2 Aquatic MOTUS ####
aquatic <- Order %>% filter(aquatic==1)
Order_rich_aqua <- aquatic[,c(4:188)] 
Order_rich_aqua <- as.data.frame(colSums(Order_rich_aqua))
colnames(Order_rich_aqua) <- "Order_rich_aqua"
Order_rich_aqua$sample_id <- row.names(Order_rich_aqua) 

## 4.3 join MOtus richness ####
a4.Order_rich<-merge(Order_rich_all, Order_rich_aqua, by="sample_id")

rm(Order, Order_rich_all, Order_rich_aqua, aquatic)


########## Family richness ###########################
# 5. Creating the Family richness dataset ####
Family<- data %>%
  dplyr::select(., c(5, 8:194)) %>% 
  group_by(family) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))

## 5.1 Terrestrial + aquatic Family ####
Family_rich_all <- Family[,c(4:188)] #select the columns corresponding to pres/abs
Family_rich_all <- as.data.frame(colSums(Family_rich_all)) # computing the colsum == # of MOTUs 
colnames(Family_rich_all) <- "Family_rich_all"
Family_rich_all$sample_id <- row.names(Family_rich_all) 

## 5.2 Aquatic MOTUS ####
aquatic <- Family %>% filter(aquatic==1)
Family_rich_aqua <- aquatic[,c(4:188)] 
Family_rich_aqua <- as.data.frame(colSums(Family_rich_aqua))
colnames(Family_rich_aqua) <- "Family_rich_aqua"
Family_rich_aqua$sample_id <- row.names(Family_rich_aqua) 

## 5.3 join MOtus richness ####
a5.Family_rich<-merge(Family_rich_all, Family_rich_aqua, by="sample_id")

rm(Family, Family_rich_all, Family_rich_aqua, aquatic)

########## Richness Dataset ###########################
# 6.1 Merge all the richness ####
datasets <- list(a1.Motus_rich, a2.Phylum_rich, a3.Class_rich, a4.Order_rich, a5.Family_rich)

# Merge all datasets by "sample_id"
Richness <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), datasets)

# 6.2 add the samples names ####
# Merge simplified samples names with richness
Richness_ID <- merge(samp, Richness, by.x = "Simplified_Name", by.y = "sample_id") %>% 
  group_by(Simplified_Name) %>% slice(1) 

# merge code names metadata
Richness_ID_Plana <- merge(Richness_ID, samples_pla[, c("X", "Samples_names")], 
                           by.x = "Sample_ID", by.y = "X", 
                           all.x = TRUE)

# 6.3 add the metadata ####
Richness_metadata <- merge(metadata, Richness_ID_Plana, by.x = "Samples_ID", by.y = "Samples_names")
#remove temporary datasets
rm(Richness, Richness_ID, Richness_ID_Plana)


#adding treatment as a variable 
Final_richness <- Richness_metadata %>%
  mutate(
    predator = case_when(
      Treatment %in% c("NP", "P") ~ 1,
      Treatment %in% c("none", "N") ~ 0
    ),
    nutrients = case_when(
      Treatment %in% c("NP", "N") ~ 1,
      Treatment %in% c("none", "P") ~ 0
    )
  )

#create the final dataset ####
#write.csv(Final_richness, "../Cleaning_data/final_datasets/3.Euka_DNA_richness_metadata.csv", row.names = F)

########## Raw data Visualization #######################
# data Loading ####
rm(list = ls()) 
df <- read.csv("../Cleaning_data/final_datasets/3.Euka_DNA_richness_metadata.csv", header=T, sep=)

# Load the graphic settings: 
load("../Graph_settings.RData")

# We need to re-level the factors to be able to have the normality: 
# Presence of predator and absence of nutrients.
df$predator <- factor(df$predator, levels = c(1, 0))  # Reference: predator present (1)
df$nutrients <- factor(df$nutrients, levels = c(0, 1))  # Reference: no nutrients (0)

# Create treatment combinations using the predefined labels
df$Treatment_Combined <- paste0(predator_labels[as.character(df$predator)], "/", 
                                nutrient_labels[as.character(df$nutrients)])

# 1.boxplots raw data ####
## 1.1 Reshape data into long format ####
data_long <- df %>%
  pivot_longer(
    cols = c(27:36),
    names_to = c("Taxonomic_Level", "Environment"),
    names_pattern = "(.*)_rich_(.*)",
    values_to = "Richness"
  )

# Clean up levels for better readability
data_long <- data_long %>%
  mutate(
    Taxonomic_Level = case_when(
      grepl("Phylum", Taxonomic_Level) ~ "Phylum",
      grepl("Class", Taxonomic_Level) ~ "Class",
      grepl("Order", Taxonomic_Level) ~ "Order",
      grepl("Family", Taxonomic_Level) ~ "Family",
      grepl("Motus", Taxonomic_Level) ~ "MOTUs",
      TRUE ~ Taxonomic_Level
    ),
    Environment = ifelse(Environment == "all", "Terrestrial and Aquatic", "Aquatic")
  )

# Reorder the taxonomic levels: 
data_long$Taxonomic_Level <- factor(data_long$Taxonomic_Level, 
                                    levels = c("Phylum", "Class", "Order", "Family", "MOTUs"))

# Reorder treatments using the combined treatment labels
data_long$Treatment_Combined <- factor(data_long$Treatment_Combined, 
                                       levels = c("Pred+/Nutr-", "Pred-/Nutr-", 
                                                  "Pred+/Nutr+", "Pred-/Nutr+"))

# Reorder environment 
data_long$Environment <- factor(data_long$Environment, 
                                levels = c("Terrestrial and Aquatic", "Aquatic"))

## 1.2 FIG.S3 Plot ####
Fig.S3 <- ggplot(data_long, aes(x = Treatment_Combined, y = Richness, fill = Environment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_point(color = "black", position = position_dodge(width = 0.8), 
             size = 1.4, alpha = 0.25) +
  facet_wrap(~ Taxonomic_Level, scales = "free_y") + 
  expand_limits(y = 0) +
  labs(y = "Richness (α-diversity)") +
  scale_fill_manual(values = c("Terrestrial and Aquatic" = "#D6BFA7", "Aquatic" = "#6BAAAE")) +
  theme_jg_bottom +  # Use predefined theme
  theme(
    # Override specific elements while keeping the base theme
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
    legend.position = c(0.8, 0.25),  # Custom legend position for this plot
    strip.text = element_text(size = 14), 
    panel.spacing = unit(1.7, "lines"), # Increase panel spacing for clarity
    plot.margin = margin(1, 1, 1, 1, "lines") # Adjust margins to fit the legend
  ) +
  guides(
    fill = guide_legend(title = ""),  # Legend title
    color = "none"  # Remove point legend
  )

Fig.S3

ggsave("figures/S3-Figure_raw_alpha_div.png", 
       plot = Fig.S3,
       width = 12, 
       height = 8, 
       units = "in",
       dpi = 300,
       bg = "white")


########## Richness Analysis ###########################
# 1. GLMM MODEL ####
##########    PHYLUM    ################
## 1.1 Phylum all ####
### 1.1.1 raw data visualization  ####
hist(df$Phylum_rich_all, 20)
table(df$Phylum_rich_all)

### 1.1.2 GLMM poisson ####
m1 <- glmmTMB(Phylum_rich_all ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

m1.2 <-glmmTMB(Phylum_rich_all ~ predator * nutrients + (1|Collection_month) + (1|Block_ID), data = df,
               family = poisson())
anova(m1, m1.2)

# Residual Diagnostics 
simulation_output <- simulateResiduals(m1)
plot(simulation_output)

# Model Summary
summary(m1)

# Type II Anova
car::Anova(m1)

### 1.1.3 Effect interaction graphs ####
effects_m1 <- effects::allEffects(m1)
plot(effects_m1)
interaction_effect1 <- as.data.frame(effects_m1$`predator:nutrients`)

## 1.2 Phylum aquatic ####
### 1.2.1 raw data visualization  ####
hist(df$Phylum_rich_aqua, 20)
table(df$Phylum_rich_aqua)

### 1.2.2 GLMM poisson ####
m2 <- glmmTMB(Phylum_rich_aqua ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

m2.2 <- glmmTMB(Phylum_rich_aqua ~ predator * nutrients + (1|Collection_month) + (1|Block_ID), data = df,
                family = poisson())
anova(m2, m2.2)

# Residual Diagnostics 
simulation_output <- simulateResiduals(m2)
plot(simulation_output)

# Model Summary
summary(m2)

# Type II Anova
car::Anova(m2)

### 1.2.3 Effect interaction graphs ####
effects_m2 <- effects::allEffects(m2)
plot(effects_m2)
interaction_effect2 <- as.data.frame(effects_m2$`predator:nutrients`)

##########     CLASS     ################
## 1.3 Class all ####
### 1.3.1 raw data visualization  ####
hist(df$Class_rich_all, 20)
table(df$Class_rich_all)

### 1.3.2 GLMM poisson ####
m3 <- glmmTMB(Class_rich_all ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m3)
plot(simulation_output)

# Model Summary
summary(m3)

# Type II Anova
car::Anova(m3)

### 1.3.3 Effect interaction graphs ####
effects_m3 <- effects::allEffects(m3)
plot(effects_m3)
interaction_effect3 <- as.data.frame(effects_m3$`predator:nutrients`)

## 1.4 Class aquatic ####
### 1.4.1 raw data visualization  ####
hist(df$Class_rich_aqua, 20)
table(df$Class_rich_aqua)

### 1.4.2 GLMM poisson ####
m4 <- glmmTMB(Class_rich_aqua ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m4)
plot(simulation_output)

# Model Summary
summary(m4)

# Type II Anova
car::Anova(m4)

### 1.4.3 Effect interaction graphs ####
effects_m4 <- effects::allEffects(m4)
plot(effects_m4)
interaction_effect4 <- as.data.frame(effects_m4$`predator:nutrients`)

##########    ORDER     ################
## 1.5 Order all ####
### 1.5.1 raw data visualization  ####
hist(df$Order_rich_all, 20)
table(df$Order_rich_all)

### 1.5.2 GLMM poisson ####
m5 <- glmmTMB(Order_rich_all ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m5)
plot(simulation_output)

# Model Summary
summary(m5)

# Type II Anova
car::Anova(m5)

### 1.5.3 Effect interaction graphs ####
effects_m5 <- effects::allEffects(m5)
plot(effects_m5)
interaction_effect5 <- as.data.frame(effects_m5$`predator:nutrients`)

## 1.6 Order aquatic ####
### 1.6.1 raw data visualization  ####
hist(df$Order_rich_aqua, 20)
table(df$Order_rich_aqua)

### 1.6.2 GLMM poisson ####
m6 <- glmmTMB(Order_rich_aqua ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m6)
plot(simulation_output)

# Model Summary
summary(m6)

# Type II Anova
car::Anova(m6)

### 1.6.3 Effect interaction graphs ####
effects_m6 <- effects::allEffects(m6)
plot(effects_m6)
interaction_effect6 <- as.data.frame(effects_m6$`predator:nutrients`)

##########    FAMILY    ################
## 1.7 Family all ####
### 1.7.1 raw data visualization  ####
hist(df$Family_rich_all, 20)
table(df$Family_rich_all)

### 1.7.2 GLMM poisson ####
m7 <- glmmTMB(Family_rich_all ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m7)
plot(simulation_output)

# Model Summary
summary(m7)

# Type II Anova
car::Anova(m7)

### 1.7.3 Effect interaction graphs ####
effects_m7 <- effects::allEffects(m7)
plot(effects_m7)
interaction_effect7 <- as.data.frame(effects_m7$`predator:nutrients`)

## 1.8 Family aquatic ####
### 1.8.1 raw data visualization  ####
hist(df$Family_rich_aqua, 20)
table(df$Family_rich_aqua)

### 1.8.2 GLMM poisson ####
m8 <- glmmTMB(Family_rich_aqua ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m8)
plot(simulation_output)

# Model Summary
summary(m8)

# Type II Anova
car::Anova(m8)

### 1.8.3 Effect interaction graphs ####
effects_m8 <- effects::allEffects(m8)
plot(effects_m8)
interaction_effect8 <- as.data.frame(effects_m8$`predator:nutrients`)

##########    MOTUs    ################
## 1.9 MOTUs all ####
### 1.9.1 raw data visualization  ####
hist(df$Motus_rich_all, 20)
table(df$Motus_rich_all)

### 1.9.2 GLMM poisson ####
m9 <- glmmTMB(Motus_rich_all ~ predator * nutrients + (1|Collection_month), data = df,
              family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m9)
plot(simulation_output)

# Model Summary
summary(m9)

# Type II Anova
car::Anova(m9)

### 1.9.3 Effect interaction graphs ####
effects_m9 <- effects::allEffects(m9)
plot(effects_m9)
interaction_effect9 <- as.data.frame(effects_m9$`predator:nutrients`)

## 1.10 MOTUs aquatic ####
### 1.10.1 raw data visualization  ####
hist(df$Motus_rich_aqua, 20)
table(df$Motus_rich_aqua)

### 1.10.2 GLMM poisson ####
m10 <- glmmTMB(Motus_rich_aqua ~ predator * nutrients + (1|Collection_month), data = df,
               family = poisson())

# Residual Diagnostics 
simulation_output <- simulateResiduals(m10)
plot(simulation_output)

# Model Summary
summary(m10)

# Type II Anova
car::Anova(m10)

### 1.10.3 Effect interaction graphs ####
effects_m10 <- effects::allEffects(m10)
plot(effects_m10)
interaction_effect10 <- as.data.frame(effects_m10$`predator:nutrients`)

############################################################################
# 2.1 COMBINED PLOTS (Terrestrial + Aquatic) - Figure S4 ####
############################################################################
## Phylum Combined Plot ####
# Prepare combined data for Phylum
interaction_effect1$environment <- "Terrestrial and Aquatic"
interaction_effect2$environment <- "Aquatic"
combined_data_phylum <- rbind(interaction_effect1, interaction_effect2)

# Factor ordering
combined_data_phylum$environment <- factor(combined_data_phylum$environment, 
                                           levels = c("Terrestrial and Aquatic", "Aquatic"))

# Convert labels for predator only (keep nutrients as 0,1 for color mapping)
combined_data_phylum$predator_label <- predator_labels[as.character(combined_data_phylum$predator)]
combined_data_phylum$predator_label <- factor(combined_data_phylum$predator_label, 
                                              levels = c("Pred+", "Pred-"))

e1 <- ggplot(combined_data_phylum, aes(x = predator_label, y = fit, 
                                       color = nutrients, 
                                       group = interaction(nutrients, environment), 
                                       linetype = environment)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "Nutrients") +
  scale_linetype_manual(values = c("Terrestrial and Aquatic" = "solid", 
                                   "Aquatic" = "dashed"), 
                        name = "Environment") + 
  labs(x = "",  
       y = "Predicted α-diversity",  
       title = "Phylum") + 
  theme_jg_right +
  theme(
    legend.box = "vertical"
  )

e1

## Class Combined Plot ####
# Prepare combined data for Class
interaction_effect3$environment <- "Terrestrial and Aquatic"
interaction_effect4$environment <- "Aquatic"
combined_data_class <- rbind(interaction_effect3, interaction_effect4)

# Factor ordering
combined_data_class$environment <- factor(combined_data_class$environment, 
                                          levels = c("Terrestrial and Aquatic", "Aquatic"))

# Convert labels for predator only
combined_data_class$predator_label <- predator_labels[as.character(combined_data_class$predator)]
combined_data_class$predator_label <- factor(combined_data_class$predator_label, 
                                             levels = c("Pred+", "Pred-"))

e2 <- ggplot(combined_data_class, aes(x = predator_label, y = fit, 
                                      color = nutrients, 
                                      group = interaction(nutrients, environment), 
                                      linetype = environment)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "Nutrients") +
  scale_linetype_manual(values = c("Terrestrial and Aquatic" = "solid", 
                                   "Aquatic" = "dashed"), 
                        name = "Environment") + 
  labs(x = "",  
       y = "",  
       title = "Class") + 
  theme_jg_right +
  theme(
    legend.position = "none"
  )

e2

## Order Combined Plot ####
# Prepare combined data for Order
interaction_effect5$environment <- "Terrestrial and Aquatic"
interaction_effect6$environment <- "Aquatic"
combined_data_order <- rbind(interaction_effect5, interaction_effect6)

# Factor ordering
combined_data_order$environment <- factor(combined_data_order$environment, 
                                          levels = c("Terrestrial and Aquatic", "Aquatic"))

# Convert labels for predator only
combined_data_order$predator_label <- predator_labels[as.character(combined_data_order$predator)]
combined_data_order$predator_label <- factor(combined_data_order$predator_label, 
                                             levels = c("Pred+", "Pred-"))

e3 <- ggplot(combined_data_order, aes(x = predator_label, y = fit, 
                                      color = nutrients, 
                                      group = interaction(nutrients, environment), 
                                      linetype = environment)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "Nutrients") +
  scale_linetype_manual(values = c("Terrestrial and Aquatic" = "solid", 
                                   "Aquatic" = "dashed"), 
                        name = "Environment") + 
  labs(x = "",  
       y = "",  
       title = "Order") + 
  theme_jg_right +
  theme(
    legend.position = "none"
  )

e3

## Family Combined Plot ####
# Prepare combined data for Family
interaction_effect7$environment <- "Terrestrial and Aquatic"
interaction_effect8$environment <- "Aquatic"
combined_data_family <- rbind(interaction_effect7, interaction_effect8)

# Factor ordering
combined_data_family$environment <- factor(combined_data_family$environment, 
                                           levels = c("Terrestrial and Aquatic", "Aquatic"))

# Convert labels for predator only
combined_data_family$predator_label <- predator_labels[as.character(combined_data_family$predator)]
combined_data_family$predator_label <- factor(combined_data_family$predator_label, 
                                              levels = c("Pred+", "Pred-"))

e4 <- ggplot(combined_data_family, aes(x = predator_label, y = fit, 
                                       color = nutrients, 
                                       group = interaction(nutrients, environment), 
                                       linetype = environment)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "Nutrients") +
  scale_linetype_manual(values = c("Terrestrial and Aquatic" = "solid", 
                                   "Aquatic" = "dashed"), 
                        name = "Environment") + 
  labs(x = "",  
       y = "Predicted α-diversity",  
       title = "Family") + 
  theme_jg_right +
  theme(
    legend.position = "none"
  )

e4

## MOTUs Combined Plot ####
# Prepare combined data for MOTUs
interaction_effect9$environment <- "Terrestrial and Aquatic"
interaction_effect10$environment <- "Aquatic"
combined_data_motus <- rbind(interaction_effect9, interaction_effect10)

# Factor ordering
combined_data_motus$environment <- factor(combined_data_motus$environment, 
                                          levels = c("Terrestrial and Aquatic", "Aquatic"))

# Convert labels for predator only
combined_data_motus$predator_label <- predator_labels[as.character(combined_data_motus$predator)]
combined_data_motus$predator_label <- factor(combined_data_motus$predator_label, 
                                             levels = c("Pred+", "Pred-"))

e5 <- ggplot(combined_data_motus, aes(x = predator_label, y = fit, 
                                      color = nutrients, 
                                      group = interaction(nutrients, environment), 
                                      linetype = environment)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "Nutrients") +
  scale_linetype_manual(values = c("Terrestrial and Aquatic" = "solid", 
                                   "Aquatic" = "dashed"), 
                        name = "Environment") + 
  labs(x = "",  
       y = "",  
       title = "MOTUs") + 
  theme_jg_right +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

e5

## Combine all plots for Figure S2 ####
# Extract legend from e1
legend <- get_legend(e1)
e1_no_legend <- e1 + theme(legend.position = "none")

# Combine plots
combined_plot_S2 <- plot_grid(e1_no_legend, e2, e3, e4, e5, ncol = 3)

Fig.S4 <- combined_plot_S2 + 
  annotation_custom(
    grob = legend, 
    xmin = -0.45, xmax = 2,
    ymin = 0.65, ymax = 0
  )

Fig.S4

# Export Figure S4
ggsave("figures/S4-Figure_S4_combined_effects.png", 
       plot = Fig.S4,
       width = 15, 
       height = 10, 
       units = "in",
       dpi = 300,
       bg = "white")


############################################################################
# 2.2 AQUATIC ONLY PLOTS - Main Figure 2 ####
############################################################################
## Class Aquatic Plot ####
# Convert labels for Class aquatic (predator only)
interaction_effect4$predator_label <- predator_labels[as.character(interaction_effect4$predator)]
interaction_effect4$predator_label <- factor(interaction_effect4$predator_label, 
                                             levels = c("Pred+", "Pred-"))

Class_aqua <- ggplot(interaction_effect4, aes(x = predator_label, y = fit, 
                                              color = nutrients, group = nutrients)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1.2) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "") +
  labs(title = "Class", 
       x = "",
       y = "Predicted α-diversity") +
  theme_jg_right +
  theme(
    legend.position = c(1.02, 1.2), 
    legend.justification = c("right", "top")
  )

Class_aqua

## Order Aquatic Plot ####
# Convert labels for Order aquatic (predator only)
interaction_effect6$predator_label <- predator_labels[as.character(interaction_effect6$predator)]
interaction_effect6$predator_label <- factor(interaction_effect6$predator_label, 
                                             levels = c("Pred+", "Pred-"))

Order_aqua <- ggplot(interaction_effect6, aes(x = predator_label, y = fit, 
                                              color = nutrients, group = nutrients)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1.2) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "") +
  labs(title = "Order", 
       x = "",
       y = "") +
  theme_jg_right +
  theme(
    legend.position = c(1.02, 1.2), 
    legend.justification = c("right", "top")
  )

Order_aqua

## Family Aquatic Plot ####
# Convert labels for Family aquatic (predator only)
interaction_effect8$predator_label <- predator_labels[as.character(interaction_effect8$predator)]
interaction_effect8$predator_label <- factor(interaction_effect8$predator_label, 
                                             levels = c("Pred+", "Pred-"))

Family_aqua <- ggplot(interaction_effect8, aes(x = predator_label, y = fit, 
                                               color = nutrients, group = nutrients)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1.2) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "") +
  labs(title = "Family", 
       x = "",
       y = "Predicted α-diversity") +
  theme_jg_right +
  theme(
    legend.position = c(1.02, 1.2), 
    legend.justification = c("right", "top")
  )

Family_aqua

## MOTUs Aquatic Plot ####
# Convert labels for MOTUs aquatic
interaction_effect10$predator_label <- predator_labels[as.character(interaction_effect10$predator)]
interaction_effect10$nutrients_label <- nutrient_labels[as.character(interaction_effect10$nutrients)]
interaction_effect10$predator_label <- factor(interaction_effect10$predator_label, 
                                              levels = c("Pred+", "Pred-"))

MOTUs_aqua <- ggplot(interaction_effect10, aes(x = predator_label, y = fit, 
                                               color = nutrients, group = nutrients_label)) +
  geom_line(linewidth = 1.2) +  
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), 
                width = 0.15, linewidth = 1.2) +  
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "") +
  labs(title = "MOTUs", 
       x = "",
       y = "") +
  theme_jg_right +
  theme(
    legend.position = c(1.02, 1.2), 
    legend.justification = c("right", "top")
  )

MOTUs_aqua


## Combine aquatic plots for Figure 2 ####
Fig2_aquatic <- plot_grid(Class_aqua, Order_aqua, Family_aqua, MOTUs_aqua, ncol = 2)
Fig2_aquatic

# Export Figure 2 (aquatic only)
ggsave("figures/M2-Figure_2_aquatic_effects.png", 
       plot = Fig2_aquatic,
       width = 10, 
       height = 10, 
       units = "in",
       dpi = 300,
       bg = "white")

ggsave("figures/M2-Figure_2_aquatic_effects.pdf", 
       plot = Fig2_aquatic,
       width = 10, 
       height = 10, 
       units = "in",
       device = "pdf",
       bg = "white")

ggsave("figures/M2-Figure_2_aquatic_effects.tiff", 
       plot = Fig2_aquatic,
       width = 10, 
       height = 10, 
       units = "in",
       device = "pdf",
       bg = "white")

ggsave("figures/M2-Figure_2_aquatic_effects_F.png", 
       plot = Fig2_aquatic,
       width = 8, 
       height = 8, 
       units = "in",
       dpi = 300,
       bg = "white")

# ###########################################################################
# 3. COEFFICIENT PLOT FIGURE S5 ####
# ###########################################################################

# Extract effect sizes using broom.mixed
model_list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)  

model_summaries <- lapply(model_list, function(model) {
  tidy(model, effects = "fixed")
})

# Combine results
model_names <- c("Phylum_all", "Phylum_aquatic", "Class_all", "Class_aquatic", 
                 "Order_all", "Order_aquatic", "Family_all", "Family_aquatic", 
                 "MOTUs_all", "MOTUs_aquatic")

results_df <- do.call(rbind, Map(cbind, model_summaries, Model = model_names))

# Data preparation using predefined labels
results_df$Taxon_Type <- ifelse(grepl("aquatic", results_df$Model), "Aquatic", "All")
results_df$Taxon_Type <- factor(results_df$Taxon_Type, 
                                levels = c("All", "Aquatic"), 
                                labels = c("Terrestrial and Aquatic", "Aquatic"))

results_df$Taxonomical_Level <- gsub("_.*", "", results_df$Model)
results_df$Taxonomical_Level <- factor(results_df$Taxonomical_Level,
                                       levels = c("Phylum", "Class", "Order", "Family", "MOTUs"))

# Update term labels using predefined labels structure
results_df$term <- factor(results_df$term, 
                          levels = c("(Intercept)", "predator0:nutrients1", "nutrients1", "predator0"),
                          labels = c("Intercept", 
                                     paste0(predator_labels["0"], " × ", nutrient_labels["1"]), 
                                     nutrient_labels["1"], 
                                     predator_labels["0"]))

# Add significance stars
results_df$stars <- ifelse(results_df$p.value < 0.001, "***", 
                           ifelse(results_df$p.value < 0.01, "**", 
                                  ifelse(results_df$p.value < 0.05, "*", "")))

# Colors for taxonomic levels (keeping your original colors as they work well)
custom_colors <- c("Phylum" = "#563B53", 
                   "Class" = "#91B3B8", 
                   "Order" = "#BFC66E", 
                   "Family" = "#BD974C", 
                   "MOTUs" = "#9E6969")

dodge_width <- 0.4


# Complete coefficient plot with thicker HORIZONTAL confidence interval lines
FigS5 <- ggplot(results_df, aes(x = estimate, y = term, 
                                        xmin = estimate - 1.96 * std.error, 
                                        xmax = estimate + 1.96 * std.error,
                                        color = Taxonomical_Level,
                                        shape = Taxon_Type)) +  
  geom_pointrange(position = position_dodge(width = dodge_width), 
                  size = 1.2, fatten = 1.2) +  
  geom_point(position = position_dodge(width = dodge_width), size = 2) +  
  geom_text(aes(label = stars), color = "grey9",
            position = position_dodge(width = dodge_width), 
            hjust = -0.3, vjust = 0.5, size = 5) +  
  facet_wrap(~Taxonomical_Level, scales = "free_x", nrow = 1) +
  scale_color_manual(values = custom_colors, guide = "none") +
  scale_shape_manual(values = c("Terrestrial and Aquatic" = 20, "Aquatic" = 18)) +  
  theme_jg_bottom +  # Changed from theme_minimal()
  theme(
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),  # Kept your original strip text styling
    panel.spacing = unit(2, "lines")
  ) +
  labs(x = "Estimate", shape = "Taxon Type") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray9") +
  coord_cartesian(clip = "off") 

FigS5

# Export Figure 1
ggsave("figures/S5-Figure_S5_coefficients_plot.png", 
       plot = FigS5,
       width = 10, 
       height = 6, 
       units = "in",
       dpi = 300,
       bg = "white")



##############################################
# TABLE FUNCTIONS FOR RICHNESS ANALYSIS ####
################################################

## TABLE S3: GLMM Model Results Summary ####

# Create tidy results for all models
create_glmm_table <- function() {
  
  # List of all models with descriptive names
  model_list <- list(
    "Phylum (Terrestrial + Aquatic)" = m1,
    "Phylum (Aquatic)" = m2,
    "Class (Terrestrial + Aquatic)" = m3,
    "Class (Aquatic)" = m4,
    "Order (Terrestrial + Aquatic)" = m5,
    "Order (Aquatic)" = m6,
    "Family (Terrestrial + Aquatic)" = m7,
    "Family (Aquatic)" = m8,
    "MOTUs (Terrestrial + Aquatic)" = m9,
    "MOTUs (Aquatic)" = m10
  )
  
  # Extract tidy results for each model
  model_results <- map_dfr(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    tidy_result <- tidy(model, effects = "fixed", conf.int = TRUE)
    tidy_result$model <- model_name
    return(tidy_result)
  })
  
  # Clean up and format the results
  model_results <- model_results %>%
    select(model, term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
    rename(
      "Taxonomic Level" = model,
      "Term" = term,
      "Estimate" = estimate,
      "Std. Error" = std.error,
      "z value" = statistic,
      "p-value" = p.value,
      "Lower CI" = conf.low,
      "Upper CI" = conf.high
    )
  
  # Update term names using predefined labels
  model_results$Term <- case_when(
    model_results$Term == "(Intercept)" ~ "Intercept",
    model_results$Term == "predator0" ~ paste0(predator_labels["0"]),
    model_results$Term == "nutrients1" ~ paste0(nutrient_labels["1"]),
    model_results$Term == "predator0:nutrients1" ~ paste0(predator_labels["0"], " × ", nutrient_labels["1"]),
    TRUE ~ model_results$Term
  )
  
  # Add significance indicators
  model_results$Significance <- case_when(
    model_results$`p-value` < 0.001 ~ '***',
    model_results$`p-value` < 0.01 ~ '**',
    model_results$`p-value` < 0.05 ~ '*',
    model_results$`p-value` < 0.1 ~ '.',
    TRUE ~ ''
  )
  
  # Combine p-value with significance
  model_results$FormattedPvalue <- paste0(
    sprintf("%.3f", model_results$`p-value`), 
    " ", 
    model_results$Significance
  )
  
  # Prepare final formatted table
  model_results_formatted <- model_results %>%
    select(-Significance, -`p-value`) %>%
    rename(`p-value` = FormattedPvalue)
  
  # Create flextable
  ft_models <- flextable(model_results_formatted) %>%
    colformat_double(j = c("Estimate", "Std. Error", "z value", "Lower CI", "Upper CI"), 
                     digits = 3) %>%
    bold(i = model_results$`p-value` < 0.05, j = "p-value") %>%
    theme_vanilla() %>%
    add_header_lines("GLMM Model Results: Species Richness Analysis") %>%
    width(j = "Taxonomic Level", width = 2) %>%
    width(j = "Term", width = 1.5) %>%
    autofit()
  
  # Save to Word
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft_models)
  print(doc, target = "table/TableS3_GLMM_richness_results.docx")
  
  # Also return the table for viewing
  return(ft_models)
}

## TABLE S4: ANOVA Results for GLMM Models ####

create_anova_table <- function() {
  
  # List of all models with descriptive names
  model_list <- list(
    "Phylum (Terrestrial + Aquatic)" = m1,
    "Phylum (Aquatic)" = m2,
    "Class (Terrestrial + Aquatic)" = m3,
    "Class (Aquatic)" = m4,
    "Order (Terrestrial + Aquatic)" = m5,
    "Order (Aquatic)" = m6,
    "Family (Terrestrial + Aquatic)" = m7,
    "Family (Aquatic)" = m8,
    "MOTUs (Terrestrial + Aquatic)" = m9,
    "MOTUs (Aquatic)" = m10
  )
  
  # Extract ANOVA results for each model
  anova_results <- map_dfr(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    anova_result <- car::Anova(model)
    anova_df <- as.data.frame(anova_result)
    anova_df$Effect <- rownames(anova_df)
    rownames(anova_df) <- NULL
    anova_df$`Taxonomic Level` <- model_name
    return(anova_df)
  })
  
  # Clean up column names
  names(anova_results)[names(anova_results) == "Chisq"] <- "Chi-Square"
  names(anova_results)[names(anova_results) == "Pr(>Chisq)"] <- "p-value"
  
  # Reorder columns
  anova_results <- anova_results[, c("Taxonomic Level", "Effect", "Chi-Square", "Df", "p-value")]
  
  # Update effect names using predefined labels
  anova_results$Effect <- case_when(
    anova_results$Effect == "predator" ~ "Predator",
    anova_results$Effect == "nutrients" ~ "Nutrients", 
    anova_results$Effect == "predator:nutrients" ~ "Predator × Nutrients",
    TRUE ~ anova_results$Effect
  )
  
  # Add significance indicators
  anova_results$Significance <- case_when(
    anova_results$`p-value` < 0.001 ~ '***',
    anova_results$`p-value` < 0.01 ~ '**',
    anova_results$`p-value` < 0.05 ~ '*',
    anova_results$`p-value` < 0.1 ~ '.',
    TRUE ~ ''
  )
  
  # Combine p-value with significance
  anova_results$FormattedPvalue <- paste0(
    sprintf("%.3f", anova_results$`p-value`), 
    " ", 
    anova_results$Significance
  )
  
  # Prepare final formatted table
  anova_results_formatted <- anova_results %>%
    select(-Significance, -`p-value`) %>%
    rename(`p-value` = FormattedPvalue)
  
  # Create flextable
  ft_anova <- flextable(anova_results_formatted) %>%
    colformat_double(j = c("Chi-Square"), digits = 3) %>%
    colformat_int(j = "Df") %>%
    bold(i = anova_results$`p-value` < 0.05, j = "p-value") %>%
    theme_vanilla() %>%
    add_header_lines("ANOVA Results: Species Richness Analysis") %>%
    width(j = "Taxonomic Level", width = 2) %>%
    width(j = "Effect", width = 1.5) %>%
    autofit()
  
  # Save to Word
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft_anova)
  print(doc, target = "table/TableS4_ANOVA_richness_results.docx")
  
  # Also return the table for viewing
  return(ft_anova)
}

## Execute the functions to create tables ####

# Create the tables
glmm_table <- create_glmm_table()
anova_table <- create_anova_table()

# Display the tables
glmm_table
anova_table
