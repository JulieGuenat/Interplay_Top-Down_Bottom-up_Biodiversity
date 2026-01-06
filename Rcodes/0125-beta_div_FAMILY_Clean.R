###############################################################################
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
library(emmeans)
library(patchwork)
library(kableExtra)
library(knitr)
library(gridExtra)
library(grid)
library(ggrepel)
library(flextable)  # For creating formatted tables
library(officer)


# 0. Initial Datasets DNA ####

data <- read.csv("../Cleaning_data/final_datasets/1.Euka_DNA_Taxo_verified_pres_abs.csv", header = T, sep = ",")
samp <- read.csv("../Cleaning_data/final_datasets/sample_mapping_names.csv", header = T, sep = ",")
samples_pla <- read.csv("../Cleaning_data/sample_names.csv", sep= ";", header=T)
metadata <- read.csv("../Cleaning_data/Samples_info_Planaqua.csv", header= T, sep=";", dec= ",")


# 1. Data formatting ####
## 1.1 Selecting aquatic family only ####
aquatic <- data %>% 
  dplyr::filter(aquatic == 1)

family_all <- aquatic %>%
  dplyr::select(c(5, 8:ncol(.))) %>%  # Select columns dynamically
  group_by(family) %>% 
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%  # Summarize numeric columns
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))  # Apply mutation to numeric columns

family_all<- as.data.frame(family_all)
family_all<- family_all[-278,] #rm NA
row.names(family_all)<-family_all$family
family_all<- family_all[,-(1:3)] #rm family_name, aquatic, terrestrial
family_all<- t(family_all)
family_all<-as.data.frame(family_all)

family_all <- family_all %>%
  rownames_to_column("sample_ID") 

## 1.2 adding the metadata ####
# 1st was add the correct names
samp_unique <- samp[!duplicated(samp$Simplified_Name), ]
family_all <- merge(samp_unique, family_all, by.x = "Simplified_Name", by.y = "sample_ID")


family_all <- merge(family_all, samples_pla[, c("X", "Samples_names")], 
                    by.x = "Sample_ID", by.y = "X", 
                    all.x = TRUE)

family_metadat <- merge(metadata, family_all, by.x = "Samples_ID", by.y = "Samples_names")

## 1.3 selecting monthly data only ####
family_metadat <- family_metadat %>%
  filter(Sampling_type == "Monthly")

## 1.4 Column predators and nutrient ####
family_metadat_prednutr <- family_metadat %>%
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

family_metadat_prednutr$predator <- factor(family_metadat_prednutr$predator, levels = c(1, 0))  # Reference: predator present (1)
family_metadat_prednutr$nutrients <- factor(family_metadat_prednutr$nutrients, levels = c(0, 1))  # Reference: no nutrients (0)

rm(aquatic, data, family_all, family_metadat, metadata, samp, samples_pla, samp_unique)


# 2. Beta-div analysis CCA ####
## 2.1 Extract family data####
family_data <- family_metadat_prednutr[, 27:(ncol(family_metadat_prednutr)-2)]

### 2.1.1 removing the Perches ####
family_data <- family_data %>%
  dplyr::select(!Percidae)

### 2.1.2 Matrix families ####
# Make sure to convert it to a matrix
family_matrix <- as.matrix(family_data)
dim(family_matrix)

## 2.2 Environmental matrix ####
env_data <- family_metadat_prednutr[, c(3, 304:305)]
dim(env_data)

## 2.3 CCA Model with all the families ####
cca<- vegan::cca(family_matrix ~ predator * nutrients + Condition (Collection_month), data=env_data)
summary(cca)

set.seed(1)
anova_result <- anova(cca, by = "terms", permutations = 999)
print(anova_result)

## 2.4 Biplot with all the families ####
# Extract the scores for species and environmental variables
species_scores <- scores(cca, display = "species")
env_scores <- scores(cca, display = "bp")

# Extract variance explained by CCA1 and CCA2
var_explained <- c(0.41593*100, 0.33741*100)

# Convert scores to data frames for ggplot
species_df <- as.data.frame(species_scores)
env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_df)

#Assign specific colors to each environmental variable
env_df$Color <- c("#C14863", "#196971", "#BD974C")  # Predators, Nutrients, Interaction (chosen #00A36C for interaction)

## 2.5 Biplot without annotation with all the families ####
# Create the biplot
biplot <- ggplot() +
  # Plot species points
  geom_point(data = species_df, aes(x = CCA1, y = CCA2), color = "black", alpha=0.5) +
  # Plot environmental arrows with custom colors
  geom_segment(data = env_df, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = Color),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.2) +
  # Add arrow labels with custom colors
  geom_text(data = env_df, aes(x = CCA1, y = CCA2, label = Variable, color = Color),
            vjust = -0.9, hjust = 0.5) +
  # Customize the plot
  theme_minimal() +
  coord_fixed() +
  # Add 0,0 cross line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey11") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey11") +
  # Add axis labels with variance explained
  labs(x = paste0("CCA1 (", round(var_explained[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained[2], 1), "%)"),
       color = "Environmental Variables") +
  scale_color_identity() + 
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 16))


print(biplot)

## 2.6 CCA Model removing sporadic families ####
### 2.6.1 investigation of the sporadic families ####
# Calculate overall presence distribution 
# Extract family columns (from column 27 to second-to-last column)
family_cols <- names(family_metadat_prednutr)[27:(ncol(family_metadat_prednutr)-2)]

# Calculate presence percentage for each family
presence_df <- data.frame(
  family = family_cols,
  presence = sapply(family_metadat_prednutr[family_cols], function(x) {
    sum(x > 0) / nrow(family_metadat_prednutr) * 100
  })
)

# Create presence distribution plot
presence_plot <- ggplot(presence_df, aes(x = presence)) +
  geom_histogram(binwidth = 1, fill = "#196971", color = "white", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +  # Add breaks every 10%
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  ) +
  labs(
    title = "Distribution of Family Presence in Samples",
    x = "Presence (%)",
    y = "Number of Families"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Calculate monthly presence
monthly_presence <- family_metadat_prednutr %>%
  group_by(Collection_month) %>%
  summarise(across(all_of(family_cols), ~sum(. > 0))) %>%
  gather(family, present, -Collection_month) %>%
  group_by(Collection_month) %>%
  summarise(families = sum(present > 0))

# Create monthly presence plot
monthly_plot <- ggplot(monthly_presence, aes(x = Collection_month, y = families)) +
  geom_bar(stat = "identity", fill = "#C14863") +
  theme_minimal() +
  labs(
    title = "Number of Families Present by Month",
    x = "Month",
    y = "Number of Families"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Combine plots using patchwork
combined_plot <- presence_plot / monthly_plot +
  plot_layout(heights = c(1, 1))

# Print the combined plot
print(combined_plot)

thresholds <- c(1, 5, 10, 20)
for(thresh in thresholds) {
  rare_families <- sum(presence_df$presence <= thresh)
  cat(sprintf("\nFamilies present in %d%% or fewer samples: %d", 
              thresh, rare_families))
}

#Families present in 1% or fewer samples: 42
#Families present in 5% or fewer samples: 108
#Families present in 10% or fewer samples: 153
#Families present in 20% or fewer samples: 190

### 2.6.2 CCA with 0% sp removed ####
# Remove families with zero presence
nonzero_cols <- colSums(family_matrix > 0) > 0
family_matrix_0 <- family_matrix[, nonzero_cols]

# Run CCA
cca_result_0 <- vegan::cca(family_matrix_0 ~ predator * nutrients + Condition(Collection_month), 
                         data = env_data)

summary(cca_result_0)

# Run permutation test
set.seed(1)
anova_result_0 <- anova(cca_result_0, by = "terms", permutations = 999)

print(anova_result_0)

# Biplot 
# Extract scores
species_scores_0 <- scores(cca_result_0, display = "species")
env_scores_0 <- scores(cca_result_0, display = "bp")

# Convert scores to data frames for ggplot
species_df_0 <- as.data.frame(species_scores_0)
env_df_0 <- as.data.frame(env_scores_0)
env_df_0$Variable <- rownames(env_df_0)
env_df_0$Variable <- as.factor(env_df_0$Variable)
levels(env_df_0$Variable) <- c("Nutrient add.", "Predator rm", "Interaction")

#Assign specific colors to each environmental variable
env_df_0$Color <- c("#C14863", "#196971", "#BD974C")

var_explained_0 <- c(0.41593*100, 0.33741*100)

# Create the biplot
biplot_0 <- ggplot() +
  # Plot species points
  geom_point(data = species_df_0, aes(x = CCA1, y = CCA2), color = "black", alpha=0.5) +
  # Plot environmental arrows with custom colors
  geom_segment(data = env_df_0, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = Color),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.2) +
  # Add arrow labels with custom colors
  geom_text(data = env_df_0, aes(x = CCA1, y = CCA2, label = Variable, color = Color),
            vjust = -0.9, hjust = 0.5) +
  ggtitle("CCA removing family present in 0% of samples")+
  # Customize the plot
  theme_minimal() +
  coord_fixed() +
  # Add 0,0 cross line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey11") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey11") +
  # Add axis labels with variance explained
  labs(x = paste0("CCA1 (", round(var_explained_0[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_0[2], 1), "%)"),
       color = "Environmental Variables") +
  scale_color_identity() + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 11))


print(biplot_0)


### 2.6.3 Removing 1% ####
# Calculate presence percentage for each family
family_presence <- colSums(family_matrix > 0) / nrow(family_matrix) * 100

# Keep families present in more than 5% of samples
keep_families <- names(family_presence[family_presence > 1])
family_matrix_1 <- family_matrix[, keep_families]

# Run CCA
cca_result_1 <- vegan::cca(family_matrix_1 ~ predator * nutrients + Condition(Collection_month), 
                           data = env_data)

summary(cca_result_1)

# Run permutation test
set.seed(1)
anova_result_1 <- anova(cca_result_1, by = "terms", permutations = 999)

print(anova_result_1)

# Biplot 
# Extract scores
species_scores_1 <- scores(cca_result_1, display = "species")
env_scores_1 <- scores(cca_result_1, display = "bp")

# Convert scores to data frames for ggplot
species_df_1 <- as.data.frame(species_scores_1)
env_df_1 <- as.data.frame(env_scores_1)
env_df_1$Variable <- rownames(env_df_1)
env_df_1$Variable <- as.factor(env_df_1$Variable)
levels(env_df_1$Variable) <- c("Nutrient add.", "Predator rm", "Interaction")

#Assign specific colors to each environmental variable
env_df_1$Color <- c("#C14863", "#196971", "#BD974C")

var_explained_1 <- c(0.4247*100, 0.32886*100)

# Create the biplot
biplot_1 <- ggplot() +
  # Plot species points
  geom_point(data = species_df_1, aes(x = CCA1, y = CCA2), color = "black", alpha=0.5) +
  # Plot environmental arrows with custom colors
  geom_segment(data = env_df_1, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = Color),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.2) +
  # Add arrow labels with custom colors
  geom_text(data = env_df_1, aes(x = CCA1, y = CCA2, label = Variable, color = Color),
            vjust = -0.9, hjust = 0.5) +
  ggtitle("CCA removing family present in 1% of samples")+
  # Customize the plot
  theme_minimal() +
  coord_fixed() +
  # Add 0,0 cross line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey11") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey11") +
  # Add axis labels with variance explained
  labs(x = paste0("CCA1 (", round(var_explained_1[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_1[2], 1), "%)"),
       color = "Environmental Variables") +
  scale_color_identity() + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 11))


print(biplot_1)

### 2.6.4 Removing 5% ####
# Conservative threshold (5%): Remove 108 families that appear in ≤5% of samples
# a. It maintains enough families for meaningful community analysis
# b. Removes the most spurious occurrences that might add noise to the analysis
# c. Preserves seasonal patterns while removing truly rare families

# Calculate presence percentage for each family
family_presence <- colSums(family_matrix > 0) / nrow(family_matrix) * 100

# Keep families present in more than 5% of samples
keep_families <- names(family_presence[family_presence > 5])
family_matrix_5 <- family_matrix[, keep_families]

# Run CCA
cca_result_5 <- vegan::cca(family_matrix_5 ~ predator * nutrients + Condition(Collection_month), 
                           data = env_data)

summary(cca_result_5)

# Run permutation test
set.seed(1)
anova_result_5 <- anova(cca_result_5, by = "terms", permutations = 999)

print(anova_result_5)

# Biplot 
# Extract scores
species_scores_5 <- scores(cca_result_5, display = "species")
env_scores_5 <- scores(cca_result_5, display = "bp")

# Convert scores to data frames for ggplot
species_df_5 <- as.data.frame(species_scores_5)
env_df_5 <- as.data.frame(env_scores_5)
env_df_5$Variable <- rownames(env_df_5)
env_df_5$Variable <- as.factor(env_df_5$Variable)
levels(env_df_5$Variable) <- c("Nutrient Enrichment", "Predator Removal", "Interaction")

#Assign specific colors to each environmental variable
env_df_5$Color <- c("#E66100", "#9BA447", "#B01D8C") 

var_explained_5 <- c(0.40332*100, 0.35843*100)


# Create the biplot
biplot_5 <- ggplot() +
  # Plot species points
  geom_point(data = species_df_5, aes(x = CCA1, y = CCA2), color = "black", alpha=0.5, size = 4) +
  # Plot environmental arrows with custom colors
  geom_segment(data = env_df_5, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = Color),
               arrow = arrow(length = unit(1, "cm")), size = 3) +
  # Add arrow labels with custom colors - improved positioning and size
  geom_text(data = env_df_5, aes(x = CCA1 * 1.1, y = CCA2 * 1.2, label = Variable, color = Color),
            size = 8, fontface = "bold") +
  # Customize the plot
  theme_minimal() +
  coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +  # Set axis limits to -1.5 to 1.5
  # Add 0,0 cross line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey11") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey11") +
  # Add axis labels with variance explained
  labs(x = paste0("CCA1 (", round(var_explained_5[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_5[2], 1), "%)"),
       color = "Environmental Variables") +
  scale_color_identity() + 
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size = 20))

print(biplot_5)


### 2.6.5 Removing 10% ####
# Calculate presence percentage for each family
family_presence <- colSums(family_matrix > 0) / nrow(family_matrix) * 100

# Keep families present in more than 5% of samples
keep_families <- names(family_presence[family_presence > 10])
family_matrix_10 <- family_matrix[, keep_families]

# Run CCA
cca_result_10 <- vegan::cca(family_matrix_10 ~ predator * nutrients + Condition(Collection_month), 
                           data = env_data)

summary(cca_result_10)

# Run permutation test
set.seed(1)
anova_result_10 <- anova(cca_result_10, by = "terms", permutations = 999)

print(anova_result_10)

# Biplot 
# Extract scores
species_scores_10 <- scores(cca_result_10, display = "species")
env_scores_10 <- scores(cca_result_10, display = "bp")

# Convert scores to data frames for ggplot
species_df_10 <- as.data.frame(species_scores_10)
env_df_10 <- as.data.frame(env_scores_10)
env_df_10$Variable <- rownames(env_df_10)
env_df_10$Variable <- as.factor(env_df_10$Variable)
levels(env_df_10$Variable) <- c("Nutrient add.", "Predator rm", "Interaction")

#Assign specific colors to each environmental variable
env_df_10$Color <- c("#C14863", "#196971", "#BD974C")

var_explained_10 <- c(0.4248*100, 0.36095*100)

# Create the biplot
biplot_10 <- ggplot() +
  # Plot species points
  geom_point(data = species_df_10, aes(x = CCA1, y = CCA2), color = "black", alpha=0.5) +
  # Plot environmental arrows with custom colors
  geom_segment(data = env_df_10, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = Color),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.2) +
  # Add arrow labels with custom colors
  geom_text(data = env_df_10, aes(x = CCA1, y = CCA2, label = Variable, color = Color),
            vjust = -0.9, hjust = 0.5) +
  ggtitle("CCA removing family present in 10% of samples")+
  # Customize the plot
  theme_minimal() +
  coord_fixed() +
  # Add 0,0 cross line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey11") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey11") +
  # Add axis labels with variance explained
  labs(x = paste0("CCA1 (", round(var_explained_1[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_1[2], 1), "%)"),
       color = "Environmental Variables") +
  scale_color_identity() + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 11))


print(biplot_10)

## 2.7 Putting cca results together ####

### 2.7.1 results tables ####
# Function to extract relevant information from anova results
extract_anova_info <- function(anova_result) {
  # Convert anova result to data frame
  df <- as.data.frame(anova_result)
  
  # Add significance symbols
  df$Significance <- ""
  df$Significance[df$`Pr(>F)` <= 0.001] <- "***"
  df$Significance[df$`Pr(>F)` > 0.001 & df$`Pr(>F)` <= 0.01] <- "**"
  df$Significance[df$`Pr(>F)` > 0.01 & df$`Pr(>F)` <= 0.05] <- "*"
  df$Significance[df$`Pr(>F)` > 0.05 & df$`Pr(>F)` <= 0.1] <- "."
  
  # Round numeric columns
  df$F <- round(df$F, 3)
  df$`Pr(>F)` <- round(df$`Pr(>F)`, 3)
  
  return(df)
}

# Function to get model summary information
get_model_summary <- function(cca_result, family_matrix) {
  summary_res <- summary(cca_result)
  total_inertia <- sum(summary_res$tot.chi)
  constrained_pct <- summary_res$constr.chi / total_inertia * 100
  unconstrained_pct <- summary_res$unconst.chi / total_inertia * 100
  n_families <- ncol(family_matrix)
  
  return(c(
    total_inertia = round(total_inertia, 3),
    constrained_pct = round(constrained_pct, 3),
    unconstrained_pct = round(unconstrained_pct, 3),
    n_families = n_families
  ))
}

# Extract information from all your anova results
results_0 <- extract_anova_info(anova_result_0)
results_1 <- extract_anova_info(anova_result_1)
results_5 <- extract_anova_info(anova_result_5)
results_10 <- extract_anova_info(anova_result_10)

# Create a combined data frame for Table 1
table1_data <- data.frame(
  Threshold = c(rep("0%", nrow(results_0)), 
                rep("1%", nrow(results_1)),
                rep("5%", nrow(results_5)),
                rep("10%", nrow(results_10))),
  Term = c(rownames(results_0), rownames(results_1), 
           rownames(results_5), rownames(results_10)),
  Df = c(results_0$Df, results_1$Df, results_5$Df, results_10$Df),
  F = c(results_0$F, results_1$F, results_5$F, results_10$F),
  PrF = c(results_0$`Pr(>F)`, results_1$`Pr(>F)`, 
          results_5$`Pr(>F)`, results_10$`Pr(>F)`),
  Significance = c(results_0$Significance, results_1$Significance, 
                   results_5$Significance, results_10$Significance)
)

# Create Table 1 - ANOVA results
table1 <- kable(table1_data, 
                col.names = c("Threshold", "Term", "Df", "F", "Pr(>F)", "Significance"),
                caption = "ANOVA Results for CCA Models with Different Family Presence Thresholds",
                format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  group_rows("0% presence", 1, nrow(results_0)) %>%
  group_rows("1% presence", nrow(results_0) + 1, nrow(results_0) + nrow(results_1)) %>%
  group_rows("5% presence", nrow(results_0) + nrow(results_1) + 1, 
             nrow(results_0) + nrow(results_1) + nrow(results_5)) %>%
  group_rows("10% presence", nrow(results_0) + nrow(results_1) + nrow(results_5) + 1,
             nrow(results_0) + nrow(results_1) + nrow(results_5) + nrow(results_10))

# Print Table 1
print(table1)

# Get summary information for each model
summary_0 <- get_model_summary(cca_result_0, family_matrix_0)
summary_1 <- get_model_summary(cca_result_1, family_matrix_1)
summary_5 <- get_model_summary(cca_result_5, family_matrix_5)  # Assuming you have these variables
summary_10 <- get_model_summary(cca_result_10, family_matrix_10)  # Assuming you have these variables

# Create a data frame for Table 2
table2_data <- data.frame(
  Threshold = c("0%", "1%", "5%", "10%"),
  TotalInertia = c(summary_0["total_inertia"], summary_1["total_inertia"], 
                   summary_5["total_inertia"], summary_10["total_inertia"]),
  ConstrainedPct = c(summary_0["constrained_pct"], summary_1["constrained_pct"], 
                     summary_5["constrained_pct"], summary_10["constrained_pct"]),
  UnconstrainedPct = c(summary_0["unconstrained_pct"], summary_1["unconstrained_pct"], 
                       summary_5["unconstrained_pct"], summary_10["unconstrained_pct"]),
  NumFamilies = c(summary_0["n_families"], summary_1["n_families"], 
                  summary_5["n_families"], summary_10["n_families"])
)

# Create Table 2 - Model summary statistics
table2 <- kable(table2_data,
                col.names = c("Threshold", "Total Inertia", "Constrained (%)", 
                              "Unconstrained (%)", "No. of Families"),
                caption = "Model Summary Statistics",
                format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE)

# Print Table 2
print(table2)


# Save tables to files if needed
write.csv(table1_data, "table/cca_comparison_anova_results.csv", row.names = FALSE)
write.csv(table2_data, "table/cca_comparison_model_summary.csv", row.names = FALSE)

### 2.7.2 CCA plots ####

combined_plot <- grid.arrange(
  biplot_0, biplot_1, biplot_5, biplot_10,
  ncol = 2,
  top = textGrob("CCA Biplots with Different Family Presence Thresholds",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)


# 3. ranking affected families ####
family_scores <- scores(cca_result_5, display = "species") 
env_scores <- scores(cca_result_5, display = "bp") 

## 3.1 axis CCA1 ####
cca1<-sort(family_scores[,1])
barplot(cca1)
high_quantile <- quantile(cca1, 0.95)
low_quantile <- quantile(cca1, 0.05)
cca1<- cca1[cca1 >= high_quantile | cca1 <= low_quantile]
 
barplot(cca1, las=2)

data_cca1 <- data.frame(
  Family = names(cca1),  # Assuming the names of the family are row names
  Score = cca1
)

data_cca1$Family <- factor(data_cca1$Family, levels = data_cca1$Family[order(data_cca1$Score, decreasing = FALSE)])

data_cca1$Color <- ifelse(data_cca1$Score >= 0, "#9CAC99", "#C45F6B")  # Green for positive, red for negative

# ggplot2 bar plot with custom colors
cca1_plot <- ggplot(data_cca1, aes(x = Score, y = Family, fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  
  theme_minimal() +
  labs(title = "CCA1 axis",
       x = "Scores",
       y = "Families") + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

cca1_plot

## 3.2 axis CCA2 ####
cca2<-sort(family_scores[,2])
barplot(cca2)
high_quantile <- quantile(cca2, 0.95)
low_quantile <- quantile(cca2, 0.05)
cca2<- cca2[cca2 >= high_quantile | cca2 <= low_quantile]

barplot(cca2, las=2)

data_cca2 <- data.frame(
  Family = names(cca2),  # Assuming the names of the family are row names
  Score = cca2
)

data_cca2$Family <- factor(data_cca2$Family, levels = data_cca2$Family[order(data_cca2$Score, decreasing = FALSE)])

data_cca2$Color <- ifelse(data_cca2$Score >= 0, "#9CAC99", "#C45F6B")  # Green for positive, red for negative

# ggplot2 bar plot with custom colors
cca2_plot <- ggplot(data_cca2, aes(x = Score, y = Family, fill = Color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +  
  theme_minimal() +
  labs(title = "CCA2 axis",
       x = "Scores",
       y = "Families") + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12))

cca2_plot



################################################
# REDO Biplot and Family classification ####
#################################################
load("../Graph_settings.RData")

# 3. Ranking families ####
## 3.1 projection computation ####
# Calculate projections onto environmental vectors (using your better methodology)
family_scores <- scores(cca_result_5, display = "species") 
env_scores <- scores(cca_result_5, display = "bp") 

# Calculate projections for each treatment using correct row names
predator_projections <- family_scores %*% env_scores["predator0", ]
nutrients_projections <- family_scores %*% env_scores["nutrients1", ]
interaction_projections <- family_scores %*% env_scores["predator0:nutrients1", ]

# Convert to data frames and get top/bottom families
predator_effects <- data.frame(
  family = rownames(family_scores),
  projection = as.numeric(predator_projections)
) %>%
  arrange(desc(projection))

nutrients_effects <- data.frame(
  family = rownames(family_scores),
  projection = as.numeric(nutrients_projections)
) %>%
  arrange(desc(projection))

interaction_effects <- data.frame(
  family = rownames(family_scores),
  projection = as.numeric(interaction_projections)
) %>%
  arrange(desc(projection))

## 3.2 Top 3 families affected per treatment ####
# Top 3 positive and negative for each treatment
# Predator 
print("Predator - Positively associated (top 3):")
head(predator_effects, 3)
print("Predator - Negatively associated (bottom 3):")
tail(predator_effects, 3)

# Nutrients
print("Nutrients - Positively associated (top 3):")
head(nutrients_effects, 3)
print("Nutrients - Negatively associated (bottom 3):")
tail(nutrients_effects, 3)

# Interaction
print("Interaction - Positively associated (top 3):")
head(interaction_effects, 3)
print("Interaction - Negatively associated (bottom 3):")
tail(interaction_effects, 3)

# Get the families for plotting (3 top + 3 bottom for each treatment)
pred_families <- c(head(predator_effects$family, 3), tail(predator_effects$family, 3))
nutr_families <- c(head(nutrients_effects$family, 3), tail(nutrients_effects$family, 3))
interact_families <- c(head(interaction_effects$family, 3), tail(interaction_effects$family, 3))

# Fix: Get all unique affected families correctly
all_affected_families <- c(pred_families, nutr_families, interact_families) %>% unique()

# Check what we have
print(all_affected_families)

## 3.3 Dataset of affected families BIPLOT ####
# Create a comprehensive table showing which treatments affect each family
affected_families_table <- data.frame(
  Family = all_affected_families,
  Predator_effect = NA,
  Nutrients_effect = NA,
  Interaction_effect = NA,
  stringsAsFactors = FALSE
)

# Fill in the effects for each family
for(i in 1:nrow(affected_families_table)) {
  family_name <- affected_families_table$Family[i]
  
  # Check predator effects
  if(family_name %in% head(predator_effects$family, 3)) {
    affected_families_table$Predator_effect[i] <- "Positive"
  } else if(family_name %in% tail(predator_effects$family, 3)) {
    affected_families_table$Predator_effect[i] <- "Negative"
  } else {
    affected_families_table$Predator_effect[i] <- "None"
  }
  
  # Check nutrients effects
  if(family_name %in% head(nutrients_effects$family, 3)) {
    affected_families_table$Nutrients_effect[i] <- "Positive"
  } else if(family_name %in% tail(nutrients_effects$family, 3)) {
    affected_families_table$Nutrients_effect[i] <- "Negative"
  } else {
    affected_families_table$Nutrients_effect[i] <- "None"
  }
  
  # Check interaction effects
  if(family_name %in% head(interaction_effects$family, 3)) {
    affected_families_table$Interaction_effect[i] <- "Positive"
  } else if(family_name %in% tail(interaction_effects$family, 3)) {
    affected_families_table$Interaction_effect[i] <- "Negative"
  } else {
    affected_families_table$Interaction_effect[i] <- "None"
  }
}

# Print the table
print("Comprehensive table of affected families:")
print(affected_families_table)

# Save the table if needed
#write.csv(affected_families_table, "table/FINALaffected_families_summary.csv", row.names = FALSE)

# FIGURE BIPLOT FINAL ####
## 1. Prepare dataset ####
# Extract scores (using your cca_result_5)
species_scores_5 <- scores(cca_result_5, display = "species")
env_scores_5 <- scores(cca_result_5, display = "bp")

# Convert to data frames
species_df_5 <- as.data.frame(species_scores_5)
species_df_5$Family <- rownames(species_df_5)

# Add highlighting information
species_df_5$is_highlighted <- species_df_5$Family %in% all_affected_families

print("Number of highlighted families:")
print(sum(species_df_5$is_highlighted))

## 2. determine effect type on families ####
species_df_5$effect_type <- "Not highlighted"

for(i in 1:nrow(species_df_5)) {
  if(species_df_5$is_highlighted[i]) {
    family_name <- species_df_5$Family[i]
    effects <- c()
    
    # Check predator effects
    if(family_name %in% c(head(predator_effects$family, 3), tail(predator_effects$family, 3))) {
      effects <- c(effects, "predator")
    }
    # Check nutrients effects  
    if(family_name %in% c(head(nutrients_effects$family, 3), tail(nutrients_effects$family, 3))) {
      effects <- c(effects, "nutrients")
    }
    # Check interaction effects
    if(family_name %in% c(head(interaction_effects$family, 3), tail(interaction_effects$family, 3))) {
      effects <- c(effects, "interaction")
    }
    
    # Classify based on combination - UPDATED WITH & SYMBOLS
    if(length(effects) == 3) {
      species_df_5$effect_type[i] <- "Pred- & Nutr+ & Pred- × Nutr+"
    } else if(all(c("predator", "nutrients") %in% effects)) {
      species_df_5$effect_type[i] <- "Pred- & Nutr+"
    } else if(all(c("predator", "interaction") %in% effects)) {
      species_df_5$effect_type[i] <- "Pred- & Pred- × Nutr+"
    } else if(all(c("nutrients", "interaction") %in% effects)) {
      species_df_5$effect_type[i] <- "Nutr+ & Pred- × Nutr+"
    } else if("predator" %in% effects) {
      species_df_5$effect_type[i] <- "Pred- only"
    } else if("nutrients" %in% effects) {
      species_df_5$effect_type[i] <- "Nutr+ only"
    } else if("interaction" %in% effects) {
      species_df_5$effect_type[i] <- "Pred- × Nutr+ only"
    }
  }
}

print(table(species_df_5$effect_type))

## 3. prepare env data ####
env_df_5 <- as.data.frame(env_scores_5)
env_df_5$Variable <- rownames(env_df_5)
env_df_5$Variable <- case_when(
  env_df_5$Variable == "predator0" ~ "Pred-",
  env_df_5$Variable == "nutrients1" ~ "Nutr+", 
  env_df_5$Variable == "predator0:nutrients1" ~ "Pred- × Nutr+",
  TRUE ~ env_df_5$Variable
)

print("Environmental variables:")
print(env_df_5)

## 4. Calculate scaling and variance explained ####

# Scale arrows appropriately
arrow_scale <- 0.9 * max(abs(range(species_scores_5[,1:2]))) / max(abs(range(env_scores_5[,1:2])))

# Extract variance explained
# values from your CCA summary
var_explained_5 <- c(40.332, 35.843, 23.825)  # From "Proportion Explained" in constrained section

print("Variance explained:")
print(paste("CCA1:", round(var_explained_5[1], 1), "%"))
print(paste("CCA2:", round(var_explained_5[2], 1), "%"))


## 5. CCA BIPLOT Fig. 3 ####
# Update with your final palette
effect_colors <- c(
  "Pred- only" = "#ffc100",                    # Yellow
  "Nutr+ only" = "#86c75f",                    # Green
  "Pred- × Nutr+ only" = "#d35400",            # Orange
  "Pred- & Nutr+" = "#ff6b6b",                 # Coral/Pink
  "Pred- & Pred- × Nutr+" = "#c20089",         # Magenta
  "Nutr+ & Pred- × Nutr+" = "#9b59b6",         # Purple
  "Pred- & Nutr+ & Pred- × Nutr+" = "#0091a8", # Teal
  "Not highlighted" = "gray60"
)

# Make sure the factor levels are in the right order
species_df_5$effect_type <- factor(species_df_5$effect_type, 
                                   levels = c("Pred- only", "Nutr+ only", "Pred- × Nutr+ only",
                                              "Pred- & Nutr+", "Pred- & Pred- × Nutr+", 
                                              "Nutr+ & Pred- × Nutr+", "Pred- & Nutr+ & Pred- × Nutr+",
                                              "Not highlighted"))
print(table(species_df_5$effect_type))

# Use the correct variance values
var_explained_5 <- c(40.332, 35.843, 23.825)

# Scale arrows appropriately
arrow_scale <- 0.9 * max(abs(range(species_scores_5[,1:2]))) / max(abs(range(env_scores_5[,1:2])))

# Create the CCA biplot
CCA_biplot_final <- ggplot() +
  
  # Plot points of families unaffected 
  geom_point(data = species_df_5 %>% filter(!is_highlighted), 
             aes(x = CCA1, y = CCA2), 
             color = "gray50", alpha = 0.5, size = 2, shape = 16) +
  
  # Plot points of families affected 
  geom_point(data = species_df_5 %>% filter(is_highlighted), 
             aes(x = CCA1, y = CCA2, color = effect_type), 
             alpha = 0.8, size = 3.5, shape = 16, stroke = 0.5) +
  
  # Affected families text settings
  geom_text_repel(
    data = species_df_5 %>% filter(is_highlighted),
    aes(x = CCA1, y = CCA2, label = Family, color = effect_type),
    size = 3.5, max.overlaps = Inf, show.legend = FALSE, fontface = "bold",
    box.padding = 1.2, point.padding = 1, force = 4, min.segment.length = 0
  ) +
  
  # Environmental arrows
  geom_segment(data = env_df_5, 
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale * 0.8, yend = CCA2 * arrow_scale * 0.8),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 1.2, alpha = 0.9, color = "black") +
  
  # Arrow labels
  # In the arrow labels section, try:
  geom_text(data = env_df_5,
            aes(x = CCA1 * arrow_scale * 0.9, y = CCA2 * arrow_scale * 0.9, 
                label = Variable),
            size = 4, fontface = "bold", color = "black") +
  
  # Use your color palette
  scale_color_manual(values = effect_colors, name = "Treatment Response") +
  
  # Coordinate system and reference lines
  coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", alpha = 0.7, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", alpha = 0.7, linewidth = 0.5) +
  
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1.2, 0.8)) +
  
  # Apply your custom theme
  theme_jg_right +
  
  # Axis labels with correct variance
  labs(x = paste0("CCA1 (", round(var_explained_5[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_5[2], 1), "%)"),
       title = "") +
  
  # Fine-tune theme elements
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    
  )

# Print the plot
print(CCA_biplot_final)

ggsave("figures/M3-Fig.3_CCA_biplot_final.png", 
       plot = CCA_biplot_final,
       width = 9, height = 6, 
       dpi = 700, 
       units = "in")

ggsave("figures/M3-Fig.3_CCA_biplot_final.pdf", 
       plot = CCA_biplot_final,
       width = 9, height = 6, 
       dpi = 300, 
       units = "in")

#################
### CCA biplot for presentation ####
#################
# Assuming you've already run the ranking families code from your original script
# and have: cca_result_5, predator_effects, nutrients_effects, interaction_effects

# Extract scores
species_scores_5 <- scores(cca_result_5, display = "species")
env_scores_5 <- scores(cca_result_5, display = "bp")

# Convert to data frames
species_df_5 <- as.data.frame(species_scores_5)
species_df_5$Family <- rownames(species_df_5)

# Get affected families
pred_families <- c(head(predator_effects$family, 3), tail(predator_effects$family, 3))
nutr_families <- c(head(nutrients_effects$family, 3), tail(nutrients_effects$family, 3))
interact_families <- c(head(interaction_effects$family, 3), tail(interaction_effects$family, 3))

# Prepare environmental data
env_df_5 <- as.data.frame(env_scores_5)
env_df_5$Variable <- rownames(env_df_5)
env_df_5$Variable <- case_when(
  env_df_5$Variable == "predator0" ~ "Pred-",
  env_df_5$Variable == "nutrients1" ~ "Nutr+", 
  env_df_5$Variable == "predator0:nutrients1" ~ "Pred- × Nutr+",
  TRUE ~ env_df_5$Variable
)

# Calculate scaling
arrow_scale <- 0.9 * max(abs(range(species_scores_5[,1:2]))) / max(abs(range(env_scores_5[,1:2])))
var_explained_5 <- c(40.332, 35.843, 23.825)

# Base theme and settings
base_plot_settings <- list(
  coord_fixed(),
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", alpha = 0.7, linewidth = 0.5),
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", alpha = 0.7, linewidth = 0.5),
  scale_x_continuous(limits = c(-1, 1)),
  scale_y_continuous(limits = c(-1.2, 0.8)),
  theme_jg_right,
  labs(x = paste0("CCA1 (", round(var_explained_5[1], 1), "%)"),
       y = paste0("CCA2 (", round(var_explained_5[2], 1), "%)"),
       title = ""),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    legend.position = "right"
  )
)

################################
# PLOT 1: All families in grey
################################
species_df_plot1 <- species_df_5 %>%
  mutate(highlight = "Not highlighted")

CCA_biplot_1_all_grey <- ggplot() +
  # All families in grey
  geom_point(data = species_df_plot1, 
             aes(x = CCA1, y = CCA2), 
             color = "gray60", alpha = 0.6, size = 2.5, shape = 16) +
  
  # Environmental arrows
  geom_segment(data = env_df_5, 
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale * 0.8, yend = CCA2 * arrow_scale * 0.8),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 1.2, alpha = 0.9, color = "black") +
  
  # Arrow labels
  geom_text(data = env_df_5,
            aes(x = CCA1 * arrow_scale * 0.9, y = CCA2 * arrow_scale * 0.9, 
                label = Variable),
            size = 4, fontface = "bold", color = "black") +
  
  base_plot_settings

print(CCA_biplot_1_all_grey)

################################
# PLOT 2: Predator families in yellow
################################
species_df_plot2 <- species_df_5 %>%
  mutate(highlight = ifelse(Family %in% pred_families, "Predator", "Not highlighted"))

CCA_biplot_2_predator <- ggplot() +
  # Grey families
  geom_point(data = species_df_plot2 %>% filter(highlight == "Not highlighted"), 
             aes(x = CCA1, y = CCA2), 
             color = "gray60", alpha = 0.5, size = 2, shape = 16) +
  
  # Predator-affected families in yellow
  geom_point(data = species_df_plot2 %>% filter(highlight == "Predator"), 
             aes(x = CCA1, y = CCA2), 
             color = "#ffc100", alpha = 0.8, size = 3.5, shape = 16) +
  
  # Labels for predator families
  geom_text_repel(
    data = species_df_plot2 %>% filter(highlight == "Predator"),
    aes(x = CCA1, y = CCA2, label = Family),
    size = 3.5, max.overlaps = Inf, fontface = "bold",
    color = "#ffc100",
    box.padding = 1.2, point.padding = 1, force = 4, min.segment.length = 0
  ) +
  
  # Environmental arrows
  geom_segment(data = env_df_5, 
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale * 0.8, yend = CCA2 * arrow_scale * 0.8),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 1.2, alpha = 0.9, color = "black") +
  
  # Arrow labels
  geom_text(data = env_df_5,
            aes(x = CCA1 * arrow_scale * 0.9, y = CCA2 * arrow_scale * 0.9, 
                label = Variable),
            size = 4, fontface = "bold", color = "black") +
  
  base_plot_settings

print(CCA_biplot_2_predator)

################################
# PLOT 3: Nutrient families only (6 families in green)
# These are ONLY the 6 families most affected by nutrients
################################
species_df_plot3 <- species_df_5 %>%
  mutate(highlight = ifelse(Family %in% nutr_families, "Nutrients", "Not highlighted"))

CCA_biplot_3_nutrients <- ggplot() +
  # Grey families
  geom_point(data = species_df_plot3 %>% filter(highlight == "Not highlighted"), 
             aes(x = CCA1, y = CCA2), 
             color = "gray60", alpha = 0.5, size = 2, shape = 16) +
  
  # Nutrient-affected families in green
  geom_point(data = species_df_plot3 %>% filter(highlight == "Nutrients"), 
             aes(x = CCA1, y = CCA2), 
             color = "#86c75f", alpha = 0.8, size = 3.5, shape = 16) +
  
  # Labels for nutrient families
  geom_text_repel(
    data = species_df_plot3 %>% filter(highlight == "Nutrients"),
    aes(x = CCA1, y = CCA2, label = Family),
    size = 3.5, max.overlaps = Inf, fontface = "bold",
    color = "#86c75f",
    box.padding = 1.2, point.padding = 1, force = 4, min.segment.length = 0
  ) +
  
  # Environmental arrows
  geom_segment(data = env_df_5, 
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale * 0.8, yend = CCA2 * arrow_scale * 0.8),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 1.2, alpha = 0.9, color = "black") +
  
  # Arrow labels
  geom_text(data = env_df_5,
            aes(x = CCA1 * arrow_scale * 0.9, y = CCA2 * arrow_scale * 0.9, 
                label = Variable),
            size = 4, fontface = "bold", color = "black") +
  
  base_plot_settings +
  theme(legend.position = "none")

print(CCA_biplot_3_nutrients)

################################
# PLOT 4: Interaction families only (6 families in orange)
# These are ONLY the 6 families most affected by interaction
################################
species_df_plot4 <- species_df_5 %>%
  mutate(highlight = ifelse(Family %in% interact_families, "Interaction", "Not highlighted"))

CCA_biplot_4_interaction <- ggplot() +
  # Grey families
  geom_point(data = species_df_plot4 %>% filter(highlight == "Not highlighted"), 
             aes(x = CCA1, y = CCA2), 
             color = "gray60", alpha = 0.5, size = 2, shape = 16) +
  
  # Interaction-affected families in orange
  geom_point(data = species_df_plot4 %>% filter(highlight == "Interaction"), 
             aes(x = CCA1, y = CCA2), 
             color = "#d35400", alpha = 0.8, size = 3.5, shape = 16) +
  
  # Labels for interaction families
  geom_text_repel(
    data = species_df_plot4 %>% filter(highlight == "Interaction"),
    aes(x = CCA1, y = CCA2, label = Family),
    size = 3.5, max.overlaps = Inf, fontface = "bold",
    color = "#d35400",
    box.padding = 1.2, point.padding = 1, force = 4, min.segment.length = 0
  ) +
  
  # Environmental arrows
  geom_segment(data = env_df_5, 
               aes(x = 0, y = 0, xend = CCA1 * arrow_scale * 0.8, yend = CCA2 * arrow_scale * 0.8),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 1.2, alpha = 0.9, color = "black") +
  
  # Arrow labels
  geom_text(data = env_df_5,
            aes(x = CCA1 * arrow_scale * 0.9, y = CCA2 * arrow_scale * 0.9, 
                label = Variable),
            size = 4, fontface = "bold", color = "black") +
  
  base_plot_settings +
  theme(legend.position = "none")

print(CCA_biplot_4_interaction)

################################
# Save all plots
################################

ggsave("CCA_defense_1_all_grey.png", CCA_biplot_1_all_grey, 
       width = 9, height = 6, dpi = 300, bg = "white")
ggsave("CCA_defense_2_predator.png", CCA_biplot_2_predator, 
       width = 9, height = 6, dpi = 300, bg = "white")
ggsave("CCA_defense_3_nutrients.png", CCA_biplot_3_nutrients, 
       width = 9, height = 6, dpi = 300, bg = "white")
ggsave("CCA_defense_4_interaction.png", CCA_biplot_4_interaction, 
       width = 9, height = 6, dpi = 300, bg = "white")


################# 
# TABLES CCA ####
#################
# Function to extract CCA PERMANOVA results 
extract_cca_anova_results <- function(anova_result, model_name, phase = "CCA") {
  # Convert anova result to data frame
  anova_df <- as.data.frame(anova_result)
  anova_df$Effect <- rownames(anova_df)
  rownames(anova_df) <- NULL
  
  # Add model and phase info
  anova_df$Model <- model_name
  anova_df$Phase <- phase
  
  # Add significance symbols
  anova_df$Significance <- ""
  anova_df$Significance[anova_df$`Pr(>F)` <= 0.001] <- "***"
  anova_df$Significance[anova_df$`Pr(>F)` > 0.001 & anova_df$`Pr(>F)` <= 0.01] <- "**"
  anova_df$Significance[anova_df$`Pr(>F)` > 0.01 & anova_df$`Pr(>F)` <= 0.05] <- "*"
  anova_df$Significance[anova_df$`Pr(>F)` > 0.05 & anova_df$`Pr(>F)` <= 0.1] <- "."
  
  return(anova_df)
}

# Function to extract CCA model summary information
extract_cca_summary <- function(cca_result, family_matrix, model_name, phase = "CCA") {
  summary_res <- summary(cca_result)
  
  # Extract key statistics
  total_inertia <- sum(summary_res$tot.chi)
  constrained_inertia <- summary_res$constr.chi
  unconstrained_inertia <- summary_res$unconst.chi
  
  # Calculate percentages
  constrained_pct <- (constrained_inertia / total_inertia) * 100
  unconstrained_pct <- (unconstrained_inertia / total_inertia) * 100
  
  # Get number of families
  n_families <- ncol(family_matrix)
  
  # Get eigenvalues and their proportions
  eigenvals <- summary_res$cont$importance
  
  # Create summary data frame
  summary_df <- data.frame(
    Model = model_name,
    Phase = phase,
    Total_Inertia = round(total_inertia, 4),
    Constrained_Inertia = round(constrained_inertia, 4),
    Unconstrained_Inertia = round(unconstrained_inertia, 4),
    Constrained_Percent = round(constrained_pct, 2),
    Unconstrained_Percent = round(unconstrained_pct, 2),
    N_Families = n_families,
    CCA1_Eigenvalue = round(eigenvals[2, 1], 4),  # Proportion of Variance
    CCA2_Eigenvalue = round(eigenvals[2, 2], 4),  # Proportion of Variance
    CCA1_Cumulative = round(eigenvals[3, 1], 4),  # Cumulative Proportion
    CCA2_Cumulative = round(eigenvals[3, 2], 4),  # Cumulative Proportion
    stringsAsFactors = FALSE
  )
  
  return(summary_df)
}

# Function to extract axis loadings/scores information
extract_cca_loadings <- function(cca_result, model_name, phase = "CCA") {
  # Get biplot scores for environmental variables
  env_scores <- scores(cca_result, display = "bp")
  
  # Convert to data frame
  loadings_df <- as.data.frame(env_scores)
  loadings_df$Variable <- rownames(loadings_df)
  rownames(loadings_df) <- NULL
  
  # Add model and phase info
  loadings_df$Model <- model_name
  loadings_df$Phase <- phase
  
  # Clean up variable names based on your current model
  loadings_df$Variable <- case_when(
    loadings_df$Variable == "predator0" ~ "Predator",
    loadings_df$Variable == "nutrients1" ~ "Nutrients", 
    loadings_df$Variable == "predator0:nutrients1" ~ "Predator × Nutrients",
    TRUE ~ loadings_df$Variable
  )
  
  # Reorder columns
  loadings_df <- loadings_df[, c("Model", "Phase", "Variable", "CCA1", "CCA2")]
  
  return(loadings_df)
}

# Function to create family projections table
extract_family_projections <- function(cca_result, model_name, top_n = 3) {
  # Calculate projections (from your existing code)
  family_scores <- scores(cca_result, display = "species") 
  env_scores <- scores(cca_result, display = "bp") 
  
  # Calculate projections for each treatment
  predator_projections <- family_scores %*% env_scores["predator0", ]
  nutrients_projections <- family_scores %*% env_scores["nutrients1", ]
  interaction_projections <- family_scores %*% env_scores["predator0:nutrients1", ]
  
  # Create data frames
  predator_effects <- data.frame(
    family = rownames(family_scores),
    projection = as.numeric(predator_projections),
    treatment = "Predator",
    model = model_name
  ) %>%
    arrange(desc(projection)) %>%
    mutate(rank = row_number(),
           effect_direction = ifelse(rank <= top_n, "Top Positive", 
                                     ifelse(rank > (n() - top_n), "Top Negative", "Other")))
  
  nutrients_effects <- data.frame(
    family = rownames(family_scores),
    projection = as.numeric(nutrients_projections),
    treatment = "Nutrients",
    model = model_name
  ) %>%
    arrange(desc(projection)) %>%
    mutate(rank = row_number(),
           effect_direction = ifelse(rank <= top_n, "Top Positive", 
                                     ifelse(rank > (n() - top_n), "Top Negative", "Other")))
  
  interaction_effects <- data.frame(
    family = rownames(family_scores),
    projection = as.numeric(interaction_projections),
    treatment = "Predator × Nutrients",
    model = model_name
  ) %>%
    arrange(desc(projection)) %>%
    mutate(rank = row_number(),
           effect_direction = ifelse(rank <= top_n, "Top Positive", 
                                     ifelse(rank > (n() - top_n), "Top Negative", "Other")))
  
  # Combine all effects
  all_effects <- bind_rows(predator_effects, nutrients_effects, interaction_effects)
  
  # Filter to only top and bottom families
  top_bottom_effects <- all_effects %>%
    filter(effect_direction != "Other") %>%
    select(model, treatment, family, projection, rank, effect_direction) %>%
    arrange(treatment, effect_direction, rank)
  
  return(top_bottom_effects)
}

# TABLE S7: CCA Model Summary Statistics ####
cca_summary_5 <- extract_cca_summary(cca_result_5, family_matrix_5, "Model_5")

cca_summary_table <- cca_summary_5 %>%
  select(Model, Total_Inertia, Constrained_Percent, Unconstrained_Percent, 
         N_Families, CCA1_Eigenvalue, CCA2_Eigenvalue) %>%
  rename(
    "Model" = Model,
    "Total Inertia" = Total_Inertia,
    "Constrained (%)" = Constrained_Percent,
    "Unconstrained (%)" = Unconstrained_Percent,
    "No. Families" = N_Families,
    "CCA1 (Prop. Var.)" = CCA1_Eigenvalue,
    "CCA2 (Prop. Var.)" = CCA2_Eigenvalue
  )

# Create flextable
ft_cca_summary <- flextable(cca_summary_table) %>%
  colformat_double(digits = 3) %>%
  theme_vanilla() %>%
  add_header_lines("CCA Model Summary Statistics") %>%
  autofit()

print(ft_cca_summary)

# TABLE S8: CCA PERMANOVA Results ####
cca_anova_5 <- extract_cca_anova_results(anova_result_5, "Model_5")

# Add formatted p-values with significance
cca_anova_5$FormattedPvalue <- paste0(
  sprintf("%.3f", cca_anova_5$`Pr(>F)`), 
  " ", 
  cca_anova_5$Significance
)

cca_anova_table <- cca_anova_5 %>%
  select(Model, Phase, Effect, Df, F, FormattedPvalue) %>%
  rename(
    "Model" = Model,
    "Analysis" = Phase,
    "Effect" = Effect,
    "df" = Df,
    "F value" = F,
    "p-value" = FormattedPvalue
  ) %>%
  arrange(Effect)

# Create flextable
ft_cca_anova <- flextable(cca_anova_table) %>%
  colformat_double(digits = 3) %>%
  bold(i = grepl("\\*", cca_anova_table$`p-value`), j = "p-value") %>%
  theme_vanilla() %>%
  add_header_lines("CCA PERMANOVA Results") %>%
  autofit()

print(ft_cca_anova)

# TABLE : CCA Environmental Variable Loadings ####
cca_loadings_5 <- extract_cca_loadings(cca_result_5, "Model_5")

cca_loadings_table <- cca_loadings_5 %>%
  select(Model, Variable, CCA1, CCA2) %>%
  rename(
    "Model" = Model,
    "Environmental Variable" = Variable,
    "CCA1 Loading" = CCA1,
    "CCA2 Loading" = CCA2
  ) %>%
  arrange(`Environmental Variable`)

# Create flextable
ft_cca_loadings <- flextable(cca_loadings_table) %>%
  colformat_double(digits = 4) %>%
  theme_vanilla() %>%
  add_header_lines("CCA Environmental Variable Loadings") %>%
  autofit()

print(ft_cca_loadings)

# TABLE S8: Top Families Affected by Each Treatment
family_projections_5 <- extract_family_projections(cca_result_5, "Model_5", top_n = 3)

family_projections_table <- family_projections_5 %>%
  rename(
    "Model" = model,
    "Treatment" = treatment,
    "Family" = family,
    "Projection" = projection,
    "Rank" = rank,
    "Effect Direction" = effect_direction
  )

# Create flextable
ft_family_projections <- flextable(family_projections_table) %>%
  colformat_double(j = "Projection", digits = 4) %>%
  theme_vanilla() %>%
  add_header_lines("Top Families Affected by Each Treatment") %>%
  autofit()

print(ft_family_projections)

# SAVE TABLES (uncomment and modify paths as needed)
# Save tables
doc_cca_summary <- read_docx()
doc_cca_summary <- body_add_flextable(doc_cca_summary, ft_cca_summary)
print(doc_cca_summary, target = "table/CCA_summary.docx")

doc_cca_anova <- read_docx()
doc_cca_anova <- body_add_flextable(doc_cca_anova, ft_cca_anova)
print(doc_cca_anova, target = "table/CCA_permanova_results.docx")

doc_cca_loadings <- read_docx()
doc_cca_loadings <- body_add_flextable(doc_cca_loadings, ft_cca_loadings)
print(doc_cca_loadings, target = "tables/CCA_loadings.docx")

doc_family_projections <- read_docx()
doc_family_projections <- body_add_flextable(doc_family_projections, ft_family_projections)
print(doc_family_projections, target = "table/CCA_family_projections.docx")

# AFFECTED FAMILIES ####
## 3.4 Affected Families RANKING PLOT ####
# Create a comprehensive data frame with all families and their scores for each treatment
all_families_scores <- data.frame(
  Family = rownames(family_scores),
  Predator_score = as.numeric(predator_projections),
  Nutrients_score = as.numeric(nutrients_projections),
  Interaction_score = as.numeric(interaction_projections),
  stringsAsFactors = FALSE
)

# Filter to only include the affected families
affected_families_scores <- all_families_scores %>%
  filter(Family %in% all_affected_families)

# Create long format for plotting
affected_families_long <- affected_families_scores %>%
  pivot_longer(cols = c(Predator_score, Nutrients_score, Interaction_score),
               names_to = "Treatment",
               values_to = "Score") %>%
  mutate(
    Treatment = case_when(
      Treatment == "Predator_score" ~ "Pred-",
      Treatment == "Nutrients_score" ~ "Nutr+",
      Treatment == "Interaction_score" ~ "Pred- × Nutr+"
    ),
    Treatment = factor(Treatment, levels = c("Pred-", "Nutr+", "Pred- × Nutr+"))
  )

# Check the data
print(head(affected_families_long))

#3.5 Adding trophic information 
trophic_groups <- read.csv("FINALaffected_families_summary.csv", 
                           header=T, sep =";",
                           stringsAsFactors = FALSE)


fam_tg<-right_join(trophic_groups, affected_families_long)
str(fam_tg)

## 3.5 datapreparation for plotting ####
# Get top 3 and bottom 3 families for each treatment
top_families_per_treatment <- fam_tg %>%
  group_by(Treatment) %>%
  arrange(desc(Score)) %>%
  slice(c(1:3, (n()-2):n())) %>%  # Top 3 and bottom 3
  ungroup() %>%
  select(Family, Treatment) %>%
  mutate(is_top_affected = TRUE)

# Add the top_affected flag to the main dataset
fam_tg_plot <- fam_tg %>%
  left_join(top_families_per_treatment, by = c("Family", "Treatment")) %>%
  mutate(
    is_top_affected = ifelse(is.na(is_top_affected), FALSE, is_top_affected),
    alpha_level = ifelse(is_top_affected, 1.0, 0.3)  # Full opacity for top families, transparent for others
  )

fam_tg_plot$grp_trophic <- factor(fam_tg_plot$grp_trophic,
                                  levels = c("primary_producer", "mixotroph", "primary_consumer", 
                                             "secondary_consumer", "decomposer", "parasite" ), 
                                  labels =c("Primary Producer", "Mixotroph", "Primary Consumer", 
                                            "Secondary Consumer", "Decomposer", "Parasite"))


# Define colors for trophic groups
trophic_colors <- c(
  "Primary Consumer" = "#5f86b4",
  "Secondary Consumer" = "#f6b300",
  "Primary Producer" = "#4a7c59",
  "Mixotroph" = "#9cbe53",
  "Decomposer" = "#c0713e",
  "Parasite" = "#c40011"
)




## 3.6 Plot affected families scores ####
# Create the modified comparison plot
Fig.4 <- ggplot(fam_tg_plot, aes(x = Score, y = reorder(Family, Score), 
                                 fill = grp_trophic, alpha = I(alpha_level))) +
  geom_col(position = "identity", width = 0.7) +
  # Add zero reference line
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~Treatment, ncol = 3) +
  # Set consistent x-axis limits for all panels
  xlim(c(-0.9, 0.9)) +
  scale_fill_manual(values = trophic_colors, name = "") +
  theme_jg_right +
  labs(x = "Projection scores", y = "") +
  theme(
    strip.text = element_text(size = 16)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1)  # Show legend without transparency
    )
  )

Fig.4

ggsave("figures/M4-treatment_response_Familytrophic.png", 
       plot = Fig.4,
       width = 10, height = 6, 
       dpi = 700, 
       units = "in")

ggsave("figures/M4-treatment_response_Familytrophic.pdf", 
       plot = Fig.4,
       width = 16, height = 8, 
       dpi = 300, 
       units = "in")

ggsave("figures/M4-treatment_response_Familytrophic.tiff", 
       plot = Fig.4,
       width = 16, height = 8, 
       dpi = 300, 
       units = "in")

##########################
# 4. Sorensen index ####
###########################
## 4.1 Data preparation ####
family_treatments <- family_metadat_prednutr[, c(3, 7, 18, 27:303)]
family_treatments <- family_treatments %>%
  dplyr::select(!Percidae)
# Rename Treatments
family_treatments$Treatment <- factor(family_treatments$Treatment, levels = c("P", "none", "NP", "N"),
                                      labels = c("P(ctrl)", "PR", "N", "NPR"))

# Add this part - Filter out rare families (present in less than 5% of samples)
# Calculate presence frequency for each family
n_samples <- nrow(family_treatments)
presence_threshold <- 0.05 * n_samples  # 5% of samples

# Get columns with taxonomic data (excluding metadata)
taxa_columns <- family_treatments %>%
  select(-Collection_month, -Treatment, -Lake_ID) %>%
  colnames()

# Calculate presence frequency for each family
presence_counts <- sapply(taxa_columns, function(col) {
  sum(family_treatments[[col]] > 0)
})

# Identify families that meet the 5% threshold
families_to_keep <- names(presence_counts[presence_counts >= presence_threshold])

# Filter the dataset to keep only common families
community_data <- family_treatments %>%
  select(all_of(c(families_to_keep)))

metadata <- family_treatments %>%
  select(Collection_month, Treatment, Lake_ID)

# Optional: Print information about filtering
cat("Total families before filtering:", length(taxa_columns), "\n")
cat("Families retained after 5% threshold:", length(families_to_keep), "\n")
cat("Families removed:", length(taxa_columns) - length(families_to_keep), "\n")

## 4.2 Function to calculate within-treatment dissimilarity ####
calculate_one_combination <- function(comm_subset) {
  if(nrow(comm_subset) == 3) {
    dissim <- vegdist(comm_subset, method = "bray", binary = TRUE)
    return(as.vector(dissim))
  }
  return(NULL)
}

## 4.3 Calculate dissimilarities ####
results_list <- list()
counter <- 1

# Process each month-treatment combination
for(month in unique(metadata$Collection_month)) {
  for(treat in unique(metadata$Treatment)) {
    # Get indices for this combination
    indices <- metadata$Collection_month == month & 
      metadata$Treatment == treat
    
    # Get data for these samples
    current_comm <- community_data[indices, ]
    current_lakes <- metadata$Lake_ID[indices]
    
    # Only process if we have exactly 3 lakes
    if(nrow(current_comm) == 3) {
      # Calculate dissimilarity
      dissim <- calculate_one_combination(current_comm)
      
      if(!is.null(dissim)) {
        # Create lake pairs
        lake_pairs <- t(combn(current_lakes, 2))
        
        # Store results
        results_list[[counter]] <- data.frame(
          Collection_month = month,
          Treatment = treat,
          Lake_pair = paste(lake_pairs[,1], lake_pairs[,2], sep="_"),
          Dissimilarity = dissim
        )
        counter <- counter + 1
      }
    }
  }
}

# Combine results
dissimilarity_results <- do.call(rbind, results_list)

## 4.4 statistical analysis ####
# Ensure proper factor levels
dissimilarity_results$Treatment <- factor(dissimilarity_results$Treatment, 
                                          levels = c("P(ctrl)", "PR", "N", "NPR"))
dissimilarity_results$Collection_month <- factor(dissimilarity_results$Collection_month)

# Create binary columns for predator and nutrients
dissimilarity_results$predator <- ifelse(dissimilarity_results$Treatment %in% c("P(ctrl)", "N"), 1, 0)
dissimilarity_results$nutrients <- ifelse(dissimilarity_results$Treatment %in% c("N", "NPR"), 1, 0)

dissimilarity_results$predator <- factor(dissimilarity_results$predator, 
                                         levels = c(1, 0))

dissimilarity_results$nutrients <- factor(dissimilarity_results$nutrients, 
                                         levels = c(0, 1))

# Fit model
model <- glmmTMB(log(Dissimilarity) ~ predator * nutrients + (1|Collection_month),
                 data = dissimilarity_results)

simulateResiduals(model, plot=T)
summary(model)
Anova(model)

anova_results <- Anova(model)



## 4.5 Vizualization ####
## Estimated Marginal Means ####
# Create estimated marginal means
emm <- emmeans(model, ~ predator * nutrients)
emm_df <- as.data.frame(emm)

## Figure 5 - Dissimilarity Interaction Plot ####
# Convert labels for dissimilarity plot
predator_labels <- c("1" = "Pred+", "0" = "Pred-")
nutrient_labels <- c("0" = "Nutr-", "1" = "Nutr+")
nutrient_colors <- c("0" = "#3F72B3", "1" = "#9BA447")

emm_df$predator_label <- predator_labels[as.character(emm_df$predator)]
emm_df$predator_label <- factor(emm_df$predator_label, 
                                levels = c("Pred+", "Pred-"))

Fig.5 <- ggplot(emm_df, aes(x = predator_label, y = emmean, 
                                         color = nutrients, group = nutrients)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.15, linewidth = 1.2) +
  scale_color_manual(values = nutrient_colors, labels = nutrient_labels, name = "") +
  labs(x = "",
       y = "Log Estimated Sorensen Index") +
  theme_jg_right +
  theme(      
    legend.position = c(1.02, 1.02), 
    legend.justification = c("right", "top")
  )

Fig.5

# Export Figure 2 (aquatic only)
ggsave("figures/M5-Figure_5_heterogeneity_plot.png", 
       plot = Fig.5,
       width = 6, 
       height = 6, 
       units = "in",
       dpi = 700,
       bg = "white")

ggsave("figures/M5-Figure_5_heterogeneity_plot.pdf", 
       plot = Fig.5,
       width = 6, 
       height = 6, 
       units = "in",
       dpi = 300,
       bg = "white")

ggsave("figures/M5-Figure_5_heterogeneity_plot.tiff", 
       plot = Fig.5,
       width = 7, 
       height = 7, 
       units = "in",
       dpi = 300,
       bg = "white")

## TABLE S9: GLMM Model Results Summary - Dissimilarity Analysis ####

# Create tidy results for dissimilarity model
create_dissimilarity_glmm_table <- function() {
  
  # Single model for dissimilarity analysis
  model_list <- list(
    "Community Dissimilarity" = model
  )
  
  # Extract tidy results for the model
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
      "Analysis" = model,
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
    add_header_lines("GLMM Model Results: Community Dissimilarity Analysis") %>%
    width(j = "Analysis", width = 2) %>%
    width(j = "Term", width = 1.5) %>%
    autofit()
  
  # Save to Word
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft_models)
  print(doc, target = "table/TableS9_GLMM_dissimilarity_results.docx")
  
  # Also return the table for viewing
  return(ft_models)
}

## TABLE S10: ANOVA Results for Dissimilarity GLMM Model ####

create_dissimilarity_anova_table <- function() {
  
  # Single model for dissimilarity analysis
  model_list <- list(
    "Community Dissimilarity" = model
  )
  
  # Extract ANOVA results for the model
  anova_results <- map_dfr(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    anova_result <- car::Anova(model)
    anova_df <- as.data.frame(anova_result)
    anova_df$Effect <- rownames(anova_df)
    rownames(anova_df) <- NULL
    anova_df$Analysis <- model_name
    return(anova_df)
  })
  
  # Clean up column names
  names(anova_results)[names(anova_results) == "Chisq"] <- "Chi-Square"
  names(anova_results)[names(anova_results) == "Pr(>Chisq)"] <- "p-value"
  
  # Reorder columns
  anova_results <- anova_results[, c("Analysis", "Effect", "Chi-Square", "Df", "p-value")]
  
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
    add_header_lines("ANOVA Results: Community Dissimilarity Analysis") %>%
    width(j = "Analysis", width = 2) %>%
    width(j = "Effect", width = 1.5) %>%
    autofit()
  
  # Save to Word
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft_anova)
  print(doc, target = "table/TableS10_ANOVA_dissimilarity_results.docx")
  
  # Also return the table for viewing
  return(ft_anova)
}

## Execute the functions to create tables ####

# Create the tables
dissimilarity_glmm_table <- create_dissimilarity_glmm_table()
dissimilarity_anova_table <- create_dissimilarity_anova_table()

# Display the tables
dissimilarity_glmm_table
dissimilarity_anova_table
