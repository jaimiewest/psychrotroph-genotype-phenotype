# Load  libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(glmnet)
library(pROC)
library(caret) # For confusionmatrix
library(purrr)
library(forcats)
#install.packages("ggtext")
library(ggtext)

#BiocManager::install("KEGGREST")
library(KEGGREST)

setwd("path/to/dir")


# Load KEGG summary. This is actually copy numbers, not binary..
KEGGs = (read.csv(file ="Derived_data/KEGG_Ortholog_genes_48_isolates.csv",
                  header=TRUE))
genIII = test = (read.csv(file ="Derived_data/GenIII_response_48_isolates.csv",
                          header=TRUE,check.names=FALSE))

# Substrate will be the response variable; grab the substrates from the column names (and remove isolate column name)
substrate_columns = colnames(genIII)|> setdiff("Isolate")

# Merge datasets by isolate
data <- merge(genIII, KEGGs, by = "Isolate")

# Define the KO columns (predictors)
KO_columns <- grep("^K", names(data), value = TRUE)  # Select all columns that start with "KO"
length(KO_columns)

### Convert the KO columns binary, 0/1 presence/absence (not copy numbers)
# This can avoid overfitting or instability due to rare high-copy events, plus we have tons of KOs
# Could revisit count data to explore effect of copy number--perhaps for particularly interesting predictors.
data[KO_columns] <- (data[KO_columns] > 0) * 1

# Conservatively, remove KOs with very low frequency (only in 1 or 2 of the genomes
# and KOs that are low variance.
# This removes predictors that cannot help differentiate responses and reduces overfitting risk.

# How many isolates contained each KO gene?
KO_frequencies <- colSums(data[, KO_columns])  # assuming 0/1, negative/positive

# Remove KO's whose frequency is only 1,2 isolates
KO_columns <- KO_columns[KO_frequencies >= 3]
length(KO_columns)

# Check out variance for KOs that are left.
# Removing the "rare" KOs also removed the low variance KOs, so only a few more get filtered out in this step.
ko_variances <- apply(data[KO_columns], 2, var)
histogram(ko_variances)
KO_columns_filtered <- KO_columns[ko_variances > 0.05]  # Adjust threshold as needed
length(KO_columns_filtered)


# Set higher number of nested expressions to allow for the (up to) 2900 KO's 
# This step is perhaps unnecessary?
options(expressions = 5000)

# Originally, I checked for substrates that have < 8 presence or absence responses.
# However, with cross validation, we need to ensure that each fold has sufficient class balance.
# So, bumped it up to min per class of 18, which seems to be sufficient.
# Class imbalance in the response variable is problematic for glmnet. 
# In logistic regression, you need enough samples (eg isolates) in each class of the response variable.
# If, for a given substrate, only 3 isolates can grow on it, and 47 cannot, you have an imbalance. 
# The suggestion is to Skip substrates where either class (presence or absence) has fewer than 8 isolates
# Create function:
check_and_print_skipped_substrates <- function(data, substrate_columns, min_per_class = 18) {
  skipped <- c()
  for (substrate in substrate_columns) {
    vals <- data[[substrate]]
    if (!all(vals %in% c(0, 1))) next
    counts <- table(vals)
    if (any(counts < min_per_class)) {
      cat("Skipping:", substrate, "- class counts:", paste(counts, collapse = ", "), "\n")
      skipped <- c(skipped, substrate)
    }
  }
  return(skipped)
}

skipthese = check_and_print_skipped_substrates(data, substrate_columns)

# # Output 6/12/25, min = 18
# Skipping: Maltose - class counts: 15, 33 
# Skipping: D-Trehalose - class counts: 15, 33 
# Skipping: Sucrose - class counts: 15, 33 
# Skipping: D-Turanose - class counts: 17, 31 
# Skipping: pH 6.0 - class counts: 4, 44 
# Skipping: N-Acetyl-D-Mannosamine - class counts: 39, 9 
# Skipping: N-Acetyl-D-Galactosamine - class counts: 36, 12 
# Skipping: N-Acetyl-Neuraminic Acid - class counts: 36, 12 
# Skipping: 1% NaCl - class counts: 4, 44 
# Skipping: 4% NaCl - class counts: 17, 31 
# Skipping: a-D-Glucose - class counts: 7, 41 
# Skipping: D-Mannose - class counts: 13, 35 
# Skipping: D-Fructose - class counts: 7, 41 
# Skipping: 3-Methyl glucose - class counts: 45, 3 
# Skipping: D-Fucose - class counts: 41, 7 
# Skipping: L-Fucose - class counts: 37, 11 
# Skipping: L-Rhamnose - class counts: 32, 16 
# Skipping: 1% Sodium Lactate - class counts: 4, 44 
# Skipping: Fusidic Acid - class counts: 39, 9 
# Skipping: D-Mannitol - class counts: 12, 36 
# Skipping: Glycerol - class counts: 14, 34 
# Skipping: D-Glucose-6-Phosphate - class counts: 36, 12 
# Skipping: D-Fructose-6-Phosphate - class counts: 35, 13 
# Skipping: D-Aspartic Acid - class counts: 43, 5 
# Skipping: D-Serine.1 - class counts: 41, 7 
# Skipping: Troleandomycin - class counts: 33, 15 
# Skipping: Minocycline - class counts: 39, 9 
# Skipping: L-Aspartic Acid - class counts: 17, 31 
# Skipping: L-Glutamic Acid - class counts: 9, 39 
# Skipping: Niaproof 4 - class counts: 37, 11 
# Skipping: D-Galacturonic Acid - class counts: 32, 16 
# Skipping: D-Gluconic Acid - class counts: 10, 38 
# Skipping: Glucuronamide - class counts: 37, 11 
# Skipping: Mucic Acid - class counts: 37, 11 
# Skipping: D-Saccharic Acid - class counts: 32, 16 
# Skipping: Vancomycin - class counts: 31, 17 
# Skipping: Tetrazolium Violet - class counts: 33, 15 
# Skipping: 4-Hydroxyphenyl Acetic Acid - class counts: 34, 14 
# Skipping: Pyruvic Acid methyl ester - class counts: 15, 33 
# Skipping: D-Lactic Acid Methyl Ester - class counts: 34, 14 
# Skipping: L-Lactic Acid - class counts: 16, 32 
# Skipping: Citric Acid - class counts: 17, 31 
# Skipping: L-Malic Acid - class counts: 8, 40 
# Skipping: Nalidixic Acid - class counts: 17, 31 
# Skipping: Potassium Tellurite - class counts: 12, 36 
# Skipping: Tween 40 - class counts: 10, 38 
# Skipping: b-Hydroxybutyric Acid - class counts: 14, 34 
# Skipping: Acetoacetic Acid - class counts: 34, 14 
# Skipping: Acetic Acid - class counts: 17, 31 
# Skipping: Aztreonam - class counts: 4, 44 

# To set flag for close call substrates that may have funky folds:
# flagged_substrates = skipthese

# Remove these substrates from list of response variables, and "Isolate" column
substrate_columns = substrate_columns|> setdiff(skipthese)
substrate_columns #  44 substrates left at min = 18


?glm
?cv.glmnet
# https://glmnet.stanford.edu/articles/glmnet.html


### Function to run cross-validated, Elastic Net regularized logistic regression (or Lasso) for each substrate. 
# # Good for small sample, high-dimensional data (n ~ 50 isolates, >2000 predictors)
# # Uses cross-validation (3-fold) to choose best lambda (regularization strength), control overfitting

## Main function:
run_regularized_logistic <- function(substrate, KO_columns_filtered, data, alpha = 0.5) {
  # Define predictor matrix and response vector
  X <- as.matrix(data[KO_columns])
  y <- as.numeric(data[[substrate]])
  
  # Fit regularized logistic regression with 3-fold CV, using AUC
  cv_model <- cv.glmnet(
    X, y,
    family = "binomial",
    alpha = alpha,
    nfolds = 3,
    type.measure = "auc",
    keep = TRUE  # Need to "keep" to get fit.preval for fold-level AUC
  )
  
  # Best lambda from CV
  lambda.min <- cv_model$lambda.min
  
  # Fit final model on full data using best lambda
  final_model <- glmnet(X, y, family = "binomial", alpha = alpha, lambda = lambda.min)
  
  # Extract non-zero coefficients
  coef_df <- data.frame(
    KO = rownames(coef(final_model)),
    coefficient = as.vector(coef(final_model))
  )
  coef_df <- coef_df[coef_df$coefficient != 0, ]
  coef_df$Substrate <- substrate  # Add substrate
  
  # ---- AUC Calculations ----
  # Fold-level AUC (cross-validated predictions)
  lambda_idx <- which(cv_model$lambda == lambda.min)
  fold_preds <- cv_model$fit.preval[, lambda_idx]
  fold_auc <- as.numeric(pROC::auc(y, as.numeric(fold_preds)))
  
  # # Training AUC (model fit on full data); not very useful
  # train_preds <- predict(final_model, newx = X, type = "response")
  # train_auc <- as.numeric(pROC::auc(y, as.numeric(train_preds)))
  
  auc_df <- data.frame(
    Substrate = substrate,
    Fold_AUC = round(fold_auc, 3),
    Alpha = alpha,
    Lambda = lambda.min
  )
  
  # Return coefficients and AUC
  return(list(coef = coef_df, auc = auc_df))
}


## Create another function to run the above function, repeatedly (to identify the repeatable results):
run_repeated_elnet_regression <- function(substrate, KO_columns, data, alpha = 0.5, n_repeats = 500) {
  all_aucs <- numeric(n_repeats)
  all_lambdas <- numeric(n_repeats)
  coef_list <- list()
  
  for (i in 1:n_repeats) {
    result <- run_regularized_logistic(substrate, KO_columns, data, alpha = alpha)
    all_aucs[i] <- result$auc$Fold_AUC
    all_lambdas[i] <- result$auc$Lambda
    coef_list[[i]] <- result$coef
  }
  
  auc_stats <- data.frame(
    Substrate = substrate,
    Mean_Fold_AUC = round(mean(all_aucs, na.rm = TRUE), 3),
    SD_Fold_AUC = round(sd(all_aucs, na.rm = TRUE), 3),
    Mean_Lambda = mean(all_lambdas, na.rm = TRUE)
  )
  
  return(list(auc_stats = auc_stats, coefs = coef_list, auc_all = all_aucs))
}


#################################################################
### Run substrates through the repeated function, 
### which calls on the regularized logistic regression function
### using filtered KO list (dropped rare and low variance KOs)
#################################################################
regression_results <- lapply(substrate_columns, function(s) {
  run_repeated_elnet_regression(s, KO_columns_filtered, data)
})
head(data)

# Ideally you won't see any warnings triggered by class imbalance.

#Save all that data
#saveRDS(regression_results, file = "Derived_data/Full_regression_results_25June2025.rds")

# To load that back in:
#regression_results <- readRDS(file = "Derived_data/Full_regression_results_25June2025.rds")

names(regression_results) <- substrate_columns

#Save all AUC results
auc_all <- do.call(rbind, lapply(regression_results, `[[`, "auc_all"))
View(auc_all)
#write.csv(auc_all, "Derived_data/KO_REPEATED_500_all_AUC_25June2025.csv", row.names = TRUE)

# Summarize AUC results
auc_summary <- do.call(rbind, lapply(regression_results, `[[`, "auc_stats"))
View(auc_summary)
#write.csv(auc_summary, "Derived_data/KO_REPEATED_regression_AUC_25June2025.csv", row.names = FALSE)



# Make lists of substrates above certain mena fold AUC values
# use 0.65 to assess borderline substrates statistically, way below
# and establish the list ~0.7 for substrates that will be included in figures etc.
# (0.698, Li Chloride, was a significant result below, that's why the funny number..it rounds up to 0.70)
high_AUC_0.65_substrates <- auc_summary %>%
  filter(Mean_Fold_AUC >= 0.65) %>%
  pull(Substrate)
high_AUC_0.65_substrates

high_AUC_0.7_substrates <- auc_summary %>%
  filter(Mean_Fold_AUC >= 0.698) %>%
  pull(Substrate)
high_AUC_0.7_substrates

ko_summary <- regression_results[high_AUC_0.7_substrates] %>%
  map("coefs") %>%                  # Get the list of coefs for each repetition
  map_dfr(bind_rows) %>%            # Bind into one data frame
  filter(KO != "(Intercept)") %>%   # remove intercepts
  group_by(Substrate, KO) %>%
  summarize(
    selection_count = n(),          # count how many times that KO was in a model for that substrate
    mean_coefficient = mean(coefficient), # calculate the mean across all reps
    .groups = "drop"
  ) %>%
  arrange(desc(selection_count))

View(ko_summary)

# Save some data
#write.csv(ko_summary, "Derived_data/HighAucSubstrate_0.65_allKOs_meanCoefficients_25June2025.csv", row.names = FALSE)

# also saved alllll data, without selecting for the high AUCs:
#write.csv(ko_summary, "Derived_data/allSubstrates_allKOs_meanCoefficients_25June2025.csv", row.names = FALSE)


# Read in data, if not continuing from above

#ko_summary = (read.csv(file ="Derived_data/HighAucSubstrate_0.65_allKOs_meanCoefficients_25June2025.csv",
#                  header=TRUE))
#auc_summary = (read.csv(file ="Derived_data/KO_REPEATED_regression_AUC_25June2025.csv",
#                      header=TRUE))


# Interpreting AUC (Area Under the Receiving Operating Characteristic (ROC) Curve):
# The ROC curve plots the True positive rate (TPR) vs False PR (FPR); 
# ideally FPR is 0 and TPR is one, thus AUC = 1
# AUC Score	Interpretation:
# 1.0	Perfect prediction
# 0.9 – 1.0	Excellent
# 0.8 – 0.9	Good
# 0.7 – 0.8	Fair
# 0.6 – 0.7	Poor
# 0.5	No better than random guessing
# < 0.5	Worse than random

# # Option to check for balance of 0's and 1's for a substrate (eg Sucrose witch has a lot of negative coefficients)
# table(data$Sucrose)


# How many KOs are predictive of each substrate?
howmany = ko_summary %>%
  filter(selection_count >= 250) %>% # Must be present in at least half the runs
  group_by(Substrate) %>%
  dplyr::summarize(n_KOs = n_distinct(KO)) %>%
  arrange(desc(n_KOs))
howmany

howmany_pos_and_neg = ko_summary %>%
  filter(selection_count >= 250) %>% # Must be present in at least half the runs
  group_by(Substrate) %>%
  dplyr::summarize(
    n_pos_KOs = n_distinct(KO[mean_coefficient > 0]),  # Count distinct positive KOs
    n_neg_KOs = n_distinct(KO[mean_coefficient < 0])   # Count distinct negative KOs
  ) %>%
  arrange(desc(n_pos_KOs))
howmany_pos_and_neg

# Substrate              n_pos_KOs n_neg_KOs
# <chr>                      <int>     <int>
#   1 Sodium Butyrate               43        37
# 2 Sodium Bromate                32        63
# 3 Lincomycin                    31        46
# 4 g-Amino-N-Butyric Acid        31        20
# 5 a-Ketoglutaric Acid           30        51
# 6 Bromosuccinic Acid            22        42
# 7 D-Raffinose                   22         9
# 8 pH 5.0                        15        63
# 9 L-Histidine                   14        10
# 10 Propionic Acid                13        22
# 11 D-Malic Acid                  12        13
# 12 8% NaCl                       11        13
# 13 L-Arginine                     7         0
# 14 L-Serine                       3         1
# 15 Lithium Chloride               1         2


# Add AUC (which is already >0.7)
howmany_pos_and_neg = merge(howmany_pos_and_neg, auc_summary, by = "Substrate") 
howmany_pos_and_neg

### SAVE!
#write.csv(howmany_pos_and_neg, "Derived_data/KO_substrate_summary_table_25June2025.csv", row.names = FALSE)


#### Let's focus on the metabolism enabling KOs.
# Remove intercepts, and filter out negative and near-zero coefficients for stricter signal threshold. 
# Positive coefficient - KO is associated with growth on the substrate. Enabler
# Negative coefficient - KO is associated with lack of growth. Inhibitor of metabolism.
# Larger absolute coefficients = greater influence on the prediction.
#kos_to_focus_on$abs_coef <- abs(kos_to_focus_on$coefficient)

kos_to_focus_on <- ko_summary %>%
  subset(KO != "(Intercept)" & 
           mean_coefficient > 0.01 & # COMMENT OUT this line to include all coefficients
           selection_count >= 250) # aim for 50% of reps?
kos_to_focus_on <- kos_to_focus_on[order(-kos_to_focus_on$mean_coefficient), ]

# # Add AUC (which is already > 0.65 or 0.7)
kos_to_focus_on = merge(kos_to_focus_on, auc_summary, by = "Substrate") 

View(kos_to_focus_on)


# Check KO frequencies to be sure rare KOs are not causing the largest coefficients.
# Interpreting frequency:
# Frequent KO (at most, we have the same KO in 75% of isolates):
# Larger coefficients are broadly predictive of the phenotype.
# Smaller coefficients may be uninformative, non-specific.

# Rare KO (e.g., in <10% of isolates):
# May be a specialized advantage.
# A large coefficient could reflect a strong effect, but may also be unstable or overfitted.
# KOs found in few isolates can inflate coefficients due to small sample size; large coefficients for rare KOs may be statistical artifacts (unless biologically supported).
# Consider flagging KOs that are very rare before drawing strong conclusions.

# How many isolates contained each KO gene?
KO_frequencies <- colSums(data[, KO_columns_filtered])  # assuming 0/1 presence/absence
kos_to_focus_on$KO_freq <- KO_frequencies[kos_to_focus_on$KO]

# Let's add in substrate freqs too (how many isolates utilized that substrate?)
substrate_frequencies <- colSums(data[, substrate_columns])  # assuming 0/1 presence/absence
substrate_frequencies
kos_to_focus_on$substrate_freq <- substrate_frequencies[kos_to_focus_on$Substrate]
View(kos_to_focus_on)

# How many KOs are assigned to multiple substrates? 
# This is more interesting to consider the growth substrates separate from sensitivity, so first need to assign substrate types 
# (could bring this in from metadata, but there are few enough that this works fine)
growth_substrates <- c("Bromosuccinic Acid","Propionic Acid","a-Ketoglutaric Acid","D-Malic Acid","g-Amino-N-Butyric Acid","D-Raffinose","L-Arginine")
sensitivity_substrates <- c("Sodium Bromate","Lincomycin","Lithium Chloride","8% NaCl")

KO_multiplicity <- kos_to_focus_on %>%
  # add substrate_type first
  mutate(
    substrate_type = case_when(
      Substrate %in% growth_substrates ~ "Growth substrates",
      Substrate %in% sensitivity_substrates ~ "Sensitivity substrates",
      TRUE ~ NA_character_)) %>%
  filter(Substrate %in% high_AUC_0.7_substrates) %>%  # only AUC>0.7
  group_by(substrate_type, KO) %>%
  summarize(n_substrates_in_growth_or_sensitivity_grouping = n_distinct(Substrate),
    .groups = "drop") %>%
  arrange(substrate_type, desc(n_substrates_in_growth_or_sensitivity_grouping))
View(KO_multiplicity)

# # Do this if you prefer to lump togehter growth and sensitivity substrates
# KO_multiplicity2 = kos_to_focus_on %>% filter(Substrate %in% high_AUC_0.7_substrates) %>% 
#   group_by(KO) %>% dplyr::summarize(n_substrates = n_distinct(Substrate)) %>% arrange(desc(n_substrates))



#### PLOT of significant substrates and KO genes mean coefficient, 
# showing KOs predictive of multiple substrates:

# Merge back into main data
plot_data <- kos_to_focus_on %>%
  # add substrate_type first
  mutate(substrate_type = case_when(
      Substrate %in% growth_substrates ~ "Growth substrates",
      Substrate %in% sensitivity_substrates ~ "Sensitivity substrates", TRUE ~ NA_character_)) %>%
  filter(Substrate %in% high_AUC_0.7_substrates) %>% # To only report on AUC>0.7
  left_join(KO_multiplicity, by = c("substrate_type", "KO")) %>%
  
  # Insert line breaks into long substrate names before building Substrate_Label
  mutate(
    Substrate_clean = case_when(
      Substrate == "g-Amino-N-Butyric Acid" ~ "g-Amino-N-<br>Butyric Acid",
      Substrate == "Bromosuccinic Acid" ~ "Bromosuccinic<br>Acid",
      Substrate == "Propionic Acid" ~ "Propionic<br>Acid",
      Substrate == "a-Ketoglutaric Acid" ~ "a-Ketoglutaric<br>Acid",
      Substrate == "Sodium Bromate" ~ "Sodium<br>Bromate",
      Substrate == "Lithium Chloride" ~ "Lithium<br>Chloride",
      TRUE ~ Substrate
    )) %>%
  mutate(
    Substrate_Label = paste0("**", Substrate_clean, "**<br>",
                             #"n = ", substrate_freq, " isolates<br>",
                             "AUC=", sprintf("%.2f", Mean_Fold_AUC))
  ) %>%
  #arrange(desc(Mean_Fold_AUC)) %>%
  mutate(Substrate_Label = factor(Substrate_Label, levels = unique(Substrate_Label)))
View(plot_data)

# exact x-order you want...a messy way
x_levels <- c(
  "**Bromosuccinic<br>Acid**<br>AUC=0.73","**Propionic<br>Acid**<br>AUC=0.71",
  "**a-Ketoglutaric<br>Acid**<br>AUC=0.74","**D-Malic Acid**<br>AUC=0.75",
  "**g-Amino-N-<br>Butyric Acid**<br>AUC=0.73","**D-Raffinose**<br>AUC=0.75",
  "**L-Arginine**<br>AUC=0.77", "**Sodium<br>Bromate**<br>AUC=0.73",
  "**Lincomycin**<br>AUC=0.79","**Lithium<br>Chloride**<br>AUC=0.70","**8% NaCl**<br>AUC=0.72"
)
plot_data$Substrate_Label <- factor(plot_data$Substrate_Label, levels = x_levels)



# Create a fill column so specific KO genes are solid points and more general KO genes (that predict multiple substrates) are open
##########################
plot_data <- plot_data %>%
  mutate(n_substrates_cat = ifelse(n_substrates_in_growth_or_sensitivity_grouping == 1, "1", "2+"),
    fill_value = ifelse(n_substrates_cat == "1", "black", NA)  )

# Now plot
plot = ggplot(plot_data, aes(x = Substrate_Label, y = mean_coefficient)) +
  geom_jitter(
    width = 0.15,
    aes(size = as.factor(n_substrates_in_growth_or_sensitivity_grouping)),
    shape = 21, color = "black", stroke = 0.8,
    fill = ifelse(plot_data$n_substrates_in_growth_or_sensitivity_grouping == 1, "black", NA)) +
  scale_size_manual(
    name = "Open points indicate\nKO genes predictive \nof multiple substrates \n(within growth vs \nsensitivity grouping):",
    values = c("1" = 2, "2" = 2, "3" = 3, "4" = 4)) +
  theme_classic() +
  labs(x = "                            Growth substrates                                                                           Sensitivity substrates",
       y = "KO gene mean coefficient") +
  theme(axis.text.x = ggtext::element_markdown(size = 9),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10))
plot
#ggsave("Figures/KOgene_coefficients_per_substrate_23Nov2025.tiff", width=12.5, height=7, units = "in", device='tiff', dpi=400)
# Have to manually make the 1 dot in legend solid!

### Version of plot with all KOs:
# ...Have to go comment out the "mean_coefficient > 0.01 &" line above when generating kos_to_focus_on

plot = ggplot(plot_data, aes(x = Substrate_Label, y = mean_coefficient)) +
  geom_hline(yintercept = 0, color="gray")+
  geom_jitter(
    width = 0.25,
    aes(size = as.factor(n_substrates_in_growth_or_sensitivity_grouping)),
    shape = 21, color = "black", stroke = 0.8, alpha = 0.8,
    fill = ifelse(plot_data$n_substrates_in_growth_or_sensitivity_grouping == 1, "black", NA)) +
  scale_size_manual(
    name = "Open points indicate\nKO genes predictive \nof multiple substrates \n(within growth vs \nsensitivity grouping):",
    values = c("1" = 1, "2" = 1, "3" = 2, "4" = 2.5)) +
  theme_classic() +
  labs(x = "                            Growth substrates                                                                           Sensitivity substrates",
       y = "KO gene mean coefficient") +
  theme(axis.text.x = ggtext::element_markdown(size = 9),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))
plot

#ggsave("Figures/ALLL_KOgene_coefficients_per_substrate_23Nov2025.tiff", width=12, height=8, units = "in", device='tiff', dpi=400)
# Have to manually make the 1 dot in legend solid!





### Annotate KOs

# Step 1: Pull list of  KOs
ko_list <- unique(kos_to_focus_on$KO)
length(ko_list) # 209 KOs at 0.01... > 460 with all KOs

# Step 2: Fetch KEGG entries (can be slow)
ko_details <- lapply(ko_list, function(ko) {
  tryCatch(keggGet(ko)[[1]], error = function(e) return(NULL))
})

# Step 3: Extract KO names
KO_name <- sapply(ko_details, function(entry) {
  if (is.null(entry)) return(NA)
  entry$NAME[1]
})
names(KO_name) <- ko_list

# Step 4: Extract Pathway(s)
get_kegg_pathway <- function(entry) {
  if (is.null(entry) || is.null(entry$PATHWAY)) return(NA)
  paste(paste(names(entry$PATHWAY), entry$PATHWAY), collapse = "; ")
}

Pathway <- sapply(ko_details, get_kegg_pathway)

# Step 5: Extract BRITE category
get_kegg_brite_category <- function(entry, ko_id) {
  if (is.null(entry) || is.null(entry$BRITE)) {
    return(NA)
  }
  brite_lines <- entry$BRITE
  brite_lines_clean <- gsub("^\\s+", "", brite_lines)
  brite_levels <- nchar(gsub("(^\\s*).*", "\\1", brite_lines))
  category_line <- brite_lines_clean[brite_levels == 2]
  if (length(category_line) > 0) category_line[1] else NA
}

Category <- mapply(get_kegg_brite_category, ko_details, ko_list, SIMPLIFY = TRUE)

# Step 6: Combine everything into one data frame
annotations <- data.frame(
  KO = ko_list,
  KO_name = KO_name,
  Pathway = Pathway,
  Category = Category,
  stringsAsFactors = FALSE
)


kos_to_focus_on <- merge(kos_to_focus_on, annotations, by = "KO")
View(kos_to_focus_on)

#write.csv(kos_to_focus_on, "Derived_data/KO_annotations_coefficient_0.01_27June2025.csv", row.names = FALSE)

#kos_to_focus_on = (read.csv(file ="Derived_data/KO_annotations_coefficient_0.01_27June2025.csv",
#                      header=TRUE))







## Some basic code to help gut check these findings. Subset the full dataframe "data" to only show isolate of interest
# and the KOs that predict utilization

Lincomycin.KOs <- as.list(kos_to_focus_on$KO[kos_to_focus_on$Substrate == "Lincomycin"])
Lincomycin = data[,c("Isolate", "Lincomycin", paste(Lincomycin.KOs))]
View(Lincomycin)

##### Plot heatmap showing substrate use and KOs that predict it
library(pheatmap)

# Make sure there's only one row per Substrate–KO pair
kos_to_focus_on_clean <- kos_to_focus_on %>%
  group_by(Substrate, KO) %>%
  summarize(value = 1, .groups = "drop")  # or just: value = 1

# Pivot to Substrate x KO matrix
heatmap_df <- kos_to_focus_on_clean %>%
  pivot_wider(names_from = KO, values_from = value, values_fill = 0)

# Set rownames and convert to matrix
heatmap_matrix <- as.matrix(heatmap_df[,-1])
rownames(heatmap_matrix) <- heatmap_df$Substrate

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)


# Define a binary blue-white color scale
binary_colors <- c("0" = "white", "1" = "#08306b")  # dark blue

# Convert matrix values to character labels for custom legend (optional)
annotation_colors <- list(`KO gene` = c("Absent" = "white", "Present" = "#08306b"))


# ensure your matrix is binary (0s and 1s only)
heatmap_matrix_bin <- heatmap_matrix
heatmap_matrix_bin[heatmap_matrix_bin != 0] <- 1

# Custom colors
col_fun <- c("1" = "#08306b", "0" = "white")  # Present first

# Convert to matrix (if not already) and transpose
m <- t(as.matrix(heatmap_matrix_bin))

# Start the TIFF device
tiff("Figures/heatmap.tiff", width = 6, height = 12, units = "in", res = 400)

# Build heatmap
Heatmap(m,
        name = "KO gene",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_rot = -60,  # rotate x-axis labels
        border = TRUE,  # Axis lines around entire heatmap
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          at = c(1, 0),
          labels = c("Present", "Absent"),
          title = "KO gene",
          border = TRUE
        ))
dev.off() # Save the plot to tiff





###############################################
# Compare results to null model to assess how robust results are.
# Run your CV loop on shuffled substrate labels. AUC of null should be significantly lower than that of model.
###############################################

# Create function for null model
run_null_model_repeat <- function(substrate, KO_columns, data, alpha = 0.5, n_repeats = 500) {
  null_aucs <- numeric(n_repeats)
  
  for (i in 1:n_repeats) {
    y <- sample(data[[substrate]])  # permute which isolates are labeled as positive or negative for that substrate, while keeping the KO features unchanged. 
    #This destroys any real association between genotype and phenotype while preserving class balance
    
    X <- as.matrix(data[KO_columns])
    
    cv_model <- tryCatch({
      cv.glmnet(X, y, family = "binomial", alpha = alpha, nfolds = 3, type.measure = "auc", keep = TRUE)
    }, error = function(e) return(NULL))
    
    if (!is.null(cv_model)) {
      lambda_idx <- which(cv_model$lambda == cv_model$lambda.min)
      preds <- cv_model$fit.preval[, lambda_idx]
      auc_val <- tryCatch({
        as.numeric(pROC::auc(y, as.numeric(preds)))
      }, error = function(e) NA)
      null_aucs[i] <- auc_val
    } else {
      null_aucs[i] <- NA
    }
  }
  
  return(null_aucs)
}


# Loop through function, perhaps 500 times for comparability to original model function
n_repeats <- 500
pval_df <- list()
all_auc_df <- list()

for (substrate in high_AUC_0.65_substrates) {
  message("Processing ", substrate)
  
  # Use stored real AUCs from auc_all
  real_aucs <- as.numeric(auc_all[substrate, ])  # This gives a numeric vector of 500 AUCs
  
  # # Run tehse lines only if you don't have AUC stored from modeling:
  # real_result <- run_repeated_elnet_regression(substrate, KO_columns_filtered, data, alpha = 0.5, n_repeats = n_repeats)
  # real_aucs <- real_result$auc_all
  
  null_aucs <- run_null_model_repeat(substrate, KO_columns_filtered, data, 
                                     alpha = 0.5, n_repeats = n_repeats)
  
  # Store real and null AUCs in combined dataframe
  df <- data.frame(
    Substrate = substrate,
    AUC = c(real_aucs, null_aucs),
    Model = rep(c("Model", "Null"), each = n_repeats)
  )
  all_auc_df[[substrate]] <- df
   
  # Empirical p-value: mean-based
  mean_p_val <- mean(null_aucs >= mean(real_aucs), na.rm = TRUE)
  
  # Empirical p-value: median-based
  median_p_val <- mean(null_aucs >= median(real_aucs), na.rm = TRUE)
  
  # Store both
  pval_df[[substrate]] <- data.frame(
    Substrate = substrate,
    mean_p_value = mean_p_val,
    median_p_value = median_p_val)
}

# Combine list of p-value data frames into one
pvals_combined <- do.call(rbind, pval_df)
colnames(pvals_combined) <- c("Substrate", "empirical_p_mean", "empirical_p_median")
View(pvals_combined)

# Apply FDR correction to both sets of p-values
pvals_combined$adj_p_mean <- p.adjust(pvals_combined$empirical_p_mean, method = "BH")
pvals_combined$adj_p_median <- p.adjust(pvals_combined$empirical_p_median, method = "BH")



# Combine all AUC data into one data frame
all_auc_data <- do.call(rbind, all_auc_df)
View(all_auc_data)
all_auc_df


# Wilcoxon rank-sum test: to compare distribution of 500 real vs 500 null AUCs, per substrate
pvals_wilcox <- all_auc_data %>%
  group_by(Substrate) %>%
  summarise(wilcox_p = wilcox.test(AUC ~ Model)$p.value) %>%
  mutate(adj_p = p.adjust(wilcox_p, method = "BH"))  # Benjamini-Hochberg correction

pvals_wilcox


# Merge Wilcoxon p-values with empirical p-values
final_pval_df <- pvals_combined %>%
  left_join(pvals_wilcox, by = "Substrate")

final_pval_df <- final_pval_df %>%
  arrange(empirical_p_mean)

View(final_pval_df)


#write.csv(final_pval_df, "Derived_data/Pvalues_empirical_Wilcoxon_model_vs_null.csv", row.names = FALSE)

# Use empirical p-values for: Primary significance assessment 
# (especially if your real AUCs are summarized as mean or median).
# Easy explanation in a methods section (“permuted labels 500×, compared mean AUCs”).
# Empirical p-value tells you if your model performs better than chance (based on the mean or median).
# ...This is the proportion of null AUCs that are as extreme or more extreme than your observed mean AUC — i.e., how often you'd get an AUC as good as your real model just by chance.
# It’s a non-parametric, simulation-based estimate of significance.
# An empirical p-value is calculated by comparing the real test statistic (like mean AUC from your model) to the distribution of test statistics from the null model, 
# typically generated by label permutation.

# Use Wilcoxon test for:Distributional comparison across runs.
# Added support, e.g., to show separation isn't due to just one or two extreme model runs.
# Wilcoxon test tells you if the distribution of model AUCs is shifted relative to the null
# i.e., it uses the full distribution.

# with both: You have more interpretive flexibility
# e.g., one shows marginal significance, the other strong, then you report both.



# Add formatted labels
final_pval_df <- final_pval_df %>%
  mutate(label = ifelse(adj_p_median < 0.001, "***",
                        ifelse(adj_p_median < 0.01, "**",
                               ifelse(adj_p_median < 0.05, "*", "ns"))))

# Merge for plotting
all_auc_data <- left_join(all_auc_data, final_pval_df, by = "Substrate")
View(all_auc_data)



#######################
# Plot modeled AUC and null model AUC
#######################

# Order substrates by mean model AUC
substrate_order <- all_auc_data %>%
  filter(Model == "Model") %>%
  group_by(Substrate) %>%
  summarise(mean_AUC = mean(AUC, na.rm = TRUE)) %>%
  arrange(desc(mean_AUC)) %>%
  pull(Substrate)

all_auc_data <- all_auc_data %>%
  mutate(Substrate = factor(Substrate, levels = substrate_order))


my_labels_df <- all_auc_data %>%
  filter(Substrate %in% high_AUC_0.7_substrates) %>%
  filter(Model == "Model") %>%
  distinct(Model, Substrate, label.adj_p_median)
my_labels_df

# Make the boxplot
p <- ggplot(filter(all_auc_data, Substrate %in% high_AUC_0.7_substrates), aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(
    color = "black", width = 0.6, outlier.shape = NA
  ) +
  #geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c("Model" = "skyblue4", "Null" = "white")) +
  geom_text(
    data = my_labels_df,
    aes(x = 1, y = 0.94, label = label.adj_p_median),
    size = 6
  ) +
  facet_wrap(~ Substrate, scales = "fixed", nrow=2) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(x = "Model Type", y = "AUC") +
  ylim(0.33, 0.97)  # give space for stars
p


#ggsave("Figures/AUC_model_vs_null.tiff", width=10, height=4, units = "in", device='tiff', dpi=400)
