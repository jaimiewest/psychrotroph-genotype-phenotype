# Sugar-Acid Preference based on KO genes
# per Gralka et al., 2023, sugar-acid preference: defined as [mean(k_sugars)-mean(k_acids)]/[mean(k_sugars)+mean(k_acids)]

# Towards bottom of this code file are correlation plots and boxplots for GC content, genome length

# Load  libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)
library(ggpubr)
library(car)       # for Levene’s test

setwd("path/to/folder")

############################################################################# 
########### Calculate observed SAP using Biolog substrate data ############
############################################################################# 

### Load in Biolog results
genIII = (read.csv(file ="Derived_data/GenIII_response_48_isolates.csv",
                   header=TRUE,check.names=FALSE))

# Pivot to long format
biolog_long = genIII %>%
  pivot_longer(
    cols = -Isolate,
    names_to = "substrate",
    values_to = "value")

### Load in plate substrate data, and merge in the S_or_A column (sugar/acid)
genIII_substrates = read.csv("GENIII_plate_substrates.csv",header=TRUE)

biolog_long_SA <- biolog_long %>%
  left_join(genIII_substrates %>% select(substrate, S_or_A), by = "substrate")

biolog_long_SA <- biolog_long_SA %>%
  mutate(
    S_or_A = case_when(
      S_or_A == "" ~ NA_character_,                   # Empty strings are NA
      S_or_A %in% c("organic acid", "amino acid") ~ "acid",  # Recode to simply acid
      TRUE ~ S_or_A))                                   # Keeps everything else as-is

# Remove rows with NA in S_or_A (these are the sensitivity assays on the plate), leaving you with 70 substrates
biolog_long_SA <- biolog_long_SA %>%  filter(!is.na(S_or_A))
unique(biolog_long_SA$substrate)


### Calculate empirical/observed SAP as log-ratio of sugar:acid use. 
# log-ratio is more statistically robust for comparing proportions
# Also, normalizing from -1 to 1 (as in Gralka paper) piles a lot of values up at 1, and we would 
# need to remove BA8 because it is an outlier that affects normalization.
Biolog_SAP <- biolog_long_SA %>%
  group_by(Isolate, S_or_A) %>%
  dplyr::summarize(S_or_A_Sum = sum(value), .groups = "drop") %>%
  pivot_wider(names_from = S_or_A, values_from = S_or_A_Sum) %>%
  mutate(
    # Normalize by total substrates per class
    S = sugar / 33,
    A = acid / 37,
    log_SAP = log2((S + 0.00001) / (A + 0.00001)),  # Add small pseudocount to avoid log(0)
  ) %>% filter(!Isolate %in% c("BA3", "BD2")) # remove outliers that were identified below (using Cook's distance)

############################################################################# 
########### Calculate predicted SAP based on KO genes in genomes ############
############################################################################# 

#############################################################################################
# ### First step - ALREADY DONE and DATABSE WAS UPDATED: 
# # We need to check for more KOs that perhaps are in our genomes but were not in Gralka's list.
# # First, Make a list of unique KOs with missing annotations
# missing_kos_df <- KO_merged %>%
#   filter(is.na(Metabolism) | is.na(Pathway) | is.na(KO_name)) %>%
#   distinct(KO)  # only need KO identifiers
# 
# # Pull annotations
# library(KEGGREST)
# 
# # Pull list of unique KOs with missing annotations
# missing_kos_df <- KO_merged %>%
#   filter(is.na(Metabolism) | is.na(Pathway) | is.na(KO_name)) %>%
#   distinct(KO)  # only need KO identifiers
# missing_kos <- unique(missing_kos_df$KO)
# 
# # Fetch KEGG entries (can be slow)
# ko_details <- lapply(missing_kos, function(ko) {
#   tryCatch(keggGet(ko)[[1]], error = function(e) return(NULL))})
# 
# # Extract KO names
# KO_name <- sapply(ko_details, function(entry) {
#   if (is.null(entry)) return(NA)
#   entry$NAME[1]})
# names(KO_name) <- missing_kos
# 
# # Extract Pathway(s)
# get_kegg_pathway <- function(entry) {
#   if (is.null(entry) || is.null(entry$PATHWAY)) return(NA)
#   paste(paste(names(entry$PATHWAY), entry$PATHWAY), collapse = "; ")}
# 
# Pathway <- sapply(ko_details, get_kegg_pathway)
# 
# # Extract BRITE category
# get_kegg_brite_category <- function(entry, ko_id) {
#   if (is.null(entry) || is.null(entry$BRITE)) {
#     return(NA)}
#   brite_lines <- entry$BRITE
#   brite_lines_clean <- gsub("^\\s+", "", brite_lines)
#   brite_levels <- nchar(gsub("(^\\s*).*", "\\1", brite_lines))
#   category_line <- brite_lines_clean[brite_levels == 2]
#   if (length(category_line) > 0) category_line[1] else NA}
# Category <- mapply(get_kegg_brite_category, ko_details, missing_kos, SIMPLIFY = TRUE)
# 
# # Combine everything into one data frame
# annotations <- data.frame(
#   KO = missing_kos,
#   KO_name = KO_name,
#   Pathway = Pathway,
#   Category = Category,
#   stringsAsFactors = FALSE)
# 
# #write.csv(annotations, "Derived_data/KOs_missing_from Gralka_database.csv", row.names = FALSE)
# # --> Manually review file for pathways that match pathways already in Gralka's list (could be automated!)
# 
# # Merge the 22 new KOs with Gralka's list
# KOs_sugar_vs_acid = (read.csv(file ="KOs_sugar_acid_Gralka2023.csv", header=TRUE))
# NEW_KOs_sugar_vs_acid = (read.csv(file ="KOs_missing_from Gralka_database.csv", header=TRUE))
# KOs_sugar_vs_acid2 = rbind(KOs_sugar_vs_acid, NEW_KOs_sugar_vs_acid)
# 
# #write.csv(KOs_sugar_vs_acid2, "KOs_sugar_acid_Gralka2023_plusDARPA.csv", row.names = FALSE)
#############################################################################################

# Load metadata
metadata <- read.csv("metadata_GenIII-6-6.csv",header=TRUE,row.names=1) 
# Merge with Biolog SAP df
Biolog_SAP2 <- Biolog_SAP %>%
  left_join(., metadata, by = "Isolate")
# Remove S and A columns - the ones that refer to Biolog substrate use
Biolog_SAP2 = Biolog_SAP2[,-c(2, 3, 4,5)]

# load KO gene data
KOs_in_isos = (read.csv(file ="Derived_data/Isolate_KEGGs_presence_absence.csv", header=TRUE))
#View(KOs_in_isos) # Isolate is first column, then >2900 columns, one for each KO

# load database of KO gene designations as sugar- or acid-preference
# This is Gralka's database plus some extra KOs added specific to the project - generated above
KOs_sugar_vs_acid = (read.csv(file ="KOs_sugar_acid_Gralka2023_plusDARPA.csv", header=TRUE))

# Pivot the KO dataframe to long format
df_long <- KOs_in_isos %>%
  pivot_longer(-Isolate, names_to = "KO", values_to = "Value") %>%
  filter(Value > 0)
# Merge with Gralka's (plus addiitonal) KO sugar/acid classification
KO_merged <- df_long %>%
  left_join(., KOs_sugar_vs_acid, by = "KO", relationship = "many-to-many")

# Remove rows with NA in S_or_A
KO_merged_clean <- KO_merged %>%
  filter(!is.na(S_or_A))
#View(KO_merged_clean) #5768 entries

# Now group and summarize proportion of S and A KO genes in each isolate
KO_SA_summarized = KO_merged_clean %>%
  group_by(Isolate, S_or_A) %>%
  dplyr::summarize(S_or_A_Sum = sum(Value)) %>%
  pivot_wider(names_from=S_or_A, values_from=S_or_A_Sum) %>%
  mutate(Total_KOs = acid + sugar) %>%
  mutate(
    S = sugar / Total_KOs,
    A = acid / Total_KOs)

# Merge with Biolog log_SAP data
Biolog_SAP3 <- Biolog_SAP2 %>%
  left_join(., KO_SA_summarized, by = "Isolate")

# Remove rows with NA in log_SAP
Biolog_SAP3 <- Biolog_SAP3 %>%
  filter(!is.na(log_SAP))

############################################################################# 
########### Calculate SAP a few ways to test different models ############
############################################################################# 

### 1. linear model
# Only use "S" (or "A") - they are collinear since S + A = 1, so one ends up NA in model if both are used.
cor(Biolog_SAP3$S, Biolog_SAP3$A)  # should be -1 since S + A = 1

#model_lm = lm(log_SAP ~ S + GC_contig1 + genome_length, data = Biolog_SAP3) 
model_lm = lm(log_SAP ~ S + GC_contig1, data = Biolog_SAP3) # genome_length excluded because it worsens the model
summary(model_lm)
AIC_lm = AIC(model_lm) # save AIC
Rsq_lm <- summary(model_lm)$r.squared # save R2


### 2. Polynomial model
#model_poly <- lm(log_SAP ~ S + I(S^2) + GC_contig1 + genome_length, data = Biolog_SAP3)
model_poly <- lm(log_SAP ~ S + I(S^2) + GC_contig1, data = Biolog_SAP3) # genome_length excluded because it worsens the model
summary(model_poly)
Rsq_poly <- summary(model_poly)$r.squared
AIC_poly <- AIC(model_poly)
#View(model_poly)

### 3. Nonlinear model (tanh, as used in Gralka et al paper)

#model_nls <- nls(log_SAP ~ tanh(s * S + b * GC_contig1 + c * genome_length),
model_nls <- nls(log_SAP ~ tanh(s * S + b * GC_contig1), # exclude genome_length because it worsens the model
                 data = Biolog_SAP3,
                 start = list(s = 0, b = 0))
summary(model_nls)
Biolog_SAP3$pred_nls <- predict(model_nls)
SStot <- sum((Biolog_SAP3$log_SAP - mean(Biolog_SAP3$log_SAP))^2)
SSres <- sum((Biolog_SAP3$log_SAP - Biolog_SAP3$pred_nls)^2)
Rsq_nls <- 1 - SSres / SStot
AIC_nls <- AIC(model_nls)


### 4. Compare models
results <- data.frame(
  Model = c("Linear", "Poly", "NLS (tanh)"),
  AIC   = c(AIC_lm, AIC_poly, AIC_nls),
  R2    = c(Rsq_lm, Rsq_poly, Rsq_nls))
print(results)


############################################################################# 
########### Test for outliers ############
############################################################################# 
# # Using the best model (Polynomial)
# library(car)
# 
# # Bonferroni-corrected test for outlying residuals
# # This tests whether the largest studentized residual is significantly different from the rest, with multiple testing correction.
# outlierTest(model_lm)
# outlierTest(model_poly)
# 
# # Visualize cook's distance to identify points that strongly influence the model fit
# # Index refers to row number of data.
# plot(cooks.distance(model_poly), type = "h", main = "Cook's Distance")
# abline(h = 4 / nrow(Biolog_SAP3), col = "red", lty = 2)  # rule of thumb threshold

# # Check against studentized residuals.
plot(rstudent(model_poly), main = "Studentized Residuals")
abline(h = c(-3, 3), col = "gray", lty = 2)

### Option to Return to top of script to exclude outliers, re-run analysis, if appropriate.

############################################################################# 
########### Test for inclusion of parameters ############
############################################################################# 
### Usign best model (Polynomial)
### For polynomial, compare the full model vs reduced (without GC content and genome length)

#model_poly_full2 <- lm(log_SAP ~ S + I(S^2) + GC_contig1 + genome_length, data = Biolog_SAP3)
#model_poly_justGenomeLength <- lm(log_SAP ~ S + I(S^2) + genome_length, data = Biolog_SAP3)
model_poly_full <- lm(log_SAP ~ S + I(S^2) + GC_contig1, data = Biolog_SAP3)
model_poly_reduced <- lm(log_SAP ~ S + I(S^2), data = Biolog_SAP3)

# Compare AIC and R²
AIC(model_poly_full)
summary(model_poly_full)$adj.r.squared

AIC(model_poly_reduced)
summary(model_poly_reduced)$adj.r.squared

# Formal model comparison (F-test)
# If F is insignificant, the models are similar.
anova(model_poly_full, model_poly_reduced)

# # Option to confirm with stepwise model selection
# step(model_poly_full, direction = "backward")

############################################################################# 
########### Plot predicted vs observed for best model ###########
############################################################################# 
# Use best model (polynomial with GC content, no genome length)
Biolog_SAP3$pred_poly <- predict(model_poly_full)

Biolog_SAP3$Phylum = ordered(Biolog_SAP3$Phylum, levels = c("Actinomycetota",
                                                            "Bacillota",
                                                            "Pseudomonadota",
                                                            "Bacteroidota"))

palette <- c(
  "Actinomycetota" = "#99cc33",
  "Bacillota" = "blue3",
  "Pseudomonadota" = "#cc3366",
  "Bacteroidota" = "gray60")

p = ggplot(Biolog_SAP3, aes(x = log_SAP, y = pred_poly, color = Phylum)) +
  geom_point() +
  stat_cor(aes(label =  after_stat(rr.label)), # Phylum-level R2
           label.y.npc = 0.935,
           size = 3.5
           ) + 
  stat_cor(data = Biolog_SAP3, aes(x = log_SAP, y = pred_poly, # OVERALL R² (black, single label)
                                   label =  after_stat(rr.label)),
                                   inherit.aes = FALSE,        # ignore color = Phylum
    color = "black",
    label.x.npc = "left",       # to tweak position 
    label.y.npc = 1
           ) +
  geom_smooth( ## OVERALL REGRESSION LINE (black)
    data = Biolog_SAP3, aes(x = log_SAP, y = pred_poly),
    inherit.aes = FALSE,
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  #geom_text(label = Biolog_SAP3$Isolate, size =3) +
  geom_abline( # 1:1 line
    slope = 1, intercept = 0, 
    linetype = "dashed", color = "lightgray", linewidth   = 0.3) +
  scale_color_manual(values = palette, name = "Phylum") +
  labs(
    #title = paste0("Polynomial Model Fit (R² = ", round(Rsq_poly, 2), ")"),
    x = "Observed SAP, GenIII substrate use (log-ratio)",
    y = "Predicted SAP (polynomial model)"
    ) +
  theme_bw() + theme(
    #panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill="white")) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))  # legend = dots only
p
#ggsave("Figures/SAP_predicted_vs_observed.tiff", width=5, height=5, units = "in", device='tiff', dpi=400)


############################################################################# 
########### Optional: identify high residuals ###########
############################################################################# 
# This ID's another 7 isolates; flag them in the figure, retain in analysis.

#library(ggrepel)
Biolog_SAP3 = Biolog_SAP3 %>%
  mutate(residuals = log_SAP - pred_poly)

p2 = ggplot(Biolog_SAP3, aes(x = pred_poly, y = residuals, color = Phylum)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = subset(Biolog_SAP3, abs(residuals) > 1.5 * sd(residuals)),
    aes(label = Isolate),
    size = 2,
    color = "black",
    min.segment.length = 0.1
  ) +
  scale_color_manual(values = palette, name = "Phylum") +
  # geom_text(aes(label = ifelse(abs(residuals) > 1.5 * sd(residuals), Isolate, "")),
  #           hjust = -0.25, vjust = 0.5, size = 2, color = "red") +
  labs(
    x = "Predicted SAP (polynomial model)",
    y = "Residuals (Observed - Predicted)"
    ) +
  theme_bw() + theme(#panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill="white"))
p2
#ggsave("Figures/SAP_residuals_observed-predicted.tiff", width=3.5, height=3, units = "in", device='tiff', dpi=400)


############################################################################# 
########### Plot GC content vs SAP ###########
############################################################################# 
### Plot GC content vs OBSERVED SAP
p = ggplot(Biolog_SAP3, aes(x = GC_contig1, y = log_SAP, color = Phylum)) +
  geom_point() +
  stat_cor(aes(label =  after_stat(rr.label)), # R2 per phylum
           label.x.npc = 0.6, label.y.npc = 0.935,
           size = 3.5) +
  stat_cor(data = Biolog_SAP3, aes(x = GC_contig1, y = log_SAP, # OVERALL R²
                                   label =  after_stat(rr.label)),
           inherit.aes = FALSE, # do this to ignore color = Phylum
           color = "black",
           label.x.npc = 0.6, label.y.npc = 1
  ) +
  geom_smooth( ## OVERALL REGRESSION LINE (black)
    data        = Biolog_SAP3,
    aes(x = GC_contig1, y = log_SAP),
    inherit.aes = FALSE,
    method      = "lm",
    se          = FALSE,
    color       = "black",
    linewidth   = 0.3
  ) +
  #geom_text(label = Biolog_SAP3$Isolate, size =3) +
  scale_color_manual(values = palette, name = "Phylum") +
  labs(
    x = "Genomic GC content",
    y = "Observed SAP, GenIII substrate use (log-ratio)"
    ) +
  theme_bw() + theme(#panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill="white"))
p
#ggsave("Figures/SAP_observed_vs_GC.tiff", width=3.35, height=5, units = "in", device='tiff', dpi=400)


### Plot GC content vs PREDICTED SAP
p = ggplot(Biolog_SAP3, aes(x = GC_contig1, y = pred_poly, color = Phylum)) +
  geom_point() +
  stat_cor(aes(label =  after_stat(rr.label)), # R2 per phylum
           label.y.npc = 0.885,
           size = 3.5) +
  stat_cor(data = Biolog_SAP3, aes(x = GC_contig1, y = pred_poly, # OVERALL R²
                                   label =  after_stat(rr.label)),
           inherit.aes = FALSE, # do this to ignore color = Phylum
           color = "black",
           label.y.npc = 0.95
  ) +

    geom_smooth( ## OVERALL REGRESSION LINE (black)
    data        = Biolog_SAP3,
    aes(x = GC_contig1, y = pred_poly),
    inherit.aes = FALSE,
    method      = "lm",
    se          = FALSE,
    color       = "black",
    linewidth   = 0.3
  ) +
  #geom_text(label = Biolog_SAP3$Isolate, size =3) +
  scale_color_manual(values = palette, name = "Phylum") +
  labs(
    x = "Genomic GC content",
    y = "Predicted SAP (polynomial model)"
  ) +
  theme_bw() + theme(#panel.grid.major = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill="white"))
p
#ggsave("Figures/SAP_predicted_vs_GC.tiff", width=3.35, height=5, units = "in", device='tiff', dpi=400)


############################################################################# 
########### Stats, boxplots describing SAP by phylum ###########
############################################################################# 

# Remove Bacteroidetes (n = 1), plus the NA rows
Biolog_SAP4 = subset(Biolog_SAP3, Phylum %in% c("Bacillota", "Pseudomonadota", "Actinomycetota"))

# Summarize SAP observed by phylum
Biolog_SAP4 %>%
  group_by(Phylum) %>%
  summarise(
    n = n(),
    mean_SAP_obs = mean(log_SAP, na.rm = TRUE),
    sd_SAP_obs = sd(log_SAP, na.rm = TRUE),
    median_SAP_obs = median(log_SAP, na.rm = TRUE),
    min_SAP_obs = min(log_SAP, na.rm = TRUE),
    max_SAP_obs = max(log_SAP, na.rm = TRUE))

# Phylum             n mean_SAP_obs sd_SAP_obs median_SAP_obs min_SAP_obs max_SAP_obs
# <ord>          <int>        <dbl>      <dbl>          <dbl>       <dbl>       <dbl>
#   1 Actinomycetota    16       -0.111      0.662         -0.194       -1.13       1.75 
# 2 Bacillota         16        0.816      1.37           0.618       -1.77       3.62 
# 3 Pseudomonadota    13       -0.847      0.817         -0.782       -2.16       0.165

# Summarize predicted SAP  by phylum
Biolog_SAP4 %>%
  group_by(Phylum) %>%
  summarise(
    n = n(),
    mean_SAP_pred = mean(pred_poly, na.rm = TRUE),
    sd_SAP_pred = sd(pred_poly, na.rm = TRUE),
    median_SAP_pred = median(pred_poly, na.rm = TRUE),
    min_SAP_pred = min(pred_poly, na.rm = TRUE),
    max_SAP_pred = max(pred_poly, na.rm = TRUE))

# Phylum             n mean_SAP_pred sd_SAP_pred median_SAP_pred min_SAP_pred max_SAP_pred
# <ord>          <int>         <dbl>       <dbl>           <dbl>        <dbl>        <dbl>
#   1 Actinomycetota    16        -0.294       0.647          -0.485       -1.05        1.41  
# 2 Bacillota         16         0.690       0.891           0.468       -0.271       3.15  
# 3 Pseudomonadota    13        -0.566       0.329          -0.687       -0.950       0.0763


######### Does SAP significantly differ by phylum?

kruskal.test(log_SAP ~ Phylum, data = Biolog_SAP4)

# Pairwise Wilcoxon with FDR correction
pairwise.wilcox.test(
  Biolog_SAP4$log_SAP,
  Biolog_SAP4$Phylum,
  p.adjust.method = "BH")

# Kruskal-Wallis rank sum test
# data:  log_SAP by Phylum
# Kruskal-Wallis chi-squared = 14.543, df = 2, p-value = 0.0006952

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  Biolog_SAP4$log_SAP and Biolog_SAP4$Phylum 
#                   Actinomycetota  Bacillota
# Bacillota         0.01555         -        
#   Pseudomonadota  0.06232         0.00093  
# P value adjustment method: BH 

kruskal.test(pred_poly ~ Phylum, data = Biolog_SAP4)

# Pairwise Wilcoxon with FDR correction
pairwise.wilcox.test(
  Biolog_SAP4$pred_poly,
  Biolog_SAP4$Phylum,
  p.adjust.method = "BH")

# Kruskal-Wallis rank sum test
# data:  pred_poly by Phylum
# Kruskal-Wallis chi-squared = 20.66, df = 2, p-value = 3.264e-05
# 
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  Biolog_SAP4$pred_poly and Biolog_SAP4$Phylum 

#               Actinomycetota    Bacillota
# Bacillota       0.00049            -        
# Pseudomonadota  0.32914            5.9e-06  
# P value adjustment method: BH 

###########################################################
# FIGURES: BOXplots for SAP
############################################################
# These are facet'ed so they look a little different than the GC and genome length boxplots..
str(Biolog_SAP4)

Biolog_SAP4_long <- Biolog_SAP4 %>%
  select(Isolate, Phylum, log_SAP, pred_poly) %>%
  pivot_longer(
    cols = c(log_SAP, pred_poly),
    names_to = "SAP_type",
    values_to = "SAP") %>%
  mutate(
    SAP_type = factor(
      SAP_type,
      levels = c("log_SAP", "pred_poly"),
      labels = c("Observed SAP", "Predicted SAP")))
Biolog_SAP4_long$Phylum = ordered(Biolog_SAP4_long$Phylum, levels = c("Bacillota", "Pseudomonadota", "Actinomycetota"))
  
p_facets <- ggplot(Biolog_SAP4_long,
                   aes(x = Phylum, y = SAP, fill = Phylum)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") + # To highlight zero
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  scale_fill_manual(values = palette, name = "Phylum") +
  scale_y_continuous(
    limits = c(-2.25, 4),
    breaks = -2:4,
    labels = -2:4
  ) +
  labs(x = NULL,y = "SAP") +
  facet_wrap(~ SAP_type) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1))
p_facets
#ggsave("Figures/SAP_boxplots_by_phylum.tiff", width=4, height=4, units = "in", device='tiff', dpi=400)



############################################################
######## TEST CORRELATIONS AMONG GC AND GENOME SIZE ########################
############################################################

cor_gc_genome <- cor.test(metadata$GC_contig1, metadata$genome_length, method = "spearman")
print(cor_gc_genome)

# Spearman's rank correlation rho
# data:  metadata$GC_contig1 and metadata$genome_length
# S = 22712, p-value = 0.03209
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.3131538 


############################################################
# FIGURES: BOXplots for GC, genome_length
############################################################
p1 = ggplot(metadata, aes(x = Phylum, y = GC_contig1, fill = Phylum)) +
    geom_boxplot(outlier.shape = NA, alpha=0.7) +
    geom_jitter(width = 0.15, size = 1, alpha=0.8) +
  scale_fill_manual(values = palette, name = "Phylum") +
    labs(x = NULL, y = "Genomic GC content") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle=25, hjust=1))
p1


p2 = ggplot(metadata, aes(x = Phylum, y = genome_length, fill = Phylum)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  geom_jitter(width = 0.15, size = 1, alpha=0.8) +
  scale_fill_manual(values = palette, name = "Phylum") +
  labs(x = NULL, y = "Genome length, Mbp") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=25, hjust=1))
p2

# Display together
p3 = ggarrange(p1, p2, ncol=2, nrow=1)
p3
#ggsave("Figures/GC_and_GenomeLength_by_phylum.tiff", width=4, height=4, units = "in", device='tiff', dpi=400)


############################################################
# FIGURE: Correlation plot for GC, genome_length, SAP
############################################################
#library(corrplot)

dat = read.csv(file = "Derived_data/DARPA_SAP_data_for_correlation.csv", header=TRUE)

# Remove the isolate name and number columns before running correlation on the df
dat = dat [,-c(1,2,3)]
dat$genome_length = as.numeric(dat$genome_length)

corr_mat = cor(dat, method = "spearman", 
               use = "pairwise.complete.obs") # This makes it work even if there are missing values in the input matrix
View(corr_mat)

# rename rows/cols of corr_mat
new_names <- c("Genome length", "GC content", "SAP, observed", "SAP, predicted")
rownames(corr_mat) <- new_names
colnames(corr_mat) <- new_names

tiff(
  filename   = "Figures/DARPA_SAP_correlation.tiff",
  width      = 5.5, # in inches
  height     = 5,
  units      = "in",
  res        = 600, # 600 dpi is journal quality
  #compression = "lzw"
)

heatmap <- corrplot(
  corr_mat,
  method = "color",
  type = "upper",
  tl.cex = 0.9,
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  diag=FALSE # drops the 1.00 diagonal
)
dev.off()

### End of script