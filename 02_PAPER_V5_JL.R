# Load Required Libraries ----

# Data Manipulation and Input/Output
library(tidyverse)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(magrittr)
library(reshape2)

# Statistical Modeling and Analysis
library(lme4)
library(lmerTest)
library(survival)
library(survminer)
library(coxme)
library(cluster)
library(rstatix)
library(pROC)
library(moments)
library(caret)
library(glmnet)
library(car)
library(performance)
library(broom.mixed)
library(randomForest)

# Data Visualization
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(viridis)
library(scales)
library(factoextra)
library(dendextend)
library(gridExtra)
library(patchwork)
library(corrplot)

# Reporting and Tables
library(gt)
library(knitr)
library(kableExtra)
library(officer)

# Imputation
library(mice)
library(naniar)

# Define Base Path ----

base_path <- file.path(".")

# Load and Combine Data ----

# Load datasets
samples <- read_excel(file.path(base_path, "samples.xlsx"))
patients <- read_excel(file.path(base_path, "patients.xlsx"))

# Remove duplicate columns before merging (except 'nhc')
common_cols <- intersect(names(samples), names(patients))
samples_reduced <- samples %>%
  select(-one_of(setdiff(common_cols, "nhc")))

# Merge datasets on 'nhc'
merged <- left_join(samples_reduced, patients, by = "nhc")

# Write merged data to an Excel file
write_xlsx(merged, file.path(base_path, "merged.xlsx"))

# Data Adjustments ----

# Convert specific 'x' values to NA
merged <- merged %>%
  mutate(across(c(v1_tn, v1_n, v1_m), ~ na_if(., "x")))

# Convert date columns
date_vars <- c("lb_date", "birthdate", "followup_end_date", "qx_r0_prev1_date", "relapse_date", "distantmetastasis_date")
merged <- merged %>%
  mutate(across(all_of(date_vars), as.Date))

# Clean and convert numeric variables
numeric_vars <- c("ldh", "cfdna_concentration", "s100")
merged <- merged %>%
  mutate(across(all_of(numeric_vars), ~ parse_number(as.character(.))))

# Convert percentage and size variables to numeric
pct_avg_sz_vars <- names(merged) %>%
  str_subset("^(pct_|avg_).*|.*_sz$")
merged <- merged %>%
  mutate(across(all_of(pct_avg_sz_vars), as.numeric))

# Convert other specific variables to numeric
merged <- merged %>%
  mutate(
    lb_order_nhc = as.numeric(lb_order_nhc),
    lb_nhc_count = as.numeric(lb_nhc_count),
    age = as.numeric(age),
    followup_time_months = ifelse(followup_time_months < 0, 0, followup_time_months)
  )

# Update 'lb_type' and 'clinical_group' for controls
control_subsets <- c("displastic_nevi", "congenital_nevus", "melanocytoma")
merged <- merged %>%
  mutate(
    lb_type = ifelse(subset %in% control_subsets, "control", lb_type),
    clinical_group = ifelse(subset %in% control_subsets, "control", clinical_group)
  )

# Convert variables to ordered factors
ordered_factors <- list(
  v1_tn = c("0", "1", "2", "3", "4"),
  v1_tl = c("a", "b"),
  v1_n = c("0", "1a", "1b", "1c", "2a", "2b", "2c", "3a", "3b", "3c"),
  v1_m = c("0", "1a", "1b", "1c"),
  disease_activity = c("AAC", "NED", "PR", "SD", "PD"),
  v1_staging_n = c("0", "I", "II", "III", "IV"),
  lb_type = c("control", "preop", "postop", "followup", "relapse", "unresect"),
  clinical_group = c("NA", "control", "A", "B", "C", "D", "E"),
  disease_cat = c("none", "tumor", "regional_micro", "regional_macro", "metastasis"),
  disease_max = c("none", "T1", "T2", "T3", "T4", "N1a", "N1b", "N1c", "N2a", "N2b", "N2c",
                  "N3a", "N3b", "N3c", "M1a", "M1b", "M1c", "M1d")
)

for (var in names(ordered_factors)) {
  merged[[var]] <- factor(merged[[var]], levels = ordered_factors[[var]], ordered = TRUE)
}

# Convert remaining variables to unordered factors
factor_vars <- c(
  "disease", "status", "clinical_status", "treatment", "subset", "primary_type",
  "primary_location", "prior_relapse", "adjuvancy_postv1", "status_braf",
  "status_nras", "status_tert", "preoperatory_lb", "relapse_status", "relapse_disease"
)
merged <- merged %>%
  mutate(across(all_of(factor_vars), as.factor))


'''
# Handle Missing Data for S100 and LDH ----

# Select variables for imputation
impute_vars <- c("s100", "ldh", "age", "sex", "clinical_status", "v1_staging_n")
impute_data <- merged %>%
  select(all_of(impute_vars)) %>%
  mutate(
    sex = as.factor(sex),
    clinical_status = as.factor(clinical_status),
    v1_staging_n = as.factor(v1_staging_n)
  )

# Initialize MICE
init <- mice(impute_data, maxit = 0, print = FALSE)
meth <- init$method
predM <- init$predictorMatrix

# Customize imputation methods
meth["s100"] <- "pmm"
meth["ldh"]  <- "pmm"

# Customize predictor matrix
predM["s100", "ldh"] <- 1
predM["ldh", "s100"] <- 1

# Prevent variables from predicting themselves
predM["s100", "s100"] <- 0
predM["ldh", "ldh"] <- 0

# Perform multiple imputation
set.seed(123)
imputed <- mice(impute_data, method = meth, predictorMatrix = predM, m = 1, seed = 123, print = FALSE)

# Extract the imputed dataset
completed <- complete(imputed, 1)

# Integrate imputed values back into the original dataset
merged <- merged %>%
  mutate(
    s100 = if_else(is.na(s100), completed$s100, s100),
    ldh  = if_else(is.na(ldh),  completed$ldh,  ldh)
  )

# Verify that there are no missing values in 's100' and 'ldh'
cat("Missing s100 after imputation:", sum(is.na(merged$s100)), "\n")
cat("Missing ldh after imputation:", sum(is.na(merged$ldh)), "\n")
'''
# Fragmentomics Filter and Imputation ----

# Filter the dataset for high-risk fragment samples
highrisk_fragment_samples <- merged %>%
  filter(
    !is.na(ngu_r20_150),                 # Exclude rows with NA in 'ngu_r20_150'
    is.na(comorbidities),                # Include only rows where 'comorbidities' is NA
    !subset %in% c("thin_mm", "merkel")  # Exclude specific subsets
  )

# Define fragmentomics variables to impute with half the minimum value
fragmentomics_to_impute_a <- c(
  "ngu_r250_320", "ngu_r320_700", "ngu_r700_1300",
  "pct_tot_r250_320", "pct_tot_r320_700", "pct_tot_r700_1300",
  "nml_r250_320", "nml_r320_700", "nml_r700_1300"
)

# Impute ngu_*, pct_tot_*, and nml_* variables with half the minimum non-zero value
for (var in fragmentomics_to_impute_a) {
  min_val <- highrisk_fragment_samples %>%
    filter(.data[[var]] > 0) %>%
    summarise(min_val = min(.data[[var]], na.rm = TRUE)) %>%
    pull(min_val)
  
  half_min <- min_val / 2
  
  highrisk_fragment_samples <- highrisk_fragment_samples %>%
    mutate(
      !!var := if_else(is.na(.data[[var]]), half_min, .data[[var]])
    )
  
  cat(paste0("Imputed ", var, " with value: ", half_min, "\n"))
}

# Define variables for imputation based on low ngu_* values
pct_cv_vars <- c("pct_cv_r250_320", "pct_cv_r320_700", "pct_cv_r700_1300")
avg_sz_vars <- c("avg_sz_r250_320", "avg_sz_r320_700", "avg_sz_r700_1300")
ngu_vars <- c("ngu_r250_320", "ngu_r320_700", "ngu_r700_1300")
ranges <- c("r250_320", "r320_700", "r700_1300")

# Calculate median pct_cv and avg_sz within low ngu_* subsets for each range
imputation_medians <- list()

for (i in seq_along(ranges)) {
  range <- ranges[i]
  ngu_var <- ngu_vars[i]
  pct_cv_var <- pct_cv_vars[i]
  avg_sz_var <- avg_sz_vars[i]
  
  # Calculate the minimum non-zero ngu value
  min_ngu <- highrisk_fragment_samples %>%
    filter(.data[[ngu_var]] > 0) %>%
    summarise(min_val = min(.data[[ngu_var]], na.rm = TRUE)) %>%
    pull(min_val)
  
  # Define threshold as 10% above the minimum ngu value
  threshold <- min_ngu * 1.1
  
  # Filter observations with ngu near the lower limit
  low_ngu_data <- highrisk_fragment_samples %>%
    filter(.data[[ngu_var]] <= threshold)
  
  # Calculate median avg_sz and pct_cv within low ngu data
  median_avg_sz <- median(low_ngu_data[[avg_sz_var]], na.rm = TRUE)
  median_pct_cv <- median(low_ngu_data[[pct_cv_var]], na.rm = TRUE)
  
  # Store the medians
  imputation_medians[[range]] <- list(
    avg_sz = median_avg_sz,
    pct_cv = median_pct_cv
  )
}

# Impute missing pct_cv and avg_sz with the calculated medians
for (i in seq_along(ranges)) {
  range <- ranges[i]
  pct_cv_var <- pct_cv_vars[i]
  avg_sz_var <- avg_sz_vars[i]
  
  median_avg_sz <- imputation_medians[[range]]$avg_sz
  median_pct_cv <- imputation_medians[[range]]$pct_cv
  
  # Impute avg_sz
  highrisk_fragment_samples <- highrisk_fragment_samples %>%
    mutate(
      !!avg_sz_var := if_else(is.na(.data[[avg_sz_var]]), median_avg_sz, .data[[avg_sz_var]])
    )
  
  # Impute pct_cv
  highrisk_fragment_samples <- highrisk_fragment_samples %>%
    mutate(
      !!pct_cv_var := if_else(is.na(.data[[pct_cv_var]]), median_pct_cv, .data[[pct_cv_var]])
    )
  
  # Print the imputed values
  cat(paste0("Imputed ", avg_sz_var, " with value: ", median_avg_sz, "\n"))
  cat(paste0("Imputed ", pct_cv_var, " with value: ", median_pct_cv, "\n"))
}

# Verify that there are no remaining missing values in the specified fragmentomics variables
missing_after_imputation <- highrisk_fragment_samples %>%
  select(all_of(c(fragmentomics_to_impute_a, pct_cv_vars, avg_sz_vars))) %>%
  summarise(across(everything(), ~ sum(is.na(.x))))

print(missing_after_imputation)

# Compute Proportional Fragment Variables ----

# Calculate 'prop_fragment' and its logarithmic transformation
highrisk_fragment_samples <- highrisk_fragment_samples %>%
  mutate(
    prop_fragment = pct_tot_r20_150 / pct_tot_r160_180,
    prop_fragment2 = pct_tot_r20_150 / pct_tot_r180_220,
    prop_fragment3 = pct_tot_r100_150 / pct_tot_r163_169,
    prop_fragment_log = log(prop_fragment + 0.01)
  )

# Define Fragmentomics Variables ----

fragmentomics_variables <- c(
  # pct_tot and props
  "pct_tot_r20_150", "pct_tot_r160_180", "pct_tot_r180_220",
  "pct_tot_r250_320", "pct_tot_r320_700", "pct_tot_r700_1300",
  "prop_fragment", "prop_fragment2", "prop_fragment3",
  
  # avg_size (global and each range)
  "avg_size", "avg_sz_r20_150", "avg_sz_r160_180",
  "avg_sz_r180_220", "avg_sz_r250_320", "avg_sz_r320_700",
  "avg_sz_r700_1300",
  
  # pct_cv
  "pct_cv_r20_150", "pct_cv_r160_180", "pct_cv_r180_220",
  "pct_cv_r250_320", "pct_cv_r320_700", "pct_cv_r700_1300"
)

# Create readable labels for variables
variable_labels <- c(
  "pct_tot_r20_150" = "Total % (20-150 bp)",
  "avg_sz_r20_150" = "Average Size (20-150 bp)",
  "pct_cv_r20_150" = "CV % (20-150 bp)",
  
  "pct_tot_r160_180" = "Total % (160-180 bp)",
  "avg_sz_r160_180" = "Average Size (160-180 bp)",
  "pct_cv_r160_180" = "CV % (160-180 bp)",
  
  "pct_tot_r180_220" = "Total % (180-220 bp)",
  "avg_sz_r180_220" = "Average Size (180-220 bp)",
  "pct_cv_r180_220" = "CV % (180-220 bp)",
  
  "pct_tot_r250_320" = "Total % (250-320 bp)",
  "avg_sz_r250_320" = "Average Size (250-320 bp)",
  "pct_cv_r250_320" = "CV % (250-320 bp)",
  
  "pct_tot_r320_700" = "Total % (320-700 bp)",
  "avg_sz_r320_700" = "Average Size (320-700 bp)",
  "pct_cv_r320_700" = "CV % (320-700 bp)",
  
  "pct_tot_r700_1300" = "Total % (700-1300 bp)",
  "avg_sz_r700_1300" = "Average Size (700-1300 bp)",
  "pct_cv_r700_1300" = "CV % (700-1300 bp)",
  
  "avg_size" = "Global Average Size (bp)",
  "prop_fragment" = "Ratio (20-150/160-180)",
  "prop_fragment2" = "Ratio (20-150/180-220)",
  "prop_fragment3" = "Ratio (100-150/163-169)"
)

# Save the Results ----

write_xlsx(highrisk_fragment_samples, file.path(base_path, "fragment.xlsx"))

# Table 1: Patient Demographics ----

# Summarize patient data
patient_summary <- highrisk_fragment_samples %>%
  arrange(nhc, lb_date) %>%
  group_by(nhc) %>%
  summarise(
    Age_at_Diagnosis = first(age),
    Sex = first(sex),
    Clinical_Group = first(clinical_group),
    Staging_N = first(v1_staging_n),
    Status_BRAF = first(status_braf),
    Primary_Location = first(primary_location),
    Primary_Type = first(primary_type),
    Adjuvancy_Post_V1 = first(adjuvancy_postv1),
    Treatment = first(treatment),
    Prior_Relapse = first(prior_relapse),
    Relapse_Status = first(relapse_status),
    .groups = 'drop'
  )

# Filter out patients with missing Clinical_Group
patient_summary_filtered <- patient_summary %>%
  filter(!is.na(Clinical_Group))

# Define Subgroups
control_group <- patient_summary_filtered %>%
  filter(Clinical_Group == "control")

resectable_group <- patient_summary_filtered %>%
  filter(Clinical_Group %in% c("A", "B", "C", "D"))

unresectable_group <- patient_summary_filtered %>%
  filter(Clinical_Group == "E")

# Helper Functions for Summary Statistics

# Function to calculate median and standard deviation of age
age_median_sd <- function(data) {
  median_age <- median(data$Age_at_Diagnosis, na.rm = TRUE)
  sd_age <- sd(data$Age_at_Diagnosis, na.rm = TRUE)
  paste0(median_age, " (±", round(sd_age, 1), ")")
}

# Function to calculate percentage of a specific sex
sex_percentage <- function(data, sex_value = "male") {
  percentage <- mean(data$Sex == sex_value, na.rm = TRUE) * 100
  paste0(round(percentage, 1), "%")
}



# Function to summarize categorical variables
categorical_summary <- function(data, variable) {
  if (nrow(data) == 0 || all(is.na(data[[variable]]))) return("N/A")
  prop_table <- prop.table(table(data[[variable]])) * 100
  prop_strings <- paste0(names(prop_table), ": ", round(prop_table, 1), "%")
  paste(prop_strings, collapse = "; ")
}

# Create the Descriptive Table

clinicopathological_table <- tibble(
  Characteristic = c(
    "Subgroup (n)", "Age at Diagnosis", "Sex (Male, %)", "Staging (N, %)", "BRAF Status (%)",
    "Primary Location (%)", "Primary Type (%)", "Adjuvancy Post V1 (%)", "Treatment (%)",
    "Prior Relapse (%)", "Relapse Status (%)"
  ),
  
  All_Patients = c(
    nrow(patient_summary_filtered),
    age_median_sd(patient_summary_filtered),
    sex_percentage(patient_summary_filtered, "male"),
    categorical_summary(patient_summary_filtered, "Staging_N"),
    categorical_summary(patient_summary_filtered, "Status_BRAF"),
    categorical_summary(patient_summary_filtered, "Primary_Location"),
    categorical_summary(patient_summary_filtered, "Primary_Type"),
    categorical_summary(patient_summary_filtered, "Adjuvancy_Post_V1"),
    categorical_summary(patient_summary_filtered, "Treatment"),
    categorical_summary(patient_summary_filtered, "Prior_Relapse"),
    categorical_summary(patient_summary_filtered, "Relapse_Status")
  ),
  
  Control = c(
    nrow(control_group),
    if(nrow(control_group) > 0) age_median_sd(control_group) else "N/A",
    if(nrow(control_group) > 0) sex_percentage(control_group, "male") else "N/A",
    "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"
  ),
  
  Resectable = c(
    nrow(resectable_group),
    age_median_sd(resectable_group),
    sex_percentage(resectable_group, "male"),
    categorical_summary(resectable_group, "Staging_N"),
    categorical_summary(resectable_group, "Status_BRAF"),
    categorical_summary(resectable_group, "Primary_Location"),
    categorical_summary(resectable_group, "Primary_Type"),
    categorical_summary(resectable_group, "Adjuvancy_Post_V1"),
    categorical_summary(resectable_group, "Treatment"),
    categorical_summary(resectable_group, "Prior_Relapse"),
    categorical_summary(resectable_group, "Relapse_Status")
  ),
  
  Unresectable = c(
    nrow(unresectable_group),
    age_median_sd(unresectable_group),
    sex_percentage(unresectable_group, "male"),
    categorical_summary(unresectable_group, "Staging_N"),
    categorical_summary(unresectable_group, "Status_BRAF"),
    categorical_summary(unresectable_group, "Primary_Location"),
    categorical_summary(unresectable_group, "Primary_Type"),
    categorical_summary(unresectable_group, "Adjuvancy_Post_V1"),
    categorical_summary(unresectable_group, "Treatment"),
    categorical_summary(unresectable_group, "Prior_Relapse"),
    "N/A"  # Relapse Status for Unresectable is N/A
  )
)

# Display the Table

print(kable(clinicopathological_table, caption = "Clinicopathological Characteristics of Patients"))

# Export Table 1 to Word
doc_table1 <- read_docx() %>%
  body_add_par("Table 1: Clinicopathological Characteristics of Patients", style = "heading 1") %>%
  body_add_table(clinicopathological_table, style = "table_template")

# Save the document
print(doc_table1, target = "Table1.docx")

# Statistical Analysis and Tables ----

# Define the function to calculate mean and standard deviation and format them
calculate_mean_sd <- function(data, variable) {
  mean_value <- mean(data[[variable]], na.rm = TRUE)
  sd_value <- sd(data[[variable]], na.rm = TRUE)
  
  # Custom rounding based on variable type
  if (grepl("avg_sz|avg_size", variable)) {
    formatted_value <- paste0(round(mean_value, 0), " ± ", round(sd_value, 0))
  } else if (variable %in% c("cfdna_concentration", "prop_fragment")) {
    formatted_value <- paste0(round(mean_value, 2), " ± ", round(sd_value, 2))
  } else {
    formatted_value <- paste0(round(mean_value, 1), " ± ", round(sd_value, 1))
  }
  
  return(formatted_value)
}

# TABLE S1a & S1b ----

# Analyze fragmentomics variables by clinical status

# Create a grouping variable 'group' for clinical status
data_clinical_status <- highrisk_fragment_samples %>%
  mutate(group = case_when(
    lb_type == "control" ~ "Control",
    clinical_status == "AWoD" ~ "AWoD",
    clinical_status == "AWD" ~ "AWD",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Calculate sample sizes for each group
clinical_status_groups <- unique(data_clinical_status$group)
sample_sizes <- setNames(sapply(clinical_status_groups, function(cs) {
  sum(data_clinical_status$group == cs)
}), clinical_status_groups)
total_sample_size <- nrow(data_clinical_status)

# Initialize result matrix
result_matrix <- matrix(ncol = length(clinical_status_groups) + 1, nrow = length(fragmentomics_variables))
colnames(result_matrix) <- c(paste0("All samples (n=", total_sample_size, ")"),
                             paste0(clinical_status_groups, " (n=", sample_sizes[clinical_status_groups], ")"))
rownames(result_matrix) <- fragmentomics_variables

# Perform statistical analysis
kruskal_p_values <- numeric(length(fragmentomics_variables))
names(kruskal_p_values) <- fragmentomics_variables
pairwise_p_values_list <- list()

for (variable in fragmentomics_variables) {
  # Mean ± SD for all samples
  result_matrix[variable, paste0("All samples (n=", total_sample_size, ")")] <- calculate_mean_sd(data_clinical_status, variable)
  
  # Mean ± SD for each group
  for (cs in clinical_status_groups) {
    data_subset <- data_clinical_status[data_clinical_status$group == cs, ]
    result_matrix[variable, paste0(cs, " (n=", sample_sizes[cs], ")")] <- calculate_mean_sd(data_subset, variable)
  }
  
  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(as.formula(paste(variable, "~ group")), data = data_clinical_status)
  kruskal_p_values[variable] <- kruskal_test$p.value
}

# Adjust p-values and format
adjusted_p_values <- p.adjust(kruskal_p_values, method = "BH")
formatted_p_values <- sapply(adjusted_p_values, function(p) {
  if (p < 0.01) "<0.01" else sprintf("%.2f", p)
})
result_matrix <- cbind(result_matrix, "Adjusted p-value" = formatted_p_values)

# Convert to data frame and add variable labels
result_df <- as.data.frame(result_matrix)
if (exists("variable_labels")) {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable") %>%
    mutate(Variable = variable_labels[Variable])
} else {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable")
}

# Print and save the main summary table
print(kable(result_df, caption = "Table S1a: Mean ± SD for All Samples and Clinical Status Groups with Adjusted p-values", align = 'c'))
doc_s1a <- read_docx() %>%
  body_add_par("Table S1a: Mean ± SD for All Samples and Clinical Status Groups with Adjusted p-values", style = "heading 1") %>%
  body_add_table(result_df, style = "table_template")
print(doc_s1a, target = "Table_S1a.docx")

# Pairwise comparisons if significant
significant_vars <- names(kruskal_p_values)[kruskal_p_values < 0.05]
if (length(significant_vars) > 0) {
  for (variable in significant_vars) {
    pairwise_test <- pairwise.wilcox.test(data_clinical_status[[variable]], data_clinical_status$group, p.adjust.method = "none")
    pairwise_p_values <- as.data.frame(as.table(pairwise_test$p.value))
    pairwise_p_values$Variable <- variable
    pairwise_p_values_list[[variable]] <- pairwise_p_values
  }
  
  pairwise_all <- do.call(rbind, pairwise_p_values_list)
  pairwise_all <- pairwise_all %>%
    filter(!is.na(Freq)) %>%
    mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      Comparison = ifelse(Var1 < Var2, paste(Var1, "vs", Var2), paste(Var2, "vs", Var1))
    ) %>%
    select(Variable, Comparison, P_value = Freq)
  
  pairwise_all$Adjusted_P_value <- p.adjust(pairwise_all$P_value, method = "BH")
  pairwise_all$Adjusted_P_value <- sapply(pairwise_all$Adjusted_P_value, function(p) {
    if (p < 0.01) "<0.01" else sprintf("%.2f", p)
  })
  
  pairwise_df <- pairwise_all %>%
    select(Variable, Comparison, Adjusted_P_value) %>%
    pivot_wider(names_from = Comparison, values_from = Adjusted_P_value)
  
  if (exists("variable_labels")) {
    pairwise_df$Variable <- variable_labels[pairwise_df$Variable]
  }
  
  print(kable(pairwise_df, caption = "Table S1b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", align = 'c'))
  doc_s1b <- read_docx() %>%
    body_add_par("Table S1b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", style = "heading 1") %>%
    body_add_table(pairwise_df, style = "table_template")
  print(doc_s1b, target = "Table_S1b.docx")
}

# TABLE 3: Diagnostic Metrics for Significant Fragmentomics Variables ----

# Proceed only if there are significant variables from TABLE S1a & S1b
if (length(significant_vars) > 0) {
  # Initialize a list to store ROC metrics
  roc_metrics_list <- list()
  
  # Loop over each significant variable to generate ROC curves
  for (var in significant_vars) {
    # Extract the variable and group from the data
    variable_data <- data_clinical_status %>%
      select(all_of(c(var, "group"))) %>%
      drop_na()
    
    # Convert 'group' to a binary variable (1 = AWD, 0 = AWoD)
    variable_data$group_binary <- ifelse(variable_data$group == "AWD", 1, 0)
    
    # Calculate ROC curve
    roc_obj <- roc(response = variable_data$group_binary, predictor = variable_data[[var]], 
                   levels = c(0, 1), direction = "auto")
    
    # Calculate AUC and its 95% CI
    auc_value <- as.numeric(auc(roc_obj))
    auc_ci <- ci.auc(roc_obj)
    auc_ci_lower <- as.numeric(auc_ci[1])
    auc_ci_upper <- as.numeric(auc_ci[3])
    
    # Find optimal cutpoint using Youden's index
    optimal_coords <- coords(roc_obj, x = "best", best.method = "youden", 
                             ret = c("threshold", "sensitivity", "specificity", "youden"), transpose = FALSE)
    
    # Extract threshold, sensitivity, specificity, Youden index
    optimal_threshold <- as.numeric(optimal_coords["threshold"])
    sensitivity <- as.numeric(optimal_coords["sensitivity"])
    specificity <- as.numeric(optimal_coords["specificity"])
    youden_index <- as.numeric(optimal_coords["youden"])
    
    # Classify observations based on the optimal threshold
    variable_data$predicted <- ifelse(variable_data[[var]] >= optimal_threshold, 1, 0)
    
    # Create confusion matrix
    confusion <- table(Predicted = variable_data$predicted, Actual = variable_data$group_binary)
    
    # Compute PPV and NPV
    TP <- ifelse("1" %in% rownames(confusion) & "1" %in% colnames(confusion), confusion["1", "1"], 0)
    TN <- ifelse("0" %in% rownames(confusion) & "0" %in% colnames(confusion), confusion["0", "0"], 0)
    FP <- ifelse("1" %in% rownames(confusion) & "0" %in% colnames(confusion), confusion["1", "0"], 0)
    FN <- ifelse("0" %in% rownames(confusion) & "1" %in% colnames(confusion), confusion["0", "1"], 0)
    PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
    
    # Store the ROC metrics
    roc_metrics_list[[var]] <- data.frame(
      Variable = var,  # Keep original variable name
      AUC = auc_value,
      AUC_Lower_CI = auc_ci_lower,
      AUC_Upper_CI = auc_ci_upper,
      Optimal_Cutpoint = optimal_threshold,
      Sensitivity = sensitivity,
      Specificity = specificity,
      PPV = PPV,
      NPV = NPV,
      Youden_Index = youden_index,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine the ROC metrics into a data frame
  roc_metrics_df <- bind_rows(roc_metrics_list)
  
  # Create a new column for variable labels
  if (exists("variable_labels")) {
    roc_metrics_df$Variable_Label <- variable_labels[roc_metrics_df$Variable]
  } else {
    roc_metrics_df$Variable_Label <- roc_metrics_df$Variable
  }
  
  # Round numeric columns for readability
  roc_metrics_df <- roc_metrics_df %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)))
  
  # Create AUC with CI column
  roc_metrics_df <- roc_metrics_df %>%
    mutate(AUC_CI = paste0(AUC, " (", AUC_Lower_CI, "-", AUC_Upper_CI, ")"))
  
  # Reorder columns (include AUC for further analysis)
  roc_metrics_df <- roc_metrics_df %>%
    select(Variable, Variable_Label, AUC, AUC_CI, Optimal_Cutpoint, Sensitivity, Specificity, PPV, NPV, Youden_Index)
  
  # Print the ROC metrics table (use labels)
  print(kable(roc_metrics_df %>% select(-Variable, -AUC), caption = "Table 3: Diagnostic Metrics for Significant Fragmentomics Variables", align = "c"))
  
  # Save the table to a Word document
  doc_table3 <- read_docx() %>%
    body_add_par("Table 3: Diagnostic Metrics for Significant Fragmentomics Variables", style = "heading 1") %>%
    body_add_table(roc_metrics_df %>% select(-Variable, -AUC), style = "table_template")
  
  # Save the document
  print(doc_table3, target = "Table3.docx")
}


# TABLE 3b: Diagnostic Metrics for Significant Fragmentomics Variables (Filtered Dataset) ----

# Apply the specified filtering criteria
# Healthy group: disease_cat == "none", exclude if relapse_status == "yes" or lb_type == "unresect"
healthy_group <- highrisk_fragment_samples %>%
  filter(
    disease_cat == "none",
    relapse_status != "yes" | is.na(relapse_status),
    lb_type != "unresect"
  )

# Sick group: disease_cat == "regional_macro" or "metastasis"
sick_group <- highrisk_fragment_samples %>%
  filter(disease_cat %in% c("regional_macro", "metastasis"))

# Combine healthy and sick groups
data_filtered <- bind_rows(healthy_group, sick_group)

# Create grouping variable 'group' for diagnostic analysis
data_filtered <- data_filtered %>%
  mutate(group = case_when(
    disease_cat == "none" ~ "Healthy",
    disease_cat %in% c("regional_macro", "metastasis") ~ "Sick",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Proceed only if there are significant variables from TABLE S1a & S1b
if (length(significant_vars) > 0) {
  # Initialize a list to store ROC metrics
  roc_metrics_list_3b <- list()
  
  # Loop over each significant variable to generate ROC metrics
  for (var in significant_vars) {
    # Extract the variable and group from the filtered data
    variable_data <- data_filtered %>%
      select(all_of(c(var, "group"))) %>%
      drop_na()
    
    # Ensure we have at least two groups with data
    if (length(unique(variable_data$group)) < 2) {
      next  # Skip this variable
    }
    
    # Convert 'group' to a binary variable (1 = Sick, 0 = Healthy)
    variable_data$group_binary <- ifelse(variable_data$group == "Sick", 1, 0)
    
    # Calculate ROC curve
    roc_obj <- roc(response = variable_data$group_binary, predictor = variable_data[[var]], 
                   levels = c(0, 1), direction = "auto")
    
    # Calculate AUC and its 95% CI
    auc_value <- as.numeric(auc(roc_obj))
    auc_ci <- ci.auc(roc_obj)
    auc_ci_lower <- as.numeric(auc_ci[1])
    auc_ci_upper <- as.numeric(auc_ci[3])
    
    # Find optimal cutpoint using Youden's index
    optimal_coords <- coords(roc_obj, x = "best", best.method = "youden", 
                             ret = c("threshold", "sensitivity", "specificity", "youden"), transpose = TRUE)
    
    # Extract threshold, sensitivity, specificity, Youden index
    optimal_threshold <- as.numeric(optimal_coords["threshold"])
    sensitivity <- as.numeric(optimal_coords["sensitivity"])
    specificity <- as.numeric(optimal_coords["specificity"])
    youden_index <- as.numeric(optimal_coords["youden"])
    
    # Classify observations based on the optimal threshold
    variable_data$predicted <- ifelse(variable_data[[var]] >= optimal_threshold, 1, 0)
    
    # Create confusion matrix
    confusion <- table(Predicted = variable_data$predicted, Actual = variable_data$group_binary)
    
    # Compute PPV and NPV
    TP <- ifelse("1" %in% rownames(confusion) & "1" %in% colnames(confusion), confusion["1", "1"], 0)
    TN <- ifelse("0" %in% rownames(confusion) & "0" %in% colnames(confusion), confusion["0", "0"], 0)
    FP <- ifelse("1" %in% rownames(confusion) & "0" %in% colnames(confusion), confusion["1", "0"], 0)
    FN <- ifelse("0" %in% rownames(confusion) & "1" %in% colnames(confusion), confusion["0", "1"], 0)
    PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
    NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
    
    # Store the ROC metrics
    roc_metrics_list_3b[[var]] <- data.frame(
      Variable = var,
      AUC = auc_value,
      AUC_Lower_CI = auc_ci_lower,
      AUC_Upper_CI = auc_ci_upper,
      Optimal_Cutpoint = optimal_threshold,
      Sensitivity = sensitivity,
      Specificity = specificity,
      PPV = PPV,
      NPV = NPV,
      Youden_Index = youden_index,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine the ROC metrics into a data frame
  roc_metrics_df_3b <- bind_rows(roc_metrics_list_3b)
  
  # Replace variable names with labels if 'variable_labels' exists
  if (exists("variable_labels")) {
    roc_metrics_df_3b$Variable <- variable_labels[roc_metrics_df_3b$Variable]
  }
  
  # Round numeric columns for readability
  roc_metrics_df_3b <- roc_metrics_df_3b %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 3)))
  
  # Create AUC with CI column
  roc_metrics_df_3b <- roc_metrics_df_3b %>%
    mutate(AUC_CI = paste0(AUC, " (", AUC_Lower_CI, "-", AUC_Upper_CI, ")"))
  
  # Reorder columns (include AUC if needed)
  roc_metrics_df_3b <- roc_metrics_df_3b %>%
    select(Variable, AUC, AUC_CI, Optimal_Cutpoint, Sensitivity, Specificity, PPV, NPV, Youden_Index)
  
  # Print the ROC metrics table (exclude the AUC column)
  print(kable(roc_metrics_df_3b %>% select(-AUC), caption = "Table 3b: Diagnostic Metrics for Significant Fragmentomics Variables (Filtered Dataset)", align = "c"))
  
  # Save the table to a Word document
  doc_table3b <- read_docx() %>%
    body_add_par("Table 3b: Diagnostic Metrics for Significant Fragmentomics Variables (Filtered Dataset)", style = "heading 1") %>%
    body_add_table(roc_metrics_df_3b %>% select(-AUC), style = "table_template")
  
  # Save the document
  print(doc_table3b, target = "Table3b.docx")
}

# TABLE S2a & S2b ----

# Analyze fragmentomics variables by disease category

# Function to abbreviate group names
abbreviate_group <- function(group_name) {
  abbreviations <- c(
    "Control" = "Con",
    "No Evidence of Disease" = "NED",
    "Primary Tumor Only" = "Tum",
    "Regional Microscopic" = "Mic",
    "Regional Macroscopic" = "Mac",
    "Metastatic" = "Met"
  )
  return(abbreviations[group_name])
}

# Create grouping variable 'group' for disease category
data_disease_cat <- highrisk_fragment_samples %>%
  mutate(group = case_when(
    lb_type == "control" ~ "Control",
    disease_cat == "none" ~ "No Evidence of Disease",
    disease_cat == "tumor" ~ "Primary Tumor Only",
    disease_cat == "regional_micro" ~ "Regional Microscopic",
    disease_cat == "regional_macro" ~ "Regional Macroscopic",
    disease_cat == "metastasis" ~ "Metastatic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Set factor levels and calculate sample sizes
desired_order <- c("Control", "No Evidence of Disease", "Primary Tumor Only",
                   "Regional Microscopic", "Regional Macroscopic", "Metastatic")
data_disease_cat$group <- factor(data_disease_cat$group, levels = desired_order)
disease_cat_groups <- levels(data_disease_cat$group)
sample_sizes <- setNames(sapply(disease_cat_groups, function(dc) {
  sum(data_disease_cat$group == dc)
}), disease_cat_groups)
total_sample_size <- nrow(data_disease_cat)

# Initialize result matrix
result_matrix <- matrix(ncol = length(disease_cat_groups) + 1, nrow = length(fragmentomics_variables))
colnames(result_matrix) <- c(paste0("All samples (n=", total_sample_size, ")"),
                             paste0(disease_cat_groups, " (n=", sample_sizes[disease_cat_groups], ")"))
rownames(result_matrix) <- fragmentomics_variables

# Perform statistical analysis
kruskal_p_values <- numeric(length(fragmentomics_variables))
names(kruskal_p_values) <- fragmentomics_variables
pairwise_p_values_list <- list()

for (variable in fragmentomics_variables) {
  # Mean ± SD for all samples
  result_matrix[variable, paste0("All samples (n=", total_sample_size, ")")] <- calculate_mean_sd(data_disease_cat, variable)
  
  # Mean ± SD for each group
  for (dc in disease_cat_groups) {
    data_subset <- data_disease_cat[data_disease_cat$group == dc, ]
    result_matrix[variable, paste0(dc, " (n=", sample_sizes[dc], ")")] <- calculate_mean_sd(data_subset, variable)
  }
  
  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(as.formula(paste(variable, "~ group")), data = data_disease_cat)
  kruskal_p_values[variable] <- kruskal_test$p.value
}

# Adjust p-values and format
adjusted_p_values <- p.adjust(kruskal_p_values, method = "BH")
formatted_p_values <- sapply(adjusted_p_values, function(p) {
  if (p < 0.01) "<0.01" else sprintf("%.2f", p)
})
result_matrix <- cbind(result_matrix, "Adjusted p-value" = formatted_p_values)

# Convert to data frame and add variable labels
result_df <- as.data.frame(result_matrix)
if (exists("variable_labels")) {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable") %>%
    mutate(Variable = variable_labels[Variable])
} else {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable")
}

# Print and save the main summary table
print(kable(result_df, caption = "Table S2a: Mean ± SD for All Samples and Disease Category Groups with Adjusted p-values", align = 'c'))
doc_s2a <- read_docx() %>%
  body_add_par("Table S2a: Mean ± SD for All Samples and Disease Category Groups with Adjusted p-values", style = "heading 1") %>%
  body_add_table(result_df, style = "table_template")
print(doc_s2a, target = "Table_S2a.docx")

# Pairwise comparisons if significant
significant_vars <- names(kruskal_p_values)[kruskal_p_values < 0.05]
if (length(significant_vars) > 0) {
  for (variable in significant_vars) {
    pairwise_test <- pairwise.wilcox.test(data_disease_cat[[variable]], data_disease_cat$group, p.adjust.method = "none")
    pairwise_p_values <- as.data.frame(as.table(pairwise_test$p.value))
    pairwise_p_values$Variable <- variable
    pairwise_p_values_list[[variable]] <- pairwise_p_values
  }
  
  pairwise_all <- do.call(rbind, pairwise_p_values_list)
  pairwise_all <- pairwise_all %>%
    filter(!is.na(Freq)) %>%
    mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      Var1_abbr = abbreviate_group(Var1),
      Var2_abbr = abbreviate_group(Var2),
      Comparison = ifelse(Var1_abbr < Var2_abbr, paste(Var1_abbr, "vs", Var2_abbr), paste(Var2_abbr, "vs", Var1_abbr))
    ) %>%
    select(Variable, Comparison, P_value = Freq)
  
  pairwise_all$Adjusted_P_value <- p.adjust(pairwise_all$P_value, method = "BH")
  pairwise_all$Adjusted_P_value <- sapply(pairwise_all$Adjusted_P_value, function(p) {
    if (p < 0.01) "<0.01" else sprintf("%.2f", p)
  })
  
  pairwise_df <- pairwise_all %>%
    select(Variable, Comparison, Adjusted_P_value) %>%
    pivot_wider(names_from = Comparison, values_from = Adjusted_P_value)
  
  if (exists("variable_labels")) {
    pairwise_df$Variable <- variable_labels[pairwise_df$Variable]
  }
  
  print(kable(pairwise_df, caption = "Table S2b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", align = 'c'))
  doc_s2b <- read_docx() %>%
    body_add_par("Table S2b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", style = "heading 1") %>%
    body_add_table(pairwise_df, style = "table_template")
  print(doc_s2b, target = "Table_S2b.docx")
}

# TABLE S3a & S3b ----

# Analyze fragmentomics variables by lb_type

# Abbreviate group names function
abbreviate_group <- function(group_name) {
  abbreviations <- c(
    "Control" = "Con",
    "Preoperative" = "Pre",
    "Postoperative" = "Post",
    "Follow-up" = "Fol",
    "Relapse" = "Rel",
    "Unresectable" = "Unr"
  )
  sapply(group_name, function(g) abbreviations[g])
}

# Create grouping variable 'group' for lb_type
data_lb_type <- highrisk_fragment_samples %>%
  mutate(group = case_when(
    lb_type == "control" ~ "Control",
    lb_type == "preop" ~ "Preoperative",
    lb_type == "postop" ~ "Postoperative",
    lb_type == "followup" ~ "Follow-up",
    lb_type == "relapse" ~ "Relapse",
    lb_type == "unresect" ~ "Unresectable",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Set factor levels and calculate sample sizes
desired_order <- c("Control", "Preoperative", "Postoperative", "Follow-up", "Relapse", "Unresectable")
data_lb_type$group <- factor(data_lb_type$group, levels = desired_order)
lb_type_groups <- levels(data_lb_type$group)
sample_sizes <- setNames(sapply(lb_type_groups, function(lb) {
  sum(data_lb_type$group == lb)
}), lb_type_groups)
total_sample_size <- nrow(data_lb_type)

# Initialize result matrix
result_matrix <- matrix(ncol = length(lb_type_groups) + 1, nrow = length(fragmentomics_variables))
colnames(result_matrix) <- c(paste0("All samples (n=", total_sample_size, ")"),
                             paste0(lb_type_groups, " (n=", sample_sizes[lb_type_groups], ")"))
rownames(result_matrix) <- fragmentomics_variables

# Perform statistical analysis
kruskal_p_values <- numeric(length(fragmentomics_variables))
names(kruskal_p_values) <- fragmentomics_variables
pairwise_p_values_list <- list()

for (variable in fragmentomics_variables) {
  # Mean ± SD for all samples
  result_matrix[variable, paste0("All samples (n=", total_sample_size, ")")] <- calculate_mean_sd(data_lb_type, variable)
  
  # Mean ± SD for each group
  for (lb in lb_type_groups) {
    data_subset <- data_lb_type[data_lb_type$group == lb, ]
    result_matrix[variable, paste0(lb, " (n=", sample_sizes[lb], ")")] <- calculate_mean_sd(data_subset, variable)
  }
  
  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(as.formula(paste(variable, "~ group")), data = data_lb_type)
  kruskal_p_values[variable] <- kruskal_test$p.value
}

# Adjust p-values and format
adjusted_p_values <- p.adjust(kruskal_p_values, method = "BH")
formatted_p_values <- sapply(adjusted_p_values, function(p) {
  if (p < 0.01) "<0.01" else sprintf("%.2f", p)
})
result_matrix <- cbind(result_matrix, "Adjusted p-value" = formatted_p_values)

# Convert to data frame and add variable labels
result_df <- as.data.frame(result_matrix)
if (exists("variable_labels")) {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable") %>%
    mutate(Variable = variable_labels[Variable])
} else {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable")
}

# Print and save the main summary table
print(kable(result_df, caption = "Table S3a: Mean ± SD for All Samples and LB Type Groups with Adjusted p-values", align = 'c'))
doc_s3a <- read_docx() %>%
  body_add_par("Table S3a: Mean ± SD for All Samples and LB Type Groups with Adjusted p-values", style = "heading 1") %>%
  body_add_table(result_df, style = "table_template")
print(doc_s3a, target = "Table_S3a.docx")

# Pairwise comparisons if significant
significant_vars <- names(kruskal_p_values)[kruskal_p_values < 0.05]
if (length(significant_vars) > 0) {
  for (variable in significant_vars) {
    pairwise_test <- pairwise.wilcox.test(data_lb_type[[variable]], data_lb_type$group, p.adjust.method = "none")
    pairwise_p_values <- as.data.frame(as.table(pairwise_test$p.value))
    pairwise_p_values$Variable <- variable
    pairwise_p_values_list[[variable]] <- pairwise_p_values
  }
  
  pairwise_all <- do.call(rbind, pairwise_p_values_list)
  pairwise_all <- pairwise_all %>%
    filter(!is.na(Freq)) %>%
    mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      Var1_abbr = abbreviate_group(Var1),
      Var2_abbr = abbreviate_group(Var2),
      Comparison = ifelse(Var1_abbr < Var2_abbr, paste(Var1_abbr, "vs", Var2_abbr), paste(Var2_abbr, "vs", Var1_abbr))
    ) %>%
    select(Variable, Comparison, P_value = Freq)
  
  pairwise_all$Adjusted_P_value <- p.adjust(pairwise_all$P_value, method = "BH")
  pairwise_all$Adjusted_P_value <- sapply(pairwise_all$Adjusted_P_value, function(p) {
    if (p < 0.01) "<0.01" else sprintf("%.2f", p)
  })
  
  pairwise_df <- pairwise_all %>%
    select(Variable, Comparison, Adjusted_P_value) %>%
    pivot_wider(names_from = Comparison, values_from = Adjusted_P_value)
  
  if (exists("variable_labels")) {
    pairwise_df$Variable <- variable_labels[pairwise_df$Variable]
  }
  
  print(kable(pairwise_df, caption = "Table S3b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", align = 'c'))
  doc_s3b <- read_docx() %>%
    body_add_par("Table S3b: Pairwise Wilcoxon Test Results for Significant Variables (BH Adjusted)", style = "heading 1") %>%
    body_add_table(pairwise_df, style = "table_template")
  print(doc_s3b, target = "Table_S3b.docx")
}

# TABLE S4 ----

# Analyze fragmentomics variables by disease activity in unresectable patients

# Filter the dataset
unresectable_samples <- highrisk_fragment_samples %>%
  filter(lb_type == "unresect" & disease_activity != "AAC")

# Ensure disease_activity is a factor
unresectable_samples$disease_activity <- as.factor(unresectable_samples$disease_activity)

# Calculate sample sizes
disease_activities <- levels(unresectable_samples$disease_activity)
sample_sizes <- setNames(sapply(disease_activities, function(activity) {
  sum(unresectable_samples$disease_activity == activity)
}), disease_activities)

# Initialize result matrix
result_matrix <- matrix(ncol = length(disease_activities), nrow = length(fragmentomics_variables))
colnames(result_matrix) <- paste0(disease_activities, " (n=", sample_sizes[disease_activities], ")")
rownames(result_matrix) <- fragmentomics_variables

# Perform statistical analysis
kruskal_p_values <- numeric(length(fragmentomics_variables))
names(kruskal_p_values) <- fragmentomics_variables

for (variable in fragmentomics_variables) {
  # Mean ± SD for each disease activity
  for (activity in disease_activities) {
    data_subset <- unresectable_samples[unresectable_samples$disease_activity == activity, ]
    if (nrow(data_subset) > 0) {
      result_matrix[variable, paste0(activity, " (n=", sample_sizes[activity], ")")] <- calculate_mean_sd(data_subset, variable)
    } else {
      result_matrix[variable, paste0(activity, " (n=", sample_sizes[activity], ")")] <- "NA ± NA"
    }
  }
  
  # Kruskal-Wallis test
  if (length(disease_activities) > 1) {
    kruskal_test <- kruskal.test(as.formula(paste(variable, "~ disease_activity")), data = unresectable_samples)
    kruskal_p_values[variable] <- kruskal_test$p.value
  } else {
    kruskal_p_values[variable] <- NA
  }
}

# Adjust p-values and format
adjusted_p_values <- p.adjust(kruskal_p_values, method = "BH")
formatted_p_values <- sapply(adjusted_p_values, function(p) {
  if (is.na(p)) "NA" else if (p < 0.01) "<0.01" else sprintf("%.2f", p)
})
result_matrix <- cbind(result_matrix, "Adjusted p-value (BH)" = formatted_p_values)

# Convert to data frame and add variable labels
result_df <- as.data.frame(result_matrix)
if (exists("variable_labels")) {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable") %>%
    mutate(Variable = variable_labels[Variable])
} else {
  result_df <- result_df %>%
    rownames_to_column(var = "Variable")
}

# Print and save the table
print(kable(result_df, caption = "Table S4: Mean ± SD for Each Variable by Disease Activity in Unresectable Patients with Adjusted p-values (BH)", align = 'c'))
doc_s4 <- read_docx() %>%
  body_add_par("Table S4: Mean ± SD for Each Variable by Disease Activity in Unresectable Patients with Adjusted p-values (BH)", style = "heading 1") %>%
  body_add_table(result_df, style = "table_template")
print(doc_s4, target = "Table_S4.docx")


# FIGURE 2: Correlation Matrix for Fragmentomics Variables ----

# Prepare the data for correlation analysis
df_cor <- highrisk_fragment_samples %>%
  select(all_of(fragmentomics_variables)) %>%
  mutate(across(everything(), as.numeric))  # Ensure all variables are numeric

# Rename variables using labels if available
if (exists("variable_labels")) {
  df_cor <- df_cor %>%
    rename_with(~ variable_labels[.x], everything())
}

# Compute the Pearson correlation matrix
cor_matrix <- cor(df_cor, use = "pairwise.complete.obs", method = "pearson")

# Generate and save the correlation plot using corrplot
jpeg("Figure2_Correlation_Matrix_Fragmentomics.jpeg", width = 12, height = 12, units = "in", res = 300)
corrplot(
  cor_matrix, 
  method = "circle", 
  type = "upper",
  addCoef.col = "black", 
  tl.col = "black", 
  tl.srt = 45,
  number.cex = 0.7, 
  col <- colorRampPalette(c("#3B4CC0", "#F5F5F5", "#B40426"))(200),
  title = "Figure 2: Correlation Matrix for Fragmentomics Variables", 
  mar = c(0,0,1,0)
)
dev.off()

# TABLE 2: GLMM Results for Associations Between Fragmentomics Variables and Clinical Status ----

# Ensure 'nhc' (patient ID) is a factor
highrisk_fragment_samples$nhc <- factor(highrisk_fragment_samples$nhc)

# Create a grouping variable 'group' for clinical status
data_clinical_status <- highrisk_fragment_samples %>%
  mutate(group = case_when(
    clinical_status == "AWoD" ~ "AWoD",
    clinical_status == "AWD" ~ "AWD",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# Convert 'group' to a factor with levels 'AWoD' and 'AWD'
data_clinical_status$group <- factor(data_clinical_status$group, levels = c("AWoD", "AWD"))

# Scale the fragmentomics variables
data_scaled <- data_clinical_status %>%
  mutate(across(all_of(fragmentomics_variables), ~ scale(.) %>% as.vector()))

# Remove rows with missing values in fragmentomics variables (if any)
data_scaled <- data_scaled %>%
  drop_na(all_of(fragmentomics_variables))

# Initialize a list to store results
results_list <- list()

'''
# Loop over each fragmentomics variable
for (var in fragmentomics_variables) {
  # Define the formula
  formula <- as.formula(paste(var, "~ group + (1 | nhc)"))
  
  # Fit the mixed-effects model
  model <- lmer(formula, data = data_scaled, REML = FALSE)
  
  # Tidy the model output
  model_tidy <- tidy(model)
  
  # Extract the coefficient for 'groupAWD'
  coef_info <- model_tidy %>%
    filter(term == "groupAWD") %>%
    select(estimate, std.error, statistic, p.value)
  
  # Store the results
  results_list[[var]] <- data.frame(
    Variable = var,
    Estimate = coef_info$estimate,
    Std_Error = coef_info$std.error,
    t_value = coef_info$statistic,
    p_value = coef_info$p.value
  )
}

# Combine results into a data frame
results_df <- bind_rows(results_list)

# Adjust p-values for multiple comparisons (Benjamini-Hochberg)
results_df <- results_df %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

# Add significance labels
results_df <- results_df %>%
  mutate(Significance = case_when(
    adj_p_value < 0.001 ~ "***",
    adj_p_value < 0.01 ~ "**",
    adj_p_value < 0.05 ~ "*",
    adj_p_value < 0.1 ~ ".",
    TRUE ~ ""
  ))

# Sort variables by adjusted p-value
results_df <- results_df %>%
  arrange(adj_p_value)

# Replace variable names with labels
results_df <- results_df %>%
  mutate(Variable_Label = variable_labels[Variable])

# **Round numeric columns to 2 decimals**
results_df <- results_df %>%
  mutate(across(c(Estimate, Std_Error, t_value), ~ round(.x, digits = 2)))

# Reorder columns for clarity
results_df <- results_df %>%
  select(Variable = Variable_Label, Estimate, Std_Error, t_value, adj_p_value, Significance)

# Format adjusted p-values
results_df <- results_df %>%
  mutate(adj_p_value = ifelse(adj_p_value < 0.001, "<0.001", sprintf("%.3f", adj_p_value)))

# Rename columns
results_df <- results_df %>%
  rename(
    "Adjusted p-value" = "adj_p_value",
    "Signif." = "Significance",
    "t-value" = "t_value",
    "Std. Error" = "Std_Error"
  )

# Print the results table
print(kable(results_df, caption = "Table 2: GLMM Results for Associations Between Fragmentomics Variables and Clinical Status", align = 'c'))

# Save the results table to a Word document
doc_table2 <- read_docx() %>%
  body_add_par("Table 2: GLMM Results for Associations Between Fragmentomics Variables and Clinical Status", style = "heading 1") %>%
  body_add_table(results_df, style = "table_template")

# Save the document
print(doc_table2, target = "Table2.docx")

# FIGURE S1: Forest Plot of Statistically Significant Fragmentomics Variables in GLMM Clinical Status Analysis ----

# Filter significant variables (adjusted p-value < 0.05)
significant_results <- results_df %>%
  mutate(adj_p_value_numeric = as.numeric(gsub("<0.001", "0.000999", `Adjusted p-value`))) %>%
  filter(adj_p_value_numeric < 0.05)

# Check if any variables are significant
if (nrow(significant_results) > 0) {
  # Prepare data for plotting
  significant_results <- significant_results %>%
    mutate(
      Direction = ifelse(Estimate > 0, "Higher in AWD", "Higher in AWoD"),
      Lower_CI = Estimate - 1.96 * `Std. Error`,
      Upper_CI = Estimate + 1.96 * `Std. Error`,
      Estimate_CI_Label = sprintf("%.2f [%.2f, %.2f]", Estimate, Lower_CI, Upper_CI),
      Variable = factor(Variable, levels = rev(unique(Variable)))
    )
  
  # Create the forest plot
  forest_plot <- ggplot(significant_results, aes(x = Estimate, y = Variable)) +
    geom_point(aes(color = Direction), size = 4) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI, color = Direction),
                   height = 0.3, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
    scale_color_manual(values = c("Higher in AWD" = "red", "Higher in AWoD" = "blue")) +
    labs(
      title = "Figure S1: Forest Plot of Statistically Significant Fragmentomics Variables",
      x = "Effect Size (Estimate)",
      y = "",
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    geom_text(
      aes(label = Estimate_CI_Label),
      hjust = 0.5,
      vjust = -1.5,
      size = 4,
      fontface = "bold",
      color = "black"
    ) +
    coord_cartesian(xlim = c(min(significant_results$Lower_CI) - 0.5, max(significant_results$Upper_CI) + 0.5))
  
  # Save the plot as JPEG
  ggsave(
    filename = "Figure_S1_Forest_Plot_GLMM_Clinical_Status.jpeg",
    plot = forest_plot,
    width = 12,
    height = 8,
    dpi = 300,
    units = "in"
  )
  
} else {
  print("No variables are significant at the specified threshold.")
}
'''
# Data modeling for diagnosis: Outcome-driven selection ----

# Identify Variables to Reduce

# Extract AUC values for significant variables
auc_values <- roc_metrics_df %>%
  select(Variable, AUC)

# Prepare data for correlation analysis
df_cor_significant <- highrisk_fragment_samples %>%
  select(all_of(significant_vars)) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix
cor_matrix <- cor(df_cor_significant, use = "pairwise.complete.obs", method = "pearson")

# Identify pairs with correlation > 0.8
cor_matrix_upper <- cor_matrix
cor_matrix_upper[lower.tri(cor_matrix_upper, diag = TRUE)] <- NA
cor_df <- melt(cor_matrix_upper, na.rm = TRUE)
high_cor_pairs <- cor_df %>%
  filter(abs(value) > 0.8) %>%
  arrange(desc(abs(value)))

# Exclude highly correlated variables
auc_values_vector <- setNames(auc_values$AUC, auc_values$Variable)
variables_to_keep <- significant_vars

for(i in seq_len(nrow(high_cor_pairs))) {
  var1 <- as.character(high_cor_pairs$Var1[i])
  var2 <- as.character(high_cor_pairs$Var2[i])
  
  if(!(var1 %in% variables_to_keep) || !(var2 %in% variables_to_keep)) {
    next
  }
  
  auc1 <- auc_values_vector[var1]
  auc2 <- auc_values_vector[var2]
  
  if (is.na(auc1) || is.na(auc2)) {
    next
  }
  
  if(auc1 >= auc2) {
    variables_to_keep <- setdiff(variables_to_keep, var2)
  } else {
    variables_to_keep <- setdiff(variables_to_keep, var1)
  }
}

# Final list of reduced variables
reduced_variables <- variables_to_keep

# Replace variable names with labels if available
if (exists("variable_labels")) {
  reduced_variable_labels <- variable_labels[reduced_variables]
  cat("Reduced list of variables after variable reduction:\n")
  print(reduced_variable_labels)
}

# Correlation Analysis Between Reduced Variables

# Extract reduced variables
df_reduced_vars <- highrisk_fragment_samples %>%
  select(all_of(reduced_variables)) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix
cor_matrix_reduced <- cor(df_reduced_vars, use = "pairwise.complete.obs", method = "pearson")

# Transform and clean the correlation matrix
cor_df_reduced <- melt(cor_matrix_reduced)
cor_df_reduced <- cor_df_reduced %>%
  mutate(
    Var1 = as.character(Var1),
    Var2 = as.character(Var2)
  ) %>%
  filter(Var1 != Var2) %>%
  mutate(
    pair = paste0(pmin(Var1, Var2), " vs ", pmax(Var1, Var2))
  ) %>%
  group_by(pair) %>%
  summarise(
    Var1 = first(Var1),
    Var2 = first(Var2),
    correlation = first(value),
    .groups = 'drop'
  )

# Sort correlations
cor_df_reduced_sorted <- cor_df_reduced %>%
  arrange(desc(abs(correlation)))

# Replace variable names with labels if available
if (exists("variable_labels")) {
  cor_df_reduced_sorted <- cor_df_reduced_sorted %>%
    mutate(
      Var1_Label = variable_labels[Var1],
      Var2_Label = variable_labels[Var2]
    )
} else {
  cor_df_reduced_sorted <- cor_df_reduced_sorted %>%
    mutate(
      Var1_Label = Var1,
      Var2_Label = Var2
    )
}

# Prepare final correlation table
cor_df_reduced_sorted <- cor_df_reduced_sorted %>%
  select(
    "Variable 1" = Var1_Label,
    "Variable 2" = Var2_Label,
    "Correlation Coefficient" = correlation
  ) %>%
  mutate(`Correlation Coefficient` = round(`Correlation Coefficient`, 3))

# Display the correlation analysis table
print(kable(cor_df_reduced_sorted, caption = "Correlation Analysis Table Between Reduced Variables", align = 'c'))

# Random Forest Model for Diagnosing AWD vs. AWoD

# Prepare the data for modeling

# Filter the dataset to include only samples with clinical_status 'AWD' or 'AWoD'
data_model <- highrisk_fragment_samples %>%
  filter(clinical_status %in% c("AWD", "AWoD")) %>%
  select(all_of(reduced_variables), nhc, clinical_status)

# Ensure predictors are numeric
data_model <- data_model %>%
  mutate(across(all_of(reduced_variables), as.numeric))

# Remove rows with missing values
data_model <- data_model %>%
  drop_na()

# Convert clinical_status to a factor with levels AWoD and AWD
data_model$clinical_status <- factor(data_model$clinical_status, levels = c("AWoD", "AWD"))

# Split the data into training and testing sets
set.seed(123)  # For reproducibility

# Stratify by patient ID to prevent data leakage
unique_nhc <- unique(data_model$nhc)
train_nhc <- sample(unique_nhc, size = floor(0.7 * length(unique_nhc)))
test_nhc <- setdiff(unique_nhc, train_nhc)

# Create training and testing datasets
train_data <- data_model %>% filter(nhc %in% train_nhc) %>% select(-nhc)
test_data <- data_model %>% filter(nhc %in% test_nhc) %>% select(-nhc)

# Train the Random Forest model
set.seed(123)
rf_model <- randomForest(
  clinical_status ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 500
)

# Evaluate model performance on the test set

# Predict probabilities for the test set
test_probabilities <- predict(rf_model, test_data, type = "prob")

# Compute ROC curve and AUC
roc_obj <- roc(test_data$clinical_status, test_probabilities[, "AWD"], levels = c("AWoD", "AWD"))
auc_value <- auc(roc_obj)

# Plot ROC curve
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))

# Generate confusion matrix
test_predictions <- predict(rf_model, test_data, type = "class")
confusion <- confusionMatrix(test_predictions, test_data$clinical_status)
print(confusion)

# Extract variable importance
variable_importance <- importance(rf_model)
var_imp_df <- data.frame(
  Variable = rownames(variable_importance),
  MeanDecreaseGini = variable_importance[, "MeanDecreaseGini"]
) %>%
  arrange(desc(MeanDecreaseGini))

# Replace variable names with labels if available
if (exists("variable_labels")) {
  var_imp_df$Variable_Label <- variable_labels[var_imp_df$Variable]
} else {
  var_imp_df$Variable_Label <- var_imp_df$Variable
}

# Display variable importance table
print(kable(var_imp_df %>% select(Variable_Label, MeanDecreaseGini), caption = "Variable Importance in Random Forest Model", align = 'c'))

# Cross-validation using caret package

# Set up training control with repeated cross-validation
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Train the model using caret
set.seed(123)
rf_caret_model <- train(
  clinical_status ~ .,
  data = train_data,
  method = "rf",
  metric = "ROC",
  trControl = train_control,
  importance = TRUE,
  ntree = 500
)

# Evaluate performance on the test set
test_probabilities_caret <- predict(rf_caret_model, test_data, type = "prob")
roc_obj_caret <- roc(test_data$clinical_status, test_probabilities_caret[, "AWD"], levels = c("AWoD", "AWD"))
auc_value_caret <- auc(roc_obj_caret)
plot(roc_obj_caret, main = paste("ROC Curve (AUC =", round(auc_value_caret, 3), ")"))

# Generate confusion matrix
test_predictions_caret <- predict(rf_caret_model, test_data)
confusion_caret <- confusionMatrix(test_predictions_caret, test_data$clinical_status)
print(confusion_caret)

# Extract variable importance from caret model
var_imp_caret <- varImp(rf_caret_model, scale = FALSE)
var_imp_caret_df <- var_imp_caret$importance %>%
  rownames_to_column(var = "Variable") %>%
  mutate(
    Overall = rowMeans(select(., -Variable))
  ) %>%
  arrange(desc(Overall))

# Replace variable names with labels if available
if (exists("variable_labels")) {
  var_imp_caret_df$Variable_Label <- variable_labels[var_imp_caret_df$Variable]
} else {
  var_imp_caret_df$Variable_Label <- var_imp_caret_df$Variable
}

# Display variable importance table from caret model
print(kable(var_imp_caret_df %>% select(Variable_Label, Overall), caption = "Variable Importance in Random Forest Model (Caret)", align = 'c'))

# Data modeling for diagnosis: Lasso-regression reduction technique ----

# Prepare the Data
data_model <- highrisk_fragment_samples %>%
  filter(clinical_status %in% c("AWD", "AWoD")) %>%
  select(all_of(fragmentomics_variables), nhc, clinical_status) %>%
  mutate(across(all_of(fragmentomics_variables), as.numeric)) %>%
  drop_na() %>%
  mutate(
    clinical_status = ifelse(clinical_status == "AWD", 1, 0)
  )

data_model_vars <- data_model %>%
  select(-nhc)

# Variable Selection Using LASSO Regression
x <- as.matrix(data_model_vars %>% select(-clinical_status))
y <- data_model_vars$clinical_status

# Ensure valid column names
colnames(x) <- make.names(colnames(x), unique = TRUE)

set.seed(123)
lasso_cv <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha = 1,
  nfolds = 5,
  type.measure = "class"
)

# Use lambda_min to select variables
lambda_min <- lasso_cv$lambda.min

lasso_coef <- coef(lasso_cv, s = lambda_min)
lasso_coef_matrix <- as.matrix(lasso_coef)

nonzero_coef_indices <- which(lasso_coef_matrix != 0)
selected_variables <- rownames(lasso_coef_matrix)[nonzero_coef_indices]
selected_variables <- selected_variables[selected_variables != "(Intercept)"]

cat("Variables selected by LASSO regression:\n")
print(selected_variables)

# Prepare data with selected variables
data_selected <- data_model %>%
  select(all_of(selected_variables), nhc, clinical_status)

data_selected$clinical_status <- factor(data_selected$clinical_status, levels = c(0, 1), labels = c("AWoD", "AWD"))

# Split the data into training and testing sets, stratified by 'nhc'
set.seed(123)
unique_nhc <- unique(data_selected$nhc)
train_nhc <- sample(unique_nhc, size = floor(0.7 * length(unique_nhc)))
test_nhc <- setdiff(unique_nhc, train_nhc)

train_data <- data_selected %>% filter(nhc %in% train_nhc) %>% select(-nhc)
test_data <- data_selected %>% filter(nhc %in% test_nhc) %>% select(-nhc)

# Ensure valid column names
colnames(train_data) <- make.names(colnames(train_data), unique = TRUE)
colnames(test_data) <- make.names(colnames(test_data), unique = TRUE)

# Set up training control
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = TRUE,
  verboseIter = FALSE
)

# Train a logistic regression model using caret
set.seed(123)
model_logistic <- train(
  clinical_status ~ .,
  data = train_data,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = train_control
)

# Evaluate model performance on the test set
test_probabilities <- predict(model_logistic, test_data, type = "prob")

# Compute ROC curve and AUC
roc_obj <- roc(test_data$clinical_status, test_probabilities[, "AWD"], levels = c("AWoD", "AWD"))
auc_value <- auc(roc_obj)

# Determine optimal threshold using Youden's Index
optimal_coords <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
optimal_threshold <- optimal_coords["threshold"][[1]]
cat("Optimal threshold selected (Youden's Index):", optimal_threshold, "\n")

# Generate final predictions at the optimal threshold
test_predictions_optimal <- ifelse(test_probabilities[, "AWD"] >= optimal_threshold, "AWD", "AWoD")
test_predictions_optimal <- factor(test_predictions_optimal, levels = c("AWoD", "AWD"))

# Generate confusion matrix at the optimal threshold
cm_optimal <- confusionMatrix(test_predictions_optimal, test_data$clinical_status, positive = "AWD")
print(cm_optimal)

# Plot ROC curve
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))

# Experiment with Variable Combinations
evaluate_model <- function(variables) {
  data_subset <- data_model %>%
    select(all_of(variables), nhc, clinical_status)
  
  data_subset$clinical_status <- factor(data_subset$clinical_status, levels = c(0, 1), labels = c("AWoD", "AWD"))
  
  train_data <- data_subset %>% filter(nhc %in% train_nhc) %>% select(-nhc)
  test_data <- data_subset %>% filter(nhc %in% test_nhc) %>% select(-nhc)
  
  colnames(train_data) <- make.names(colnames(train_data), unique = TRUE)
  colnames(test_data) <- make.names(colnames(test_data), unique = TRUE)
  
  set.seed(123)
  model <- train(
    clinical_status ~ .,
    data = train_data,
    method = "glm",
    family = "binomial",
    metric = "ROC",
    trControl = train_control
  )
  
  test_probabilities <- predict(model, test_data, type = "prob")
  test_predictions <- ifelse(test_probabilities[, "AWD"] >= optimal_threshold, "AWD", "AWoD")
  test_predictions <- factor(test_predictions, levels = c("AWoD", "AWD"))
  
  cm <- confusionMatrix(test_predictions, test_data$clinical_status, positive = "AWD")
  
  return(list(
    Model = model,
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Accuracy = cm$overall["Accuracy"]
  ))
}

# Example: Try removing one variable at a time
variable_combinations <- combn(selected_variables, length(selected_variables) - 1, simplify = FALSE)
performance_results <- data.frame()

for (variables in variable_combinations) {
  model_result <- evaluate_model(variables)
  performance_results <- rbind(performance_results, data.frame(
    Variables = paste(variables, collapse = ", "),
    Sensitivity = model_result$Sensitivity,
    Specificity = model_result$Specificity,
    Accuracy = model_result$Accuracy
  ))
}

print(kable(performance_results, caption = "Model Performance with Different Variable Combinations", align = 'c'))

best_model_index <- which.max(performance_results$Sensitivity)
best_variables <- variable_combinations[[best_model_index]]
cat("Best variable combination based on sensitivity:\n")
print(best_variables)

# Train final model with the best variable combination

final_model_result <- evaluate_model(best_variables)
final_model <- final_model_result$Model
final_sensitivity <- final_model_result$Sensitivity
final_specificity <- final_model_result$Specificity
final_accuracy <- final_model_result$Accuracy

cat("\nFinal model performance:\n")
cat("Sensitivity:", final_sensitivity, "\n")
cat("Specificity:", final_specificity, "\n")
cat("Accuracy:", final_accuracy, "\n")

