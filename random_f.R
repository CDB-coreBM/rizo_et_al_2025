library(randomForest)
library(ranger)
library(ROCR)
library(tidyverse)
library(pROC)

set.seed(42)
# Filter the dataset to include only samples with clinical_status 'AWD' or 'AWoD'
data_model <- highrisk_fragment_samples %>%
  filter(clinical_status %in% c("AWD", "AWoD")) %>%
  select(all_of(fragmentomics_variables), clinical_status, ldh, s100, sex, age)

# names of features
features <- setdiff(names(data_model), "clinical_status")

library(purrr)
# Balance data by age
data_balanced <- data_model %>%
  group_by(clinical_status, age_group = cut(age, breaks = seq(0, 100, by = 5))) %>%
  group_split() %>%
  map_df(~ if (nrow(.x) >= min(table(data_model$clinical_status))) {
    sample_n(.x, min(table(data_model$clinical_status)))
  } else {
    .x # Keep smaller groups as is
  }) %>%
  ungroup()

# Convert clinical_status to a factor with levels AWoD and AWD
data_balanced$clinical_status <- factor(data_balanced$clinical_status, levels = c("AWoD", "AWD"))

data_balanced <- data_balanced %>% select(-age) %>% filter(!is.na(ldh)) %>% filter(!is.na(s100))
mtry <- tuneRF(select(data_balanced, -clinical_status), data_balanced$clinical_status, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)

best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry)
print(best.m)

# Now use optimal mtry (the one with minimum bag(OOB) error.) in new model:
rf <-randomForest(clinical_status~., data=data_balanced, mtry=best.m, importance=TRUE,ntree=500)

#Conditional=True, adjusts for correlations between predictors.
i_scores <- caret::varImp(rf, conditional = TRUE) 

#Gathering rownames in 'var'  and converting it to the factor
#to provide 'fill' parameter for the bar chart. 
i_scores <- i_scores %>% tibble::rownames_to_column("var") 
i_scores$var<- i_scores$var %>% as.factor()
i_scores <- arrange(i_scores, desc(AWoD))
i_scores_plot <- i_scores[1:20,]

i_bar <- i_scores_plot %>% ggplot(aes(x = AWoD, y= reorder(var, AWoD), fill = var)) + 
  geom_bar(stat = "identity", show.legend = FALSE, width = 1) + 
  labs(x = NULL, y = NULL, title = "Feature importance for shannon_entropy")+ 
  theme_minimal()
plot(i_bar)

# Export dataset here.
write.table(data_balanced, sep=',', file='~/data_balanced_not_imputed.csv')



library(ISLR)
library(ROCR)

# Load a binary classification dataset from ISLR package
mydata <- ISLR::Default

# Set seed
set.seed(1234)

# 70% of dataset goes to training data and remaining 30% to test data
train_idx  <- sample(c(TRUE, FALSE), nrow(data_balanced), replace=TRUE, prob=c(0.8,0.2))
train <- data_balanced[train_idx, ]
test <- data_balanced[!train_idx, ]

# Build logistic regression model
model <- glm(clinical_status~ldh+s100, family="binomial", data=train)

# Calculate predicted probability of default of test data
predicted <- predict(model, test, type="response")

# Storing Model Performance Scores
pred  <- prediction(predicted, test$clinical_status)

# Calculating Area under Curve
perf <- performance(pred,"auc")
auc <- as.numeric(perf@y.values)
auc

# Plot ROC curve
roc_curve <- performance(pred, "tpr", "fpr")
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)
abline(0, 1, col = "gray", lty = 2, lwd = 1)
text(0.5, 0.3, paste("AUC =", round(auc, 2)), adj = c(0.5, 0.5), col = "black", cex = 1.5)
