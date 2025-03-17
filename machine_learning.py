# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, confusion_matrix
from imblearn.over_sampling import SMOTE
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

from statistics import mean
from matplotlib import pyplot
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import ConfusionMatrixDisplay


# Load the dataset
#data = pd.read_csv("../data/data_balanced.csv", sep=',')

data = pd.read_csv("../data/data_balanced_not_imputed.csv", sep=',')

# Label encode 'sex'
data["sex"] = LabelEncoder().fit_transform(data["sex"])
# Label encode 'age_group'
data["age_group"] = LabelEncoder().fit_transform(data["age_group"])

# Check data structure
print(data.info())
print(data.head())

# Split the features and target
X = data.drop("clinical_status", axis=1)
y = data["clinical_status"]

# Drop variables
X = X.drop("avg_size", axis=1)
X = X.drop("prop_fragment2", axis=1)
X = X.drop("prop_fragment3", axis=1)
#X = X.drop("avg_sz_r20_150", axis=1)
#X = X.drop("avg_sz_r160_180", axis=1)
#X = X.drop("avg_sz_r180_220", axis=1)
#X = X.drop("avg_sz_r250_320", axis=1)
#X = X.drop("avg_sz_r320_700", axis=1)
#X = X.drop("avg_sz_r700_1300", axis=1)
#X = X.drop("pct_cv_r20_150", axis=1)
#X = X.drop("pct_cv_r160_180", axis=1)
#X = X.drop("pct_cv_r180_220", axis=1)
#X = X.drop("pct_cv_r250_320", axis=1)
#X = X.drop("pct_cv_r320_700", axis=1)
#X = X.drop("pct_cv_r700_1300", axis=1)
#X = X.drop("pct_tot_r160_180", axis=1)
#X = X.drop("ldh", axis=1)
#X = X.drop("s100", axis=1)

# Ensure the target variable is binary and encoded (e.g., 0 for "AWoD", 1 for "AWD")
y = y.map({"AWoD": 0, "AWD": 1})  # Adjust as necessary

# Check clinical_status distribution
print("Clinical status:\n", y.value_counts())


# Randomly split dataset to test and train set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, shuffle=True, random_state=42)
print(y_test.value_counts())

# Change the scale of critical variables
X_train['s100'] = [np.log(val) if not val <= 0 else val for val in X_train['s100']]
X_test['s100'] = [np.log(val) if not val <= 0 else val for val in X_test['s100']]

# for now, I will keep the X_train and y_train together
X_train['clinical_status'] = y_train

'''
# Now let's check the distribution of each variable regarding the class.
for c in X_train.columns[:-1]:
    sns.catplot(x='clinical_status', y=c, data=X_train, kind='violin')
    plt.show()
'''

# Correlation matrix    
plt.figure(figsize=(25,16))
sns.heatmap(X_train.corr(), annot=True)
plt.show()

# Show features with high correlation
X_train.corr()['clinical_status'][((X_train.corr()['clinical_status'] > 0.15 ) | (X_train.corr()['clinical_status'] < -0.15))]

'''
# first, shuffle the records
X_train = X_train.sample(frac=1)
disease = X_train.loc[X_train['clinical_status'] == 1]
n_disease = disease.shape[0]
print("Number of disease in train data: " + str(n_disease))

# now, let's select the non-fraud records and keep the same number as fraud ones
non_disease = X_train.loc[X_train['clinical_status'] == 0][:n_disease]
print("Number of disease in new data: " + str(non_disease.shape[0]))

# let's merge the two datasets together and shuffle
sub_X_train = pd.concat([disease, non_disease])
sub_X_train = sub_X_train.sample(frac=1)
print(sub_X_train.head())

print(sub_X_train.corr()['clinical_status'][((sub_X_train.corr()['clinical_status'] >= 0.3 ) | (sub_X_train.corr()['clinical_status'] < -0.3))])
most_correlated = sub_X_train.corr()['clinical_status'][((sub_X_train.corr()['clinical_status'] >= 0.3 ) | (sub_X_train.corr()['clinical_status'] < -0.3))].index

for c in most_correlated:
    sns.catplot(x="clinical_status", y=c, data=sub_X_train, kind='box')
    plt.show()
'''

from sklearn.linear_model import LogisticRegressionCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from xgboost import XGBClassifier
from sklearn.linear_model import LinearRegression
import tensorflow as tf
from tensorflow.keras import layers, optimizers, Sequential
from scikeras.wrappers import KerasClassifier, KerasRegressor


'''
def KerasMLP(hidden=1000, activation='relu', learning_rate=0.01):
    model = Sequential([
        layers.Dense(hidden, activation=activation),
        layers.Dense(1, activation='sigmoid')
    ])
    opt = optimizers.Adam(learning_rate=learning_rate)
    model.compile(
        optimizer=opt,
        loss='binary_crossentropy',
        metrics=["accuracy"]
    )
    return model
 
classifiers = {
    "KNearest": KNeighborsClassifier(),
    "Support Vector Classifier": SVC(),
    "GradBoost": GradientBoostingClassifier(),
    "XGBCBoost" : XGBClassifier(),
    "RandomForestClassifier": RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=6),
}

# prepare data
sub_y_train = sub_X_train['clinical_status']
sub_X_train.drop('clinical_status', axis=1, inplace=True)

from sklearn.model_selection import cross_val_score

cv_scores_mean = []
cv_scores_std = []

for k,c in zip(classifiers.keys(), classifiers.values()):
    print(k)
    if k == "MLP":
        n_jobs = 1
    else:
        n_jobs= 4  
    cv_scores = cross_val_score(c, sub_X_train, sub_y_train, n_jobs=n_jobs, cv=5)
    print(cv_scores)
    cv_scores_mean.append(np.mean(cv_scores))
    cv_scores_std.append(np.std(cv_scores))
#print(cv_scores_mean)
#print(cv_scores_std)



from sklearn.model_selection import GridSearchCV
fine_tuning_best_model = []
fine_tuning_best_score = []
selected = ["GradBoost", "XGBCBoost", "RandomForestClassifier"]
gboost_param_grid = {'loss' : ["exponential", "log_loss"],
              'n_estimators' : [100,200,300],
              'learning_rate': [0.1, 0.05, 0.01],
              'max_depth': [4, 8],
              'min_samples_leaf': [100,150],
              'max_features': [0.3, 0.1] 
              }
xgBoost_param_grid = {'booster':('gbtree', 'gblinear', 'dart')}
rf_param_grid = {'n_estimators' : [20, 50, 100, 200, 500, 1000],
                 'max_depth' : [5, 10, 30, 50, 100],
                 'criterion': ['gini', 'entropy']}
param_grids = [gboost_param_grid, xgBoost_param_grid, rf_param_grid]


def ParamGridSearch(classifier, param_grid, X, y, n_jobs):
    gs = GridSearchCV(classifier, param_grid = param_grid, cv=5, scoring="accuracy", verbose = 1, n_jobs=n_jobs)
    gs.fit(X, y)
    return gs.best_estimator_, gs.best_score_


for i,c_name in enumerate(selected):
    print(c_name)
    c = classifiers[c_name]
    n_jobs = 8
    best_model, best_score = ParamGridSearch(c, param_grids[i], sub_X_train, sub_y_train, n_jobs)
    fine_tuning_best_model.append(best_model)
    fine_tuning_best_score.append(best_score)
    print(best_score)



best_results = pd.DataFrame({"CrossValMeans":fine_tuning_best_score, "Classifier":selected})
best_results = best_results.sort_values('CrossValMeans')
sns.barplot(x="CrossValMeans", y="Classifier", data = best_results)

from sklearn.model_selection import learning_curve

def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                        n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 5)):
    """Generate a simple plot of the test and training learning curve"""
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")
    plt.legend(loc="best")
    return plt

for i,c in enumerate(fine_tuning_best_model):
    f = plot_learning_curve(c, selected[i], sub_X_train, sub_y_train, n_jobs=1)
    f.show()
'''


# check version number
import imblearn
from imblearn.over_sampling import SMOTE
from collections import Counter

y_train_smote = X_train['clinical_status']
X_train_smote = X_train.sample(frac=1)
X_train_smote.drop('clinical_status', axis=1, inplace=True)

# transform the dataset
# oversample = SMOTE()
# X_train_smote, y_train_smote = oversample.fit_resample(X_train_smote, y_train_smote)

# summarize the new class distribution
counter = Counter(y_train_smote)
print(counter)

dataset_smote = X_train_smote.copy()
dataset_smote['clinical_status']=y_train_smote
dataset_smote = dataset_smote.sample(frac=1)
dataset_smote.head()
print(dataset_smote.corr()['clinical_status'][((dataset_smote.corr()['clinical_status'] >= 0.5 ) | (dataset_smote.corr()['clinical_status'] < -0.6))])
most_correlated_smote = dataset_smote.corr()['clinical_status'][((dataset_smote.corr()['clinical_status'] >= 0.5 ) | (dataset_smote.corr()['clinical_status'] < -0.6))].index


# we plot the distribution for all variables to check outliers
for c in dataset_smote.columns[:-1]:
    sns.catplot(x="clinical_status", y=c, data=dataset_smote, kind='box')
    plt.show()



def IQR_outlier_detection(data, feature, className, classVal, threshold):
    feature_data = data[feature].loc[data[className] == classVal].values
    # get the first and third quartile
    q25, q75 = np.percentile(feature_data, 25), np.percentile(feature_data, 75)
    # get interquartile range
    iqr = q75 - q25
    # multiply range by threshold
    limit = iqr * threshold
    # get upper and lower limit
    lower = q25 - limit
    upper = q75 + limit
    outliers = [x for x in feature_data if x < lower or x > upper]
    print("Number of outliers detected: " + str(len(outliers)))
    data = data.drop(data[(data[feature] > upper) | (data[feature] < lower)].index)
    return data


# Detect and remove outliers
for c in dataset_smote.columns[:-1]:
    dataset_smote = IQR_outlier_detection(dataset_smote, c, 'clinical_status', 1, 1.5)

print(dataset_smote['clinical_status'].value_counts())


# Model configuration
def KerasMLP(hidden=1000, activation='relu', learning_rate=0.01):
    model = Sequential([
        layers.Dense(hidden, activation=activation),
        layers.Dense(1, activation='sigmoid')
    ])
    opt = optimizers.Adam(learning_rate=learning_rate)
    model.compile(
        optimizer=opt,
        loss='binary_crossentropy',
        metrics=["accuracy"]
    )
    return model

classifiers = {
    "KNearest": KNeighborsClassifier(),
    "Support Vector Classifier": SVC(),
    "GradBoost": GradientBoostingClassifier(),
    "XGBCBoost" : XGBClassifier(),
    "RandomForestClassifier": RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=6),
}

# prepare data
y_train = X_train['clinical_status']
X_train.drop('clinical_status', axis=1, inplace=True)

from sklearn.model_selection import cross_val_score

cv_scores_mean = []
cv_scores_std = []

for k,c in zip(classifiers.keys(), classifiers.values()):
    print(k)
    if k == "MLP":
        n_jobs = 1
    else:
        n_jobs= 4  
    cv_scores = cross_val_score(c, X_train, y_train, n_jobs=n_jobs, cv=10)
    print(cv_scores)
    cv_scores_mean.append(np.mean(cv_scores))
    cv_scores_std.append(np.std(cv_scores))


from sklearn.model_selection import GridSearchCV

#Defining optimizable parameters in a grid for each model
fine_tuning_best_model = []
fine_tuning_best_score = []
selected = ["GradBoost", "XGBCBoost", "RandomForestClassifier"]
gboost_param_grid = {'loss' : ["exponential", "log_loss"],
              'n_estimators' : [100,200,300],
              'learning_rate': [0.1, 0.05, 0.01],
              'max_depth': [4, 8],
              'min_samples_leaf': [100,150],
              'max_features': [0.3, 0.1] 
              }
xgBoost_param_grid = {'booster':('gbtree', 'gblinear', 'dart')}
rf_param_grid = {'n_estimators' : [10, 100, 200],
                 'max_depth' : [3, 5, 10, None],
                 'max_features' : [3,4,5,6,7],
                 'min_samples_leaf' : [1,2,3],
                 'min_samples_split' : [2,3]}
param_grids = [gboost_param_grid, xgBoost_param_grid, rf_param_grid]


def ParamGridSearch(classifier, param_grid, X, y, n_jobs):
    gs = GridSearchCV(classifier, param_grid = param_grid, cv=5, scoring="accuracy", verbose = 1, n_jobs=n_jobs)
    gs.fit(X, y)
    return gs.best_estimator_, gs.best_score_

# Perform the parameter optimization
for i,c_name in enumerate(selected):
    print(c_name)
    c = classifiers[c_name]
    n_jobs = 8
    best_model, best_score = ParamGridSearch(c, param_grids[i], X_train, y_train, n_jobs)
    fine_tuning_best_model.append(best_model)
    fine_tuning_best_score.append(best_score)
    print(best_score)


# Plot the best result for each model
best_results = pd.DataFrame({"CrossValMeans":fine_tuning_best_score, "Classifier":selected})
best_results = best_results.sort_values('CrossValMeans')
sns.barplot(x="CrossValMeans", y="Classifier", data = best_results)
plt.show()

from sklearn.model_selection import learning_curve

def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                        n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 5)):
    """Generate a simple plot of the test and training learning curve"""
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")
    plt.legend(loc="best")
    return plt

for i,c in enumerate(fine_tuning_best_model):
    f = plot_learning_curve(c, selected[i], X_train, y_train, n_jobs=1)
    f.show()


from sklearn import metrics
import numpy as np

def calc_metrics(labels_test, test_probs, threshold = 0.5):
    scores = [1 if x>=threshold else 0 for x in test_probs]
    auc = metrics.roc_auc_score(labels_test, test_probs)
    kappa = metrics.cohen_kappa_score(labels_test,scores)
    confusion = metrics.confusion_matrix(labels_test,scores, labels=list(set(labels_test)))
    print('thresh: %.2f, kappa: %.3f, AUC test-set: %.3f'%(threshold, kappa, auc))
    print(confusion)
    print(metrics.classification_report(labels_test,scores))
    return

for i, c in enumerate(fine_tuning_best_model[0:]):
    print(c)
    # Get prediction probabilities for the test set
    test_probs = c.predict_proba(X_test)[:,1]
    # Print confusion matrix and classification metrics
    calc_metrics(y_test, test_probs, threshold = 0.5)




from sklearn.metrics import RocCurveDisplay
for i, c in enumerate(fine_tuning_best_model[0:]):
    #preds = c.predict(X_train)
    #fpr, tpr, thresold = roc_curve(y_train, log_reg_pred)
    RocCurveDisplay.from_estimator(c, X_test, y_test)
    plt.show()






# Feature importance visualization
feature_importances = pd.DataFrame({
    "Feature": X_train.columns,
    "Importance": c.feature_importances_
}).sort_values(by="Importance", ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(data=feature_importances.head(20), x="Importance", y="Feature", palette="viridis")
plt.title("Top 20 Feature Importances")
plt.tight_layout()
plt.show()


'''
# Predict probabilities on the test set
y_pred_proba = SMOTE_SRF.predict_proba(over_X_test)[:, 1]

# Calculate AUC-ROC
roc_auc = roc_auc_score(over_y_test, y_pred_proba)
print(f"AUC-ROC: {roc_auc}")

# Plot the ROC curve
fpr, tpr, thresholds = roc_curve(over_y_test, y_pred_proba)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color="blue", label=f"AUC = {roc_auc:.2f}")
plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()
'''




#Use SMOTE to oversample the minority class in the train set
oversample = SMOTE(random_state=42)
over_X_train, over_y_train = oversample.fit_resample(X_train, y_train)

#Build SMOTE SRF model
SMOTE_SRF = RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=6, max_features=7)
#Create Stratified K-fold cross validation
cv = RepeatedStratifiedKFold(n_splits=25, n_repeats=12, random_state=42)
scoring = ('f1', 'recall', 'precision')
#Evaluate SMOTE SRF model
scores = cross_validate(SMOTE_SRF, X_test, y_test, scoring=scoring, cv=cv)
#Get average evaluation metrics
print('Mean f1: %.3f' % mean(scores['test_f1']))
print('Mean recall: %.3f' % mean(scores['test_recall']))
print('Mean precision: %.3f' % mean(scores['test_precision']))

#Randomly spilt dataset to test and train set
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
#Train SMOTE SRF
SMOTE_SRF.fit(over_X_train, over_y_train)
#SMOTE SRF prediction result
y_pred = SMOTE_SRF.predict(over_X_test)


# Check new class distribution after SMOTE
print("Class distribution after SMOTE:\n", pd.Series(over_y_test).value_counts())


# Feature importance visualization
feature_importances = pd.DataFrame({
    "Feature": over_X_train.columns,
    "Importance": SMOTE_SRF.feature_importances_
}).sort_values(by="Importance", ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(data=feature_importances.head(20), x="Importance", y="Feature", palette="viridis")
plt.title("Top 20 Feature Importances")
plt.tight_layout()
plt.show()

# Predict probabilities on the test set
y_pred_proba = SMOTE_SRF.predict_proba(over_X_test)[:, 1]

# Calculate AUC-ROC
roc_auc = roc_auc_score(over_y_test, y_pred_proba)
print(f"AUC-ROC: {roc_auc}")

# Plot the ROC curve
fpr, tpr, thresholds = roc_curve(over_y_test, y_pred_proba)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color="blue", label=f"AUC = {roc_auc:.2f}")
plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()

# Classification report and confusion matrix
print("\nClassification Report:\n", classification_report(over_y_test, y_pred))


conf_matrix = confusion_matrix(over_y_test, y_pred)
plt.figure(figsize=(6, 6))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", xticklabels=["AWoD", "AWD"], yticklabels=["AWoD", "AWD"])
plt.title("Confusion Matrix")
plt.ylabel("True Label")
plt.xlabel("Predicted Label")
plt.tight_layout()
plt.show()

