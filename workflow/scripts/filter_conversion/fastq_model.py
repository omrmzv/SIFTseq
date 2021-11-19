

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from sklearn.model_selection import cross_validate
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import plot_roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing 
from sklearn import linear_model



df = pd.read_csv('FASTQ_dataset.csv')
df.astype({'seqtype': 'int64'}).dtypes

print(df.shape)
print(df.head())

X = df[['C_count']]  #,'T_count', 'length','avg_prob_error']]
print(X.head())

#0=bisulfite treated, 1=standard which corresponds in this case to 0=uncontaminted, 1=contaminated
y = df['seqtype'] 

# Optional: select k best features to include
#X_new = SelectKBest(chi2, k=2).fit_transform(X, y)
#X_train, X_test, y_train, y_test = train_test_split(X_new, y, random_state = 0, test_size=0.2)

# Train-test split 80-20 and scale data:
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 0, test_size=0.2)

# Print size of X to show which features:
print(X_train.shape)

# Implement and train a few different models:
svc = make_pipeline(preprocessing.StandardScaler(), SVC(kernel='linear')).fit(X_train, y_train)
rf = RandomForestClassifier(random_state=0, ccp_alpha=0.05).fit(X_train, y_train)
sgd = make_pipeline(preprocessing.StandardScaler(), linear_model.SGDClassifier()).fit(X_train, y_train)

# Print score:
print("SVM score, train: ",svc.score(X_train,y_train))
print("SVM score, test: ",svc.score(X_test,y_test))



# Do cross-validation on the whole dataset multiple times:
cross_val = np.zeros(3)
for i in range(3):
	cross_val[i] = np.mean(cross_validate(svc, X, y, cv=5)['test_score'])

print("avg cross validation score: ",np.mean(cross_val))
'''

# Plot of feature importances for the random forest
feat_importances = pd.Series(rf.feature_importances_, index=X.columns)
feat_importances.nlargest(20).plot(kind='barh')
plt.title("Random Forest Feature Importances")
plt.savefig("feature_importance_rf.pdf") 


# Plot of feature importances for the SVM
feat_importances = pd.Series(abs(svc.steps[1][1].coef_[0]), index=X.columns)
feat_importances.nlargest(20).plot(kind='barh')
plt.title("SVM Feature Importances")
plt.savefig("feature_importance_svm.pdf")


# plot values, colored by prediction
plt.scatter(X['C_count'],X['avg_prob_error'], c = svc.predict(X), alpha=0.5)
plt.xlabel("Total C count")
plt.ylabel("Average Probability of Prediction Error")
plt.title("SVM Prediction by C, Sequencing Error")
plt.savefig('pred_plot.pdf') 


# Plot and show ROC curve:
#svc_disp = plot_roc_curve(sgd, X_test, y_test)
#plt.savefig('sgd_ROC.pdf') 
'''


