

import pandas as pd
import numpy as np
import argparse
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
from sklearn.ensemble import GradientBoostingClassifier
from joblib import dump


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-k',"--num_folds", help="number of folds for cross-validation", type=int, nargs=1)
	parser.add_argument('-n',"--reps", type=int,help="number of cross-validation repetitions to run", nargs = 1)
	parser.add_argument('-s',"--split", type=float, help="percentage of data used to train (between 0 and 1)", nargs=1)
	parser.add_argument('-i',"--input_data", type=str,help="input data file name, should be a tab-separated dataset", nargs = 1)
	parser.add_argument('-o','--output_name',type=str, help = 'name of pickle file to output model to', nargs = 1)

	args = parser.parse_args()
	k = args.num_folds[0]
	n = args.reps[0]
	s = args.split[0]
	input_data = args.input_data[0]
	output = args.output_name[0]


	df = pd.read_csv(input_data)
	print("Shape of dataset: ",df.shape)

	X = df[['C_count','T_count', 'length', 'avg_prob_error']]

	#0=bisulfite treated, 1=standard which corresponds in this case to 0=uncontaminted, 1=contaminated
	y = df['seqtype'] 

	# Train-test split:
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=s)


	# Implement and train a few different models:
	svc = make_pipeline(preprocessing.StandardScaler(), SVC(kernel='linear')).fit(X_train, y_train)
	rf = RandomForestClassifier(random_state=0, ccp_alpha=0.05).fit(X_train, y_train)
	sgd = make_pipeline(preprocessing.StandardScaler(), linear_model.SGDClassifier()).fit(X_train, y_train)
	gb = GradientBoostingClassifier(random_state=0).fit(X_train, y_train)


	# Do cross-validation on the whole dataset and print average score:
	cross_val = np.zeros(3)
	for i in range(n):
		cross_val[i] = np.mean(cross_validate(svc, X, y, cv=k)['test_score'])

	print("Average cross validation score: ",np.mean(cross_val))
	
	# Save model as pickle object
	# change the variable here to change the model type
	dump(svc, output+'.joblib') 
	print("Saved model to "+output+".joblib")
	
