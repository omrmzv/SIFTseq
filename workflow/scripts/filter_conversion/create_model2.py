#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import time
from sklearn.model_selection import cross_validate
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import plot_roc_curve, accuracy_score, roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn import linear_model
from sklearn.ensemble import GradientBoostingClassifier
from joblib import dump
from joblib import Parallel, delayed
import pickle


def get_reads(info):
	sample_id, taxid, label, n_reads = info

	filepath = './sample_output/filter_unconverted/grammy/' + sample_id + '.filtered'

	df = pd.read_csv(filepath, sep = '\t')

	df = df[["species", "C_count", "T_count"]]
	df = df.loc[df['species'] == int(taxid)]
	if n_reads != "all":
		n_reads = int(n_reads)
		if len(df.index) > n_reads:
			df = df.head(n_reads)
	df['label'] = int(label)
	df = df.drop(columns = ['species'])
	return(df)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i',"--input", type=str,help="input data file name, should be a tab-separated dataset", nargs = 1)
	parser.add_argument('-o','--output',type=str, help = 'name of pickle file to output model to', nargs = 1)
	parser.add_argument('-f', '--figure', type=str)
	parser.add_argument('-p', '--parallel', type=int, help = 'threads', default = 1, nargs=1)

	args = parser.parse_args()

	input_data = args.input[0]
	output = args.output[0]
	threads = args.parallel[0]
	fig = args.figure

	training_info = []
	with open(input_data) as f:
		next(f)
		for line in f:
			training_info.append(line.strip().split('\t'))

	reads = Parallel(n_jobs = threads)(delayed(get_reads)(i) for i in training_info)
	reads = pd.concat(reads)

	X = reads[['C_count', 'T_count']]
	y = reads[['label']].values.ravel()

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 1)

	svc = make_pipeline(preprocessing.StandardScaler(), SVC(kernel='linear')).fit(X_train, y_train)

	y_pred = svc.predict(X_test)
	print(accuracy_score(y_test, y_pred))

	fpr, tpr, _ = roc_curve(y_test, y_pred)
	roc_auc = auc(fpr, tpr)
	plt.figure()
	plt.plot(fpr, tpr, color = 'darkorange', label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.legend(loc="lower right")
	plt.savefig(fig, dpi=300)
	pickle.dump(svc, open(output, 'wb'))
