

import pandas as pd
import numpy as np
import argparse
import pickle



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-id',"--input_data", type=str,help="input data file path, should be a CSV or txt file", nargs = 1)
	parser.add_argument('-m',"--model", type=str,help="path to .joblib pickle object containing trained model", nargs = 1)
	parser.add_argument('-o','--output_data',type=str, help = 'path for dataset containing classifications of data to be saved to', nargs = 1)

	args = parser.parse_args()
	input_data = args.input_data[0]
	model = args.model[0]
	output = args.output_data[0]


	df = pd.read_csv(input_data, sep ='\t')

	X = df[['C_count','T_count']] #must be the same shape as was used to train the model

	# Load model:
	svc = pickle.load(open(model, 'rb'))

	# Make predictions and save to new dataset:
	out = svc.predict(X)
	df['SVC_classification'] = out #0 = true, 1 = contaminant

	species = df[["species", "SVC_classification"]]
	species = species.groupby(['species'])['SVC_classification'].mean().reset_index()
	species = species.rename(columns = {'species' : 'species', 'SVC_classification' : 'species_contamination_proba'})
	df = pd.merge(df, species, on = 'species')

	df.to_csv(output, sep = '\t', index=None)
