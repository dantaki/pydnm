#!/usr/bin/env python3
import pandas as pd
from sklearn.externals import joblib
from pyDNM.Backend import get_path
import os,sys
def classify(ofh=None,keep_fp=None):
	snv_clf = get_path()+'/pydnm.snv.clf.joblib'
	if not os.path.isfile(snv_clf):
		sys.stderr.write('FATAL ERROR: {} CLASSIFIER NOT FOUND\n'.format(snv_clf))
		sys.exit(1)
		
	clf = joblib.load(snv_clf)
	df = pd.read_csv(ofh,sep="\t")

	snv = df.loc[(df['ref'].str.len()==1) & (df['alt'].str.len()==1)]
	#indel = df.drop(snv.index)

	snv = snv.dropna(axis=0,subset=snv.columns[12:36])
	X_test = snv[snv.columns[12:36]].values

	snv['pred'] = clf.predict(X_test)
	snv['prob'] = clf.predict_proba(X_test)[:, 1]

	df = snv # TO DO: ADD IN THE INDELS!!!
	#df = pd.concat([snv,indel])
	if keep_fp == False:
		df = df.loc[df['pred']==1]
	df.to_csv(ofh,sep="\t",index=False)