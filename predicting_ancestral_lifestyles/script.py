import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import cross_val_score, KFold
from sklearn.feature_selection import RFECV

df=pd.read_csv('Orthogroups.GeneCount.csv', sep='\t').set_index('Unnamed: 0').drop(columns='Total').T
pos={int(c[2:])+1:c for c in df.columns}
gl=pd.read_csv('gainLossMP.2.00099.AncestralReconstructSankoff.txt', comment='#',sep='\t')
gl['protein']=gl.POS.map(pos)

gl=gl[['Node','protein','State']]
gl=gl[~gl['Node'].isin(df.index)]
gl=gl.pivot(index='Node', columns='protein', values='State')

df[df>1]=1
meta=pd.read_csv('120 genomes Endophyte NewFungalLifestyles_v2_Fantin16Oct19.csv')
df['Lifestyle']=df.index.map(meta.set_index('JGI ID')['New lifestyle'].to_dict())

#df['Lifestyle_factor']=
#print(pd.factorize(df['Lifestyle']).to_dict()[0])
#factorMapping=df.set_index('Lifestyle_factor')['Lifestyle'].to_dict()
#print(df[['Lifestyle','Lifestyle_factor']])
#df=df.drop(columns=['Lifestyle']).rename(columns={'Lifestyle_factor':'Lifestyle'})


features = df.columns[:-1]
clf = RandomForestClassifier(n_jobs=40, random_state=0)
#rfe = RFECV(clf, step=10, cv=KFold(n_splits=df.shape[0]),min_features_to_select=10,n_jobs=40)
clf.fit(df.drop(columns=['Lifestyle']), df['Lifestyle'])
scores = cross_val_score(clf,df.drop(columns=['Lifestyle']) ,df['Lifestyle'], cv=KFold(n_splits=df.shape[0]))
print('score= ', sum(scores)/len(scores))
factorMapping={pd.factorize(df['Lifestyle'])[0][i]: df['Lifestyle'][i] for i in range(len(df))}
#print(clf.classes_)



estimations=pd.DataFrame(clf.predict_proba(gl[features])).rename(index=str, columns=factorMapping)
estimations.index=gl.index
estimations.columns=clf.classes_
estimations.to_csv('gloome_predicted_lifestyles_og.csv')

