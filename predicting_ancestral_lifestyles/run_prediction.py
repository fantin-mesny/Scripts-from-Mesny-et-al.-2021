import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, KFold

## Parse orthogroup table
df=pd.read_csv('Orthogroups.GeneCount.csv', sep='\t').set_index('Unnamed: 0').drop(columns='Total').T

# Parse Lifestyles
meta=pd.read_csv('lifestyles.csv')
df['Lifestyle']=df.index.map(meta.set_index('Fungus')['Lifestyle'].to_dict())
df=df.sort_values(by='Lifestyle')

## Parse Wagner-parsimony predictions
pred=pd.read_csv('wagner_parsimony_ancestral_reconstruction.csv',sep='\t',skiprows=2)
pred=pred.set_index('# Family')
nodeNames=pd.read_csv('nodeNameMapping.csv',dtype=str).set_index('COUNT')['ETE'].to_dict()
pred=pred.rename(columns=nodeNames)
pred=pred[pred.columns[124:]].T

#Remove MyM from training
mym=df[df['Lifestyle']=='Arabidopsis mycobiota member']
df=df[df['Lifestyle']!='Arabidopsis mycobiota member']
pred=pd.concat([pred,mym],sort=True).drop(columns=['Lifestyle'])

## Model training
features = df.columns[:-1]
print('Training Random Forests...')
clf = RandomForestClassifier(n_jobs=40,n_estimators=1000,random_state=1)
clf.fit(df.drop(columns=['Lifestyle']), df['Lifestyle'])
scores = cross_val_score(clf,df.drop(columns=['Lifestyle']) ,df['Lifestyle'], cv=KFold(n_splits=df.shape[0]))
print('score= ', sum(scores)/len(scores))

estimations=pd.DataFrame(clf.predict_proba(pred[features]))
estimations.index=pred.index
estimations.columns=clf.classes_
estimations.to_csv('predicted_lifestyles.csv')

