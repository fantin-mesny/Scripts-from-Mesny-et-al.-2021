from sklearn import svm
from sklearn import metrics
from sklearn.feature_selection import RFE, RFECV, f_classif, SelectFdr
from sklearn.pipeline import make_pipeline
import numpy as np
import pandas as pd
import sys
import subprocess
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score, KFold

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
    
################################################### PARAMETERS ###################################################
Dir='/biodata/dep_psl/grp_hacquard/Fantin/SVM_enrichedOGs/SVM-RFE/git/'
ogCount=Dir+'Orthogroups.GeneCount.csv'
ls=Dir+'data.csv'
#################################################### PARSING #####################################################
df0=pd.read_csv(ogCount, sep='\t').set_index('Unnamed: 0').drop(columns='Total').T.astype(float)
meta=pd.read_csv(ls)

df0['Lifestyle']=df0.index.map(meta.set_index('id')['Lifestyle'].to_dict())
df0['Lifestyle']=df0['Lifestyle'].str.replace('Arabidopsis mycobiota member', 'Endophyte')

for i in df0.index:
    if df0.loc[i,'Lifestyle']!='Endophyte':
        df0.loc[i,'Lifestyle']='Other'

scaler = StandardScaler()
df = pd.DataFrame(scaler.fit_transform(df0.drop(columns='Lifestyle')))
df.index=df0.index
df.columns=df0.columns[:-1]
df=df.merge(df0[['Lifestyle']], left_index=True, right_index=True)

df['Lifestyle']=df['Lifestyle'].map({'Endophyte':1,'Other':0})
        
###################################################  MODEL ####################################################
print('\n ########## MODEL ##########')
print('Built on '+str(len(df.columns)-1))


anova_filter = SelectFdr(f_classif, alpha=0.05)
clf = svm.SVC(kernel='linear')
rfe = RFECV(clf, step=10, cv=KFold(n_splits=df.shape[0]),min_features_to_select=10,n_jobs=40)
anova_svm_rfe = make_pipeline(anova_filter, rfe)
anova_svm_rfe.fit(df.drop(columns=['Lifestyle']),df['Lifestyle'])
print('number of orthogroups kept in the model: ',anova_svm_rfe.steps[1][1].n_features_)
print('R2 score= ', sum(anova_svm_rfe.steps[1][1].grid_scores_)/len(anova_svm_rfe.steps[1][1].grid_scores_))

################################################## ANALYSING OGs ##################################################
print('\n ########## ANALYSING THE OGs ##########')
anova=pd.DataFrame()
anova['pvalue']=anova_svm_rfe.steps[0][1].pvalues_
anova['scores']=anova_svm_rfe.steps[0][1].scores_
anova['ogs']=df.columns[:-1]

new_features=[]
for bool, feature in zip(anova_svm_rfe.steps[0][1].get_support(), df.columns[:-1]):
    if bool:
        new_features.append(feature)
anova=anova[anova['ogs'].isin(new_features)]

cols=anova_svm_rfe.steps[1][1].get_support(indices=True)
cols = [new_features[i] for i in cols]
coeffs=pd.DataFrame()
coeffs['svm_coefficient']=list(anova_svm_rfe.steps[1][1].estimator_.coef_[0])
coeffs['ogs']=cols
coeffs=coeffs.merge(anova, on='ogs', how='left').set_index('ogs')[['scores','pvalue','svm_coefficient']]
means=df.groupby('Lifestyle').mean()
for i in coeffs.index:
    coeffs.loc[i,'enrichedInEndophytes']=means.loc[1,i]>means.loc[0,i]
    
coeffs=coeffs[['enrichedInEndophytes','scores','pvalue','svm_coefficient']]
coeffs.to_csv(Dir+'retainedOGs.csv')
