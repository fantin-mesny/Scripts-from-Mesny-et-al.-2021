#!/home/mesny/anaconda3/bin/python

import os
import sys
import pandas as pd

MAIN_DIR='/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
MAPPING_DIR=MAIN_DIR+'/mappings'
INDEX_DIR=MAIN_DIR+'/indexes'
OUTPUT_DIR=MAIN_DIR+'/readCounts'
try:
    os.mkdir(OUTPUT_DIR)
except:
    pass

conditions={c:[] for c in set([index.split('.')[0] for index in os.listdir(INDEX_DIR)])}
for i in os.listdir(MAPPING_DIR):
    if '.counts.tsv' in i and 'summary' not in i:
        conditions['_'.join(i.split('.')[0].split('_')[2:])].append(i)
        
for c in conditions:
    dfs=[]
    for f in conditions[c]:
        df=pd.read_csv(MAPPING_DIR+'/'+f,sep='\t',comment='#')
        df=df[[df.columns[0],df.columns[-1]]].rename(columns={df.columns[-1]:df.columns[-1].split('/')[-1].split('.')[0]}).set_index('Geneid')
        #print(df)
        dfs.append(df.T)
    Df=pd.concat(dfs).T
    Df[sorted([col for col in Df.columns if 'Mock' in col])+sorted([col for col in Df.columns if 'Mock' not in col])].to_csv(OUTPUT_DIR+'/'+c+'.csv')