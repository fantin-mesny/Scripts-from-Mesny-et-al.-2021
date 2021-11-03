import os
import pandas as pd

MAIN_DIR='/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
FUNG_DIR=MAIN_DIR+'/mappingsFungi'
PLANT_DIR=MAIN_DIR+'/mappings'

for fung in os.listdir(FUNG_DIR):
	dfs=[]
	print(fung)
	for counts in os.listdir(FUNG_DIR+'/'+fung):
		if '.counts.tsv' in counts and 'summary' not in counts and counts.split('.')[0]!='':
			letter=counts.split('.')[0]
			cond='c'
			df=pd.read_csv(FUNG_DIR+'/'+fung+'/'+counts,sep='\t',comment='#')
			df=df.set_index('Geneid')
			df=df[sorted([c for c in df.columns if '/' in c])]
			dfs.append(df.rename(columns={col:cond+letter for col in df.columns}))
	for counts in os.listdir(PLANT_DIR):
		if '.counts.tsv' in counts and 'summary' not in counts and 'Ath_'+fung in counts and fung+'_' in counts:
			letter=counts.split('_')[1]
			cond='p'
			df=pd.read_csv(PLANT_DIR+'/'+counts,sep='\t',comment='#')
			df=df.set_index('Geneid')
			df=df[sorted([c for c in df.columns if '/' in c])]
			dfs.append(df.rename(columns={col:cond+letter for col in df.columns}))



	Df=pd.concat(dfs, axis=1)
	#Df=Df.rename(columns={c:c.split('/')[-1] for c in Df.columns})
	Df=Df[sorted([col for col in Df.columns])]
	Df.to_csv(FUNG_DIR+'/'+fung+'/matrix.csv')

