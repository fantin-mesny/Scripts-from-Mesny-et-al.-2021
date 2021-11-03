import pandas as pd
import os
import seaborn as sns

Dir='mappings/'
Dir_idMapping='gff3/'

ogs=pd.read_csv('Orthogroups.csv',sep='\t').set_index('Unnamed: 0')
ogs=ogs[list(os.listdir(Dir))].fillna('')

meanLFCperOG=[]
for fung in os.listdir(Dir):
    print(fung)
    df=pd.read_csv(Dir+fung+'/'+fung+'.all.csv',sep=' ')
    df.loc[df['padj'] > 0.05, ['log2FoldChange']]=0
    df.loc[df['padj'].isna(), ['log2FoldChange']]=0
    idMapping=pd.read_csv(Dir_idMapping+[f for f in os.listdir(Dir_idMapping) if fung in f][0],sep='\t',comment='#',header=None)
    idMapping=idMapping[idMapping[2]=='gene']
    idMapping['protID']=idMapping[8].str.split('proteinId=').str[1].str.split(';').str[0]
    idMapping['geneID']=idMapping[8].str.split('ID=').str[1].str.split(';').str[0]
    geneToProt=idMapping.set_index('geneID')['protID'].to_dict()
    df['proteinId']=df.index.map(geneToProt).astype(int)
    df=df.set_index('proteinId')
    for og in ogs.index:
        for prot in ogs.loc[og,fung].split(', '):
            if '|' in prot:
                df.loc[int(prot.split('|')[2]),'OG']=og
    df=df[~df['OG'].isna()]
    meanLFCperOG.append(df[['OG','log2FoldChange']].rename(columns={'log2FoldChange':fung}).groupby('OG').mean())

out=pd.concat(meanLFCperOG,axis=1)
out['cumulated']=out.sum(axis=1)
out=out.sort_values(by='cumulated',ascending=False)
filteredDf=out[['Zalva1','Phapo1','Chame1','Truan1','Macpha1','Parch1']]
toKeep=(filteredDf>0) & (~filteredDf.isna())
filteredDf['nbSpDEing']=toKeep.sum(axis=1)
filteredDf['cumulatedLFC']=filteredDf.drop(columns=['nbSpDEing']).sum(axis=1)
filteredDf=filteredDf[filteredDf['nbSpDEing']==6]

filteredDf.to_csv('commonly_overexpressed_gene_families.csv)
