from sklearn import svm
from sklearn.feature_selection import SelectFdr, f_classif, RFECV
from sklearn.pipeline import make_pipeline
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, KFold
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

################################################### PHENOTYPES ###################################################

phenotypes={ #lowP
    'Alternaria_alternata_MPI-PUGE-AT-0064':'Not Pathogen',
    'Alternaria_rosae_MPI-PUGE-AT-0040':'Not Pathogen',
    'Boeremia_exigua_MPI-SDFR-AT-0100':'Not Pathogen',
    'Chaetomium_funicola_MPI-SDFR-AT-0129':'Not Pathogen',
    'Chaetomium_globosum_MPI-SDFR-AT-0079':'Not Pathogen',
    'Chaetomium_megalocarpum_MPI-CAGE-AT-0009':'Not Pathogen',
    'Coprinopsis_phaeopunctatus_MPI-PUGE-AT-0042':'Not Pathogen',
    'Cylindrocarpon_olidum_MPI-CAGE-CH-0241':'Not Pathogen',
    'Dactylonectria_estremocensis_MPI-CAGE-AT-0021':'Pathogen',
    'Dactylonectria_macrodidyma_MPI-CAGE-AT-0147':'Pathogen',
    'Dendryphion_nanum_MPI-CAGE-CH-0243':'Not Pathogen',
    'Fusarium_commune_MPI-SDFR-AT-0072':'Not Pathogen',
    'Fusarium_equiseti_MPI-CAGE-AA-0113':'Not Pathogen',
    'Fusarium_oxysporum_MPI-CAGE-CH-0212':'Not Pathogen',
    'Fusarium_oxysporum_MPI-SDFR-AT-0094':'Pathogen',
    'Fusarium_redolens_MPI-CAGE-AT-0023':'Not Pathogen',
    'Fusarium_solani_MPI-SDFR-AT-0091':'Not Pathogen',
    'Fusarium_tricinctum_MPI-SDFR-AT-0044':'Pathogen',
    'Fusarium_tricinctum_MPI-SDFR-AT-0068':'Not Pathogen',
    'Fusarium_venenatum_MPI-CAGE-CH-0201':'Not Pathogen',
    'Ilyonectria_europaea_MPI-CAGE-AT-0026':'Pathogen',
    'Leotiomycetes_sp._MPI-SDFR-AT-0126':'Not Pathogen',
    'Leptodontidium_orchidicola_MPI-SDFR-AT-0119':'Not Pathogen',
    'Macrophomina_phaseolina_MPI-SDFR-AT-0080':'Not Pathogen',
    'Microdochium_trichocladiopsis_MPI-CAGE-CH-0230':'Not Pathogen',
    'Mortierella_elongata_MPI-CAGE-AA-0104':'Not Pathogen',
    'Neonectria_radicicola_MPI-CAGE-AT-0134':'Pathogen',
    'Oliveonia_pauxilla_MPI-PUGE-AT-0066':'Not Pathogen',
    'Paraphoma_chrysanthemicola_MPI-GEGE-AT-0034':'Not Pathogen',
    'Paraphoma_chrysanthemicola_MPI-SDFR-AT-0120':'Not Pathogen',
    'Phaeosphaeria_poagena_MPI-PUGE-AT-0046c':'Not Pathogen',
    'Plectosphaerella_cucumerina_MPI-CAGE-AT-0016':'Not Pathogen',
    'Plectosphaerella_cucumerina_MPI-SDFR-AT-0117':'Pathogen',
    'Pyrenochaeta_lycopersici_MPI-SDFR-AT-0127':'Not Pathogen',
    'Rhexocercosporidium_sp._MPI-PUGE-AT-0058':'Not Pathogen',
    'Sordaria_humana_MPI-SDFR-AT-0083':'Not Pathogen',
    'Stachybotrys_elegans_MPI-CAGE-CH-0235':'Not Pathogen',
    'Thanatephorus_cucumeris_MPI-SDFR-AT-0096':'Pathogen',
    'Truncatella_angustata_MPI-SDFR-AT-0073':'Not Pathogen',
    'Verticillium_dahliae_MPI-CAGE-AT-0001':'Not Pathogen',
    'Zalerion_varium_MPI-CAGE-AT-0135':'Not Pathogen'
}

phenotypesH={ #highP
    'Alternaria_alternata_MPI-PUGE-AT-0064':'Pathogen',
    'Alternaria_rosae_MPI-PUGE-AT-0040':'Pathogen',
    'Boeremia_exigua_MPI-SDFR-AT-0100':'Pathogen',
    'Chaetomium_funicola_MPI-SDFR-AT-0129':'Not Pathogen',
    'Chaetomium_globosum_MPI-SDFR-AT-0079':'Not Pathogen',
    'Chaetomium_megalocarpum_MPI-CAGE-AT-0009':'Not Pathogen',
    'Coprinopsis_phaeopunctatus_MPI-PUGE-AT-0042':'Not Pathogen',
    'Cylindrocarpon_olidum_MPI-CAGE-CH-0241':'Pathogen',
    'Dactylonectria_estremocensis_MPI-CAGE-AT-0021':'Pathogen',
    'Dactylonectria_macrodidyma_MPI-CAGE-AT-0147':'Pathogen',
    'Dendryphion_nanum_MPI-CAGE-CH-0243':'Not Pathogen',
    'Fusarium_commune_MPI-SDFR-AT-0072':'Pathogen',
    'Fusarium_equiseti_MPI-CAGE-AA-0113':'Not Pathogen',
    'Fusarium_oxysporum_MPI-CAGE-CH-0212':'Not Pathogen',
    'Fusarium_oxysporum_MPI-SDFR-AT-0094':'Pathogen',
    'Fusarium_redolens_MPI-CAGE-AT-0023':'Pathogen',
    'Fusarium_solani_MPI-SDFR-AT-0091':'Not Pathogen',
    'Fusarium_tricinctum_MPI-SDFR-AT-0044':'Pathogen',
    'Fusarium_tricinctum_MPI-SDFR-AT-0068':'Not Pathogen',
    'Fusarium_venenatum_MPI-CAGE-CH-0201':'Not Pathogen',
    'Ilyonectria_europaea_MPI-CAGE-AT-0026':'Pathogen',
    'Leotiomycetes_sp._MPI-SDFR-AT-0126':'Not Pathogen',
    'Leptodontidium_orchidicola_MPI-SDFR-AT-0119':'Not Pathogen',
    'Macrophomina_phaseolina_MPI-SDFR-AT-0080':'Not Pathogen',
    'Microdochium_trichocladiopsis_MPI-CAGE-CH-0230':'Not Pathogen',
    'Mortierella_elongata_MPI-CAGE-AA-0104':'Not Pathogen',
    'Neonectria_radicicola_MPI-CAGE-AT-0134':'Pathogen',
    'Oliveonia_pauxilla_MPI-PUGE-AT-0066':'Not Pathogen',
    'Paraphoma_chrysanthemicola_MPI-GEGE-AT-0034':'Pathogen',
    'Paraphoma_chrysanthemicola_MPI-SDFR-AT-0120':'Not Pathogen',
    'Phaeosphaeria_poagena_MPI-PUGE-AT-0046c':'Not Pathogen',
    'Plectosphaerella_cucumerina_MPI-CAGE-AT-0016':'Pathogen',
    'Plectosphaerella_cucumerina_MPI-SDFR-AT-0117':'Pathogen',
    'Pyrenochaeta_lycopersici_MPI-SDFR-AT-0127':'Not Pathogen',
    'Rhexocercosporidium_sp._MPI-PUGE-AT-0058':'Not Pathogen',
    'Sordaria_humana_MPI-SDFR-AT-0083':'Not Pathogen',
    'Stachybotrys_elegans_MPI-CAGE-CH-0235':'Not Pathogen',
    'Thanatephorus_cucumeris_MPI-SDFR-AT-0096':'Pathogen',
    'Truncatella_angustata_MPI-SDFR-AT-0073':'Not Pathogen',
    'Verticillium_dahliae_MPI-CAGE-AT-0001':'Not Pathogen',
    'Zalerion_varium_MPI-CAGE-AT-0135':'Not Pathogen'
}

################################################### PARAMETERS ###################################################
Dir='./'
jgi41=pd.read_csv('41fungiMetadata.csv').rename(index=str,columns={'folder_name':'Treatment'})
IDmapping=jgi41.set_index('Treatment')['jgi_id'].to_dict()
phenotypes_renamed={IDmapping[p]:phenotypes[p] for p in phenotypes}
phenotypesH_renamed={IDmapping[p]:phenotypesH[p] for p in phenotypesH}
ogCount='Orthogroups.GeneCount.csv'

#################################################### PARSING #####################################################
df0=pd.read_csv(ogCount, sep='\t').set_index('Unnamed: 0').drop(columns='Total').T
df0=df0[df0.index.isin(phenotypes_renamed)]
df0=df0.loc[(df0.sum(axis=1) != 0), (df0.sum(axis=0) != 0)]

scaler = StandardScaler()
df = pd.DataFrame(scaler.fit_transform(df0))
df.index=df0.index
df.columns=df0.columns
#df=df.merge(ses[['Hedges']], left_index=True, right_index=True)
df.to_csv(Dir+'std_df.csv')

############################################## MODEL LOW PHOSPHATE ###############################################
df['Phenotype']=df.index.map(phenotypes_renamed)
print('\n########## MODEL ON PHENOTYPE  LowP ##########')
print('\n ########## MODEL ##########')
print('Built on '+str(len(df.columns)-1))
print(df.shape[0])
anova_filter = SelectFdr(f_classif, alpha=0.05)
clf = svm.SVC(kernel='linear')
rfe = RFECV(clf, step=10, cv=KFold(n_splits=df.shape[0]),min_features_to_select=1,n_jobs=40)
anova_svm_rfe = make_pipeline(anova_filter, rfe)
anova_svm_rfe.fit(df.drop(columns=['Phenotype']),df['Phenotype'])
print('number of kept features',anova_svm_rfe.steps[1][1].n_features_)
print('score= ', sum(anova_svm_rfe.steps[1][1].grid_scores_)/len(anova_svm_rfe.steps[1][1].grid_scores_))

print('\n ## ANALYSING THE OGs ##')
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
means=df.groupby('Phenotype').mean()

for i in coeffs.index:
    coeffs.loc[i,'meanGreaterInPathogens']=means.loc['Pathogen',i]>means.loc['Not Pathogen',i]
    
coeffs=coeffs[['meanGreaterInPathogens','scores','pvalue','svm_coefficient']]

print(coeffs)
############################################## MODEL HIGH PHOSPHATE ###############################################

#df=df.drop(columns=['Phenotype'])
#print(df)
#df['Phenotype']=df.index.map(phenotypesH_renamed)
#print(df)
#print('\n########## MODEL ON PHENOTYPE  HighP ##########')
#print('\n ########## MODEL ##########')
#print('Built on '+str(len(df.columns)-1))
#anova_filter = SelectFdr(f_classif, alpha=0.05)
#clf = svm.SVC(kernel='linear')
#rfe = RFECV(clf, step=10, cv=KFold(n_splits=df.shape[0]),min_features_to_select=1,n_jobs=40)
#anova_svm_rfe = make_pipeline(anova_filter, rfe)
#anova_svm_rfe.fit(df.drop(columns=['Phenotype']),df['Phenotype'])
#print('number of kept features',anova_svm_rfe.steps[1][1].n_features_)
#print('score= ', sum(anova_svm_rfe.steps[1][1].grid_scores_)/len(anova_svm_rfe.steps[1][1].grid_scores_))

