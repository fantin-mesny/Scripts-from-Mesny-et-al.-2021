library(vegan)
library(RVAideMemoire)

for name in c('cazymes','fcwde','lipases','pcwde','proteases','ssp','total_proteins_(OGs)') {
	cazy<-read.csv(paste(name,'.csv',sep=''), row.names='X')
	data<-read.csv('data.csv', row.names='X')
	cazyDist <- vegdist(cazy, method='jaccard')
	cazDistMatrix <- as.data.frame(as.matrix(cazyDist))
	perm <- adonis2(formula(paste('cazyDist~',names(data)[1],'+',names(data)[2],'+',names(data)[3],'+',names(data)[4],'+Lifestyle','+Lifestyle:',names(data)[1],'+Lifestyle:',names(data)[2],'+Lifestyle:',names(data)[3],'+Lifestyle:',names(data)[4],sep='')), data=data, permutations = 9999)
	capture.output(perm, file=paste(name,'_permanova.txt',sep=''))
	write.csv(cazDistMatrix,paste(name,'_distMatrix.csv',sep=''), row.names=TRUE)

	ls<-data['Lifestyle']
	pca<-dbrda(formula(paste('cazyDist~Condition(',names(data)[1],'+',names(data)[2],'+',names(data)[3],'+',names(data)[4],')+Lifestyle',sep='')),data=data)
	areatest<-ordiareatest(pca, ls$Lifestyle, area ="ellipse", kind = "sd", permutations = 9999)
	b<-betadisper(cazyDist,ls$Lifestyle,type="centroid")
	c<-vegdist(b$centroids, method="euclidian")
	capture.output(b, file=paste(name,'_betadisp.txt',sep=''))
	capture.output(areatest, file=paste(name,'_areatest.txt',sep=''))
	write.csv(as.matrix(c), paste(name,'_centroidDistances.txt',sep=''))

	permManova<-pairwise.perm.manova(cazyDist,data$Lifestyle,,nperm=99999)
	permManova <- as.data.frame(permManova[3])
	write.csv(permManova,paste(name,'_pairwise.csv',sep=''), row.names=TRUE)
}
