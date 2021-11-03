#!/bin/Rscript

suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))

FUNG_DIR <- 'mappings/'

fungi <- c('Chame1','Truan1','Phapo1','Macpha1','Parch1','Zalva1')
#i<-'Chame1'
for (i in fungi) {
	outputDir<-paste(FUNG_DIR,i,'/',sep='')
	file<-read.csv(file=paste('mappings/',i,'/matrix.csv',sep=''), header=TRUE, sep=",", row.names=1)
	file<-file[grep("gene*",rownames(file)),]
	sampleTable <- data.frame(condition = factor(c("c","c","c","p","p","p")))
	dds<-DESeqDataSetFromMatrix(countData = file,colData = sampleTable,design = ~condition)
	dds$condition <- relevel(dds$condition, ref = "c")
	analysis <- DESeq(dds)
	res <- results(analysis)
	results.all <- lfcShrink(analysis, coef="condition_p_vs_c", type="apeglm")
	write.table(results.all, paste(outputDir,i,'.all.csv',sep=''))
	results.de <- subset(results.all, padj < 0.05)
	write.table(results.de, paste(outputDir,i,'.de.csv',sep=''))
	results.oe <- subset(results.all, padj < 0.05 & log2FoldChange > 0)
	write.table(results.oe, paste(outputDir,i,'.oe.csv',sep=''))
	results.ue <- subset(results.all, padj < 0.05 & log2FoldChange < 0)
	write.table(results.ue, paste(outputDir,i,'.ue.csv',sep=''))
	pdf(paste(outputDir,i,'.volcano.pdf',sep=''))
	plotMA(results.all, ylim=c(-2,2), alpha=0.05)
	dev.off()
	pdf(paste(outputDir,i,'.pca.pdf',sep=''))
	rld <- rlog(dds)
	plotPCA(rld)
	dev.off()
	print(paste(i,' ok',sep=''))
}
