#!/bin/Rscript

inputDir <-'readCounts/'
outputDir <- 'DESeqOutput/'
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))

for (i in list.files(inputDir)){
	file<-read.csv(file=paste(inputDir,i,sep=''), header=TRUE, sep=",", row.names=1)
	sampleTable <- data.frame(condition = factor(c("Mock","Mock","Mock","Fungus","Fungus","Fungus")))
	file<-file[grep("AT*",rownames(file)),]
	dds<-DESeqDataSetFromMatrix(countData = file,colData = sampleTable,design = ~condition)
	dds$condition <- relevel(dds$condition, ref = "Mock")
	analysis <- DESeq(dds)
	res <- results(analysis)
	results.all <- lfcShrink(analysis, coef="condition_Fungus_vs_Mock", type="apeglm")
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
}
