#!/bin/bash

## Download fungal genomic data from JGI (gff3 and fasta)

./jgi_dl.sh Chame1 Parch1 Phapo1 Truan1 Macpha1 Sorhu1 Zalva1
rm *cookies*; rm */*.xml ; rm */*.sh; rm */*.gff.gz; rm */*MitoAssembly*
mv */* ./
find . -type d -empty -delete
gunzip *.gz
ls *.fasta | while read line ; do mv $line $(python -c "print('$line'.split('_')[0])").fasta ; done

## Download Arabidopsis genomic data from TAIR (gff3 and fasta)

curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > Ath.gff3
curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > Ath.fasta
python -c "F=open('Ath.fasta').read().replace('>','>Chr').replace('>Chrmitochondria','>ChrM').replace('>Chrchloroplast','>ChrC'); print(F)" > Ath.modified.fasta 
mv ./Ath.modified.fasta ./Ath.fasta
