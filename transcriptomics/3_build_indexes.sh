#!/bin/bash

MAIN_DIR=$(python -c "print('/'.join('${BASH_SOURCE[0]}'.split('/')[:-2]))")
FILES_DIR=$MAIN_DIR/files
INDEX_DIR=$MAIN_DIR/indexes

pairs=(
    "Ath_Chame1"
    "Ath_Macpha1"
    "Ath_Parch1"
    "Ath_Phapo1"
    "Ath_Sorhu1"
    "Ath_Truan1"
    "Ath_Zalva1"
)


#GENERATE EXONS AND SPLICING SITES FILES
ls $FILES_DIR/*.gtf | while read line ; do hisat2_extract_splice_sites.py $line > $line.ss ; hisat2_extract_exons.py $line > $line.ex ; done

mkdir $INDEX_DIR

#GENERATE INDEX FOR EACH PAIR OF GENOMES
for i in "${pairs[@]}"; do
    p1=$(python -c "print('$i'.split('_')[0])")
    p2=$(python -c "print('$i'.split('_')[1])")
    cat $FILES_DIR/$p1.gtf.ss $FILES_DIR/$p2.gtf.ss > $FILES_DIR/$i.ss
    cat $FILES_DIR/$p1.gtf.ex $FILES_DIR/$p2.gtf.ex > $FILES_DIR/$i.ex
    hisat2-build -p 60 --ss $FILES_DIR/$i.ss --exon $FILES_DIR/$i.ex $FILES_DIR/$p1.fasta,$FILES_DIR/$p2.fasta $INDEX_DIR/$i
done



