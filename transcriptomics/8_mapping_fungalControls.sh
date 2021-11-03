#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>mapping.log 2>&1

CORENUM=60

MAIN_DIR=$(python -c "print('/'.join('${BASH_SOURCE[0]}'.split('/')[:-2]))")
INDEX_DIR=$MAIN_DIR/indexes
FILES_DIR=$MAIN_DIR/files
MAPPING_DIR=$MAIN_DIR/mappingsFungi

mkdir $MAIN_DIR/mappingsFungi

ls $MAIN_DIR/reads | while read fungus ; do
	echo $fungus
	echo '  '
	index=Ath_$fungus
	FUNG_DIR=$MAIN_DIR/reads/$fungus
	python -c "import os; print(''.join(list(set([f.split('_')[1] for f in os.listdir('$FUNG_DIR') if 'trimmed.fastq' in f]))))" | sed -e 's/\(.\)/\1\n/g' | while read letter ; do 
		mkdir $MAPPING_DIR/$fungus
		cat $FUNG_DIR/*_$letter\_*.trimmed.fastq > $FUNG_DIR/$letter.all.trimmed.fastq;
		hisat2 -p $CORENUM -x $INDEX_DIR/$index -U $FUNG_DIR/$letter.all.trimmed.fastq -S $MAPPING_DIR/$fungus/$fungus\_$letter\_$index.sam;
		samtools sort -@ $CORENUM -o $MAPPING_DIR/$fungus/$fungus\_$letter\_$index.bam $MAPPING_DIR/$fungus/$fungus\_$letter\_$index.sam;
		samtools index $MAPPING_DIR/$fungus/$fungus\_$letter\_$index.bam;
		featureCounts -T $CORENUM -a $FILES_DIR/$index.gtf -o $MAPPING_DIR/$fungus/$letter.counts.tsv -M -O $MAPPING_DIR/$fungus/$fungus\_$letter\_$index.bam;
	done
done

