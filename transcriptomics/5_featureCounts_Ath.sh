#!/bin/bash

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>featureCounts.log 2>&1

CORENUM=60
MAIN_DIR=$(python -c "print('/'.join('${BASH_SOURCE[0]}'.split('/')[:-2]))")
FILES_DIR=$MAIN_DIR/files
INDEX_DIR=$MAIN_DIR/indexes
MAPPING_DIR=$MAIN_DIR/mappings
READS_DIR=$MAIN_DIR/reads

ls $MAPPING_DIR/*.bam | while read line ; do 
    INDEX=$(python -c "print('_'.join('$line'.split('/')[-1].split('.')[0].split('_')[-2:]))");
    featureCounts -T $CORENUM -a $FILES_DIR/$INDEX.gtf -o $line.counts.tsv -M -O $line ;
    echo $line'   DONE '
done
