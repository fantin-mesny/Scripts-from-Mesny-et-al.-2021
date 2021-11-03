#!/bin/bash

ls *.gff3 | while read line ; do gffread -T $line -o $(python -c "print('$line'.split('.')[0].split('_')[0])").gtf ; done

rm *.gff3
