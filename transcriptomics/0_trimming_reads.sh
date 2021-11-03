cd raw/
ls *.fastq | while read line ; do trimmomatic SE $line $line.trimmed.fastq TRAILING:20 AVGQUAL:20 HEADCROP:10 MINLEN:100 ; done > trimming.log &
