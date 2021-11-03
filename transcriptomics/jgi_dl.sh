#!/bin/bash


curl "https://signon.jgi.doe.gov/signon/create" --data-urlencode "login=**********@******.***" --data-urlencode "password=********" -c cookies > /dev/null


for i in $*; do # access each parameter
    if [ ! -d "$i" ]; then
        mkdir "$i"
    fi
    curl "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=$i" -b cookies > $i/$i.xml
    python3 /biodata/dep_psl/grp_hacquard/Fantin/scripts/jgi_dl_get_curl_files.py -xml $i/$i.xml -o $i
    chmod +x $i/download.sh
    $i/download.sh
    #gunzip -c $i/\*.gz > $i/proteins.fasta
    echo "$i data have been succesfully downloaded"
done


