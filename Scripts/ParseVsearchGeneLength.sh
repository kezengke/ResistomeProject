#!/bin/bash

input_fasta=$1
output_csv=$2

echo "Header,Length" > $output_csv

while read line; do
    if [[ $line =~ ^\> ]]; then
        # Get header without the '>' sign
        header=$(echo $line | sed 's/^>//')
    else
        # Get length of the sequence
        length=${#line}
        # Output to csv
        echo "\"$header\",$length" >> $output_csv
    fi
done < $input_fasta
