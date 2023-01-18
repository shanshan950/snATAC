#!/bin/bash
scATAClib=scATAC_wrapped/
fqgz1=$1
fqgz2=$2
fqgz3=$3
fqgz4=$4
xlsxFile=$5

## Unzip the fastq.gz files
for file in $fqgz1 $fqgz2 $fqgz3 $fqgz4;do
	cat $file | gunzip > ${file/.gz/}
done
##1  Generate BCindex
Rscript $scATAClib/makedemultiplex.R $xlsxFile
##2 Demultiplexing and Make single cell tag file
$scATAClib/GetTag.pl ${fqgz1/.gz/} ${fqgz4/.gz/} ${fqgz2/.gz/} ${fqgz3/.gz/} `ls *BCindex`
