#!/bin/bash

Peakfile=$1
genome=$2  
scATAClib=$3

##  to calculate for each single cell the number of fragment on TSS and on peaks
curdir=`pwd`
echo -e "Name\tFragments\tFragTSS\tFragPeak" > $curdir/Demit.SPLIT.${genome}/Summary.TSSPeak
for file in ./Demit.SPLIT.${genome}/*.bed;do
	monofragment=`wc -l $file| cut -d ' ' -f1`
	fragmentonTSS=`bedtools intersect -a $file -b $scATAClib/${genome}.refFlat.TSS.1000round -wa -u| wc -l`
	fragmentononPeak=`bedtools intersect -a $file -b $Peakfile -wa -u| wc -l`
	echo -e "${file/"./Demit.SPLIT.${genome}/"/}\t$monofragment\t$fragmentonTSS\t$fragmentononPeak" >> $curdir/Demit.SPLIT.${genome}/Summary.TSSPeak
done


