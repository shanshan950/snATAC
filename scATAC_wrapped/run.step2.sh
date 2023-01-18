#!/bin/bash
picard=/mnt/rstor/genetics/JinLab/cxw486/mysoftware/picard-tools-1.93
scATAClib=/mnt/rstor/genetics/JinLab/ssz20/zshanshan/differentiation_ATAC/processed/Columbia_scATAC_1029_lane1/scATAC_wrapped
genome=$1

if [ $genome == hg19 ];then
	bowtie2Index=/mnt/rstor/genetics/JinLab/ssz20/zshanshan/new_ChIPSeq/bowtie2_hg19_index/hg19
	chromSize=/mnt/rstor/genetics/JinLab/xxl244/Reference_Indexes/hg19/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chrom.sizes
fi
if [ $genome == mm10 ];then
	bowtie2Index=/mnt/rstor/genetics/JinLab/cxw486/Genome/bowtie2index/mm10/mm10
	chromSize=/mnt/rstor/genetics/JinLab/xxl244/Reference_Indexes/mm10_bowtie_index/mm10.chrom.sizes
fi

currpath=`pwd`
## Run basic pipline
for folder in `ls | grep Folder.`;do
	echo $folder
	cd ./$folder
	name=${folder/Folder./}
	$scATAClib/scATAC.basic.sh $scATAClib $genome $bowtie2Index ${name} 1000
	$scATAClib/Bam2bw.sh $picard $scATAClib $genome $chromSize ${name}
	$scATAClib/CalculatPEAKTSS.sh ${name}.${genome}"_"peaks.narrowPeak $genome $scATAClib
	$scATAClib/ReadsQC.sh $name $genome
	cd $currpath
	echo $folder `cat $folder/RawQC | awk '{if($1!="")print}'| tail -1 | awk '{print $1}'` `cat $folder/Demit.SPLIT.${genome}/bedsummary | awk '{if($1>=1000)print}' | wc -l` >> summaray.table
done
