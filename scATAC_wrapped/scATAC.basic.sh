#!/bin/bash
scATAClib=$1
genome=$2
Bowtie2Index=$3
name=$4
cutoff=$5 #1000

bowtie2 -X 1200  --very-sensitive -p 5 -x $Bowtie2Index -1 $name.R1.fq -2 $name.R2.fq |  samtools view -u -  |  samtools sort -  > $name.$genome.PE.bam

samtools view -H $name.$genome.PE.bam  > $name.$genome.PE.header

samtools  view -F 1804 -q30  $name.$genome.PE.bam |$scATAClib/AddTag.pl $name.barcode -| grep -v NNNNNN | rev | sort -k1,1 | rev | cat $name.$genome.PE.header - | samtools view -S -o $name.$genome.PE.sctag.bam

#Remove mitochondria
samtools view $name.$genome.PE.sctag.bam | grep -v "chrM" | cat $name.$genome.PE.header - | samtools view -S -o $name.$genome.PE.sctag.nomit.bam

# Split bam by cell tag
mkdir Demit.SPLIT.$genome
cd ./Demit.SPLIT.$genome
samtools view ../$name.$genome.PE.sctag.nomit.bam | $scATAClib/SplitBam.pl $cutoff ../$name.$genome.PE.header ## To split bam into many small ones, eahcis a cell

# remove PCR duplicates
rm -rf tmp.uniqReads
touch tmp.uniqReads
for dedup in `ls *.dedup.bam`;do
  UniqReads=`samtools view -c $dedup`
  echo -e "${dedup/.dedup.bam/}\t$UniqReads" >>tmp.uniqReads
  samtools sort -n $dedup |samtools view -bf 2  | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,".",1,$9}' >${dedup/.dedup.bam/paired.monoclonal.bed}
done

sort -k1,1 tmp.uniqReads > tmp.uniqReads.sort
sort -k1,1 readssummary >readssummary.sort
cut -f2 tmp.uniqReads.sort | paste readssummary.sort - >Summary
rm -rf readssummary
rm -rf tmp.uniqReads.sort
rm -rf tmp.uniqReads
rm -rf readssummary.sort

touch bedsummary
for bedfile in `ls *.bed`;do
	echo `wc -l $bedfile` >>bedsummary
done

