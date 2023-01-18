#! /bin/bash
picard=$1
scATAClib=$2
genome=$3
chromSize=$4
name=$5

##################################################### Make aggregate track
java -jar -Xmx4g $picard/MarkDuplicates.jar  I=$name.$genome.PE.bam O=$name.$genome.PE.dupMark.bam M=$name.$genome.PE.matrix REMOVE_DUPLICATES=true

samtools sort -n $name.$genome.PE.dupMark.bam | samtools view -bf 2  | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,".",1,$9}' |sort -k1,1 -k2,2n -k3,3n>$name.$genome.PE.monoclonal.bed

bedtools genomecov -bg -i $name.$genome.PE.monoclonal.bed -g $chromSize > $name.$genome.PE.monoclonal.bed.bg

bedGraphToBigWig $name.$genome.PE.monoclonal.bed.bg $chromSize $name.$genome.PE.monoclonal.bed.bw

macs2 callpeak -f BAMPE -t $name.$genome.PE.dupMark.bam -g hs  -n $name.$genome

less $name.$genome'_'peaks.narrowPeak | cut -f1-4 > $name.$genome.narrowpeak.bed

$scATAClib/bedToBigBed $name.$genome.narrowpeak.bed $chromSize $name.$genome.narroepeak.bb
