#!/bin/bash
name=$1
genome=$2
Rawreads=`samtools view -c  $name.$genome.PE.bam`
UniqmapReads=`samtools view -c  $name.$genome.PE.sctag.bam`
DmitReads=`samtools view -c  $name.$genome.PE.sctag.nomit.bam`
echo -e "Rawreads\tUniqmap$genome\tDemit$genome\n"  >RawQC
echo -e "$Rawreads\t$UniqmapReads\t$DmitReads\n" >>RawQC
