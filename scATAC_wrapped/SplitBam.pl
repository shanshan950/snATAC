#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: ./SplitBam.pl <cutoff> <header><STD> \n".
"This script is to split the bam into sub-bam, each of which is for a cell\n".
"Output a large number of sam\n";
my $cutoff=shift @ARGV;
my $headerfile=shift @ARGV;
my $preCell=<>;
chomp $preCell;
my @Preline=split(/\s/,$preCell);
my $PrecellfieldN=@Preline;
my $PreIdentity="$Preline[$PrecellfieldN-1]";
my $output;
open $output,">",$PreIdentity or die();
print $output "$preCell\n";
my $CycleN=0;
my $readssummary;
open $readssummary,">","readssummary";
while(<>){
  $CycleN=$CycleN+1;
  chomp $_;
  my @curline=split(/\s/,$_);
  my $fieldsnum=scalar @curline;
  my $CellIdentity="$curline[$fieldsnum-1]";
  if ($CellIdentity eq $PreIdentity){
    print $output "$_\n";
  }else{
    close $output;
    if ($CycleN<$cutoff){
      unlink $PreIdentity;
    }else{
      my $dedupename="$PreIdentity.dedup.bam";
      my $command1="cat " .$headerfile . " " . $PreIdentity . " | samtools sort  >tmp.bam";
      my $command2="java -jar -Xmx4g /mnt/rstor/genetics/JinLab/cxw486/mysoftware/picard-tools-1.93/MarkDuplicates.jar I=tmp.bam O=" . $dedupename ." M=dump REMOVE_DUPLICATES=true";
      system($command1);
      system($command2);
      # my $uniqreadsN=system("samtools view -c " . $dedupename);
      print $readssummary "$PreIdentity\t$CycleN\n";
      unlink $PreIdentity;
    }
    open $output,">",$CellIdentity or die();
    print $output "$_\n";
    $PreIdentity=$CellIdentity;
    $CycleN=0;
  }
}

if ($CycleN<$cutoff){
  unlink $PreIdentity;
}else{
  my $dedupename="$PreIdentity.dedup.bam";
  my $command1="cat " .$headerfile . " " . $PreIdentity . " | samtools sort  >tmp.bam";
  my $command2="java -jar -Xmx4g /mnt/rstor/genetics/JinLab/cxw486/mysoftware/picard-tools-1.93/MarkDuplicates.jar I=tmp.bam O=" . $dedupename ." M=dump REMOVE_DUPLICATES=true";
  system($command1);
  system($command2);
  # my $uniqreadsN=system("samtools view -c " . $dedupename);
  print $readssummary "$PreIdentity\t$CycleN\n";
  unlink $PreIdentity;
}
close $output;
close $readssummary;
unlink "tmp.bam";
unlink "dump";
