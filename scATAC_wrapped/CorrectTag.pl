#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: ./CorrectTag.pl <STD> \n".
          "This script Correct the barcode\n".
          "Output is the corrected barcode for futher bam file annotation\n";

sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

#### Input the barcode reference, the order of input should be T7(Tagmentation tagi7) P7(PCRtagi7) T5(Tagmentation tagi5) P5(PCRtagi5)
open T7,"/mnt/NFS/homeGene/JinLab/cxw486/Chip-seq/ATAC/scATAC/scATAC19.4.26/test/T7.ref" or die();
open P7,"/mnt/NFS/homeGene/JinLab/cxw486/Chip-seq/ATAC/scATAC/scATAC19.4.26/test/P7.ref" or die();
open T5,"/mnt/NFS/homeGene/JinLab/cxw486/Chip-seq/ATAC/scATAC/scATAC19.4.26/test/T5.ref" or die();
open P5,"/mnt/NFS/homeGene/JinLab/cxw486/Chip-seq/ATAC/scATAC/scATAC19.4.26/test/P5.ref" or die();

my @T7ref;
my @P7ref;
my @T5ref;
my @P5ref;
while(<T7>){
  chomp $_;
  push(@T7ref,$_);
}
while(<P7>){
  chomp $_;
  push(@P7ref,$_);
}
while(<T5>){
  chomp $_;
  push(@T5ref,$_);
}

while(<P5>){
  chomp $_;
  push(@P5ref,$_);
}

while(<>){
chomp $_;
my ($Name,$T7BC,$P7BC,$T5BC,$P5BC)=split(/\s/,$_);
print "$Name\t";
my @distances;
my $closest=6;
my $closestBC;
my $secondclosest=6;
foreach  (@T7ref)
{
  my $dis=hd($T7BC,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  print "T7:Z:$closestBC\t";
}else{
  print "T7:Z:ND\t";
}


$closest=6;
$secondclosest=6;
foreach  (@P7ref)
{
  my $dis=hd($P7BC,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  print "P7:Z:$closestBC\t";
}else{
  print "P7:Z:ND\t";
}


$closest=6;
$secondclosest=6;
foreach  (@T5ref)
{
  my $dis=hd($T5BC,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  print "T5:Z:$closestBC\t";
}else{
  print "T5:Z:ND\t";
}


$closest=6;
$secondclosest=6;
foreach  (@P5ref)
{
  my $dis=hd($P5BC,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  print "P5:Z:$closestBC\n";
}else{
  print "P5:Z:ND\n";
}



# my $i7inner=substr $curline[0],0,6;
# # print "$i7inner\n";
# my @distances;
# my $closest=6;
# my $closestBC;
# my $secondclosest=6;

}





#
# my $file1;
# my $filename1=shift@ARGV;
# my @expecti7inner=@ARGV;
# open $file1, $filename1 or die $!;
# while(<$file1>){
# chomp $_;
# my @curline=split(/\s/,$_);
# my $i7inner=substr $curline[0],0,6;
# # print "$i7inner\n";
# my @distances;
# my $closest=6;
# my $closestBC;
# my $secondclosest=6;
# foreach  (@expecti7inner)
# {
#   my $dis=hd($i7inner,$_);
#         if($dis<$closest){
#     $closest=$dis;
#     $closestBC=$_;
#   }elsif($dis>$closest && $dis<$secondclosest){
#     $secondclosest=$dis;
#   }
# }
# if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
#   print "TA:Z:$closestBC\n";
# }else{
#   print "TA:Z:ND\n";
# }
# }
#
