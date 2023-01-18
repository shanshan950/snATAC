#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: ./GetTag.pl <Read1.fq><Read2.fq><i7 index fastq > <i7 index fastq > <debarcode file1><debarcode file2><debarcode file3>...\n".
          "Based on debarcoding file, to write fq1, fq2 and barcode file into different directories\n".
          "The barcode file format is name\tT7P7T5P5".
          "out put is 3 files into each directory\n";
my $R1fqname=shift@ARGV;
my $R2fqname=shift@ARGV;
my $i7fqname=shift@ARGV;
my $i5fqname=shift@ARGV;
my @Barcodenames=@ARGV;
open R1,$R1fqname or die($usage);
open R2,$R2fqname or die($usage);
open i7,$i7fqname or die($usage);
open i5,$i5fqname or die($usage);


#################################################  define a subrotine hd for hamming distance
sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
################################################### Make reference for barcode correction
#### Input the barcode reference, the order of input should be T7(Tagmentation tagi7) P7(PCRtagi7) T5(Tagmentation tagi5) P5(PCRtagi5)
my @P5ref = ('TAGATCGC','CTCTCTAT','TATCCTCT','AGAGTAGA','GTAAGGAG','ACTGCATA','AAGGAGTA','CTAAGCCT','CGTCTAAT','TCTCTCCG','TCGACTAG','TTCTAGCT','CCTAGAGT','GCGTAAGA','CTATTAAG','AAGGCTAT','GAGCCTTA','TTATGCGA','TTTACATC','AAACAAAC','AGCAAATT','CGTTGGTG','GACTTCTT','GGTCACTC','CAATCTGC','CCGAGGAA','CGAACGGT','CGGTATCA','CTTGCGAC','GAGTGTAG','GCTGAACC','GTATACTG','GTGGTGCA','GTTTCGTT','TACGCCAC','TCAAACAT','TCCCATGT','TCGTGCTC','AGGGTCTG','ATGCTTAA');
my @P7ref = ('TCGCCTTA','CTAGTACG','TTCTGCCT','GCTCAGGA','AGGAGTCC','CATGCCTA','GTAGAGAG','CCTCTCTG','AGCGTAGC','CAGCCTCG','TGCCTCTT','TCCTCTAC','TCATGAGC','CCTGAGAT','TAGCGAGT','GTAGCTCC','TACTACGC','AGGCTCCG','GCAGCGTA','CTGCGCAT','GAGCGCTA','CGCTCAGT','GTCTTAGG','ACTGATCG','TAGCTGCA','CTCATTGT','GAGTGTAG','GCTAGGTT','GGAAGGAT','GGTCACTC','GTTTCGTT','TCAAACAT','TCCCATGT','TGTTTCCC','AAACAAAC','TTTACATC');
my @T5ref = ('AAGATC','AATTCC','CTAAGC','CTTGTG','GAATAG','GATGCT','TCCAAG','TGGGTT','AAACGA','TGTTGG','AGAGCT','ATTCGG','GAGCAA','GGCTTT','GTGGTA','TCAACT');
my @T7ref = ('AAACAC','AAAGGG','AACCGA','CCAATA','CCTCTT','CGTGAT','GCTGAA','GGACAA','GGCATT','TGTCTA','TTACCC','TTTCAG');

###################################################MAke hash for demultiplexing
my %DecodeHASH;
my %FHR1;
my %FHR2;
my %FHbarcode;
# my @FHbarcode;
foreach my $barcodegroup (@Barcodenames){
  system("mkdir Folder.$barcodegroup");
  open $FHR1{$barcodegroup},">","Folder.$barcodegroup/$barcodegroup.R1.fq" or die();
  open $FHR2{$barcodegroup},">","Folder.$barcodegroup/$barcodegroup.R2.fq" or die();
  open $FHbarcode{$barcodegroup},">","Folder.$barcodegroup/$barcodegroup.barcode" or die();
  # print "$barcodegroup\n";
  my $curbarfile;
  open $curbarfile,$barcodegroup or die($usage);
  my @P7;
  my @P5;
  my $n=1;
  while(<$curbarfile>){
    chomp $_;
    my @curline=split(/\s/,$_);
    if($n<=12){
    push @P7, $curline[1];
  }else{
    push @P5, $curline[1];
    }
    $n=$n+1;
  }
  foreach my $N7 (@P7){
     foreach my $N5 (@P5){
       # print "$N7$N5\n";
       $DecodeHASH{"$N7$N5"}=$barcodegroup;
     }
  }

}

# for (keys %DecodeHASH){
#   print "$_\t$DecodeHASH{$_}\n";
#
# }


#######################################################################Start to read the fastq
my $failN=0;
my $total=0;
while(!eof(i7) and !eof(i5) and !eof(R1) and !eof(R2)){
  my $R1line1=<R1>;
  chomp $R1line1;
  my $R1line2=<R1>;
  chomp $R1line2;
  my $R1line3=<R1>;
  chomp $R1line3;
  my $R1line4=<R1>;
  chomp $R1line4;

  my $R2line1=<R2>;
  chomp $R2line1;
  my $R2line2=<R2>;
  chomp $R2line2;
  my $R2line3=<R2>;
  chomp $R2line3;
  my $R2line4=<R2>;
  chomp $R2line4;


  my $f1line1=<i7>;
  chomp $f1line1;
  my @name=split(/\s/,$f1line1);
  my $f1line2=<i7>;
  chomp $f1line2;
  my $f1line3=<i7>;
  chomp $f1line3;
  my $f1line4=<i7>;
  chomp $f1line4;

  my $f2line1=<i5>;
  chomp $f2line1;
  my $f2line2=<i5>;
  chomp $f2line2;
  my $f2line3=<i5>;
  chomp $f2line3;
  my $f2line4=<i5>;
  chomp $f2line4;
my $T7=substr $f1line2,0,6;
my $P7=substr $f1line2,21,8;
my $T5=substr $f2line2,0,6;
my $P5=substr $f2line2,20,8;

$T7=reverse $T7;
$T7=~ tr/ATGCatgc/TACGtacg/;
$P7=reverse $P7;
$P7=~ tr/ATGCatgc/TACGtacg/;
$T5=reverse $T5;
$T5=~ tr/ATGCatgc/TACGtacg/;
$P5=reverse $P5;
$P5=~ tr/ATGCatgc/TACGtacg/;

#########################################################Correct T7
my @distances;
my $closest=6;
my $closestBC;
my $secondclosest=6;
foreach  (@T7ref)
{
  my $dis=hd($T7,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  $T7=$closestBC;
  # print "T7:Z:$closestBC\t";
}else{
  $T7="NNNNNN";
  # print "T7:Z:ND\t";
}

#########################################################Correct P7
$closest=6;
$secondclosest=6;
foreach  (@P7ref)
{
  my $dis=hd($P7,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  $P7=$closestBC;
  # print "P7:Z:$closestBC\t";
}else{
  $P7="NNNNNNNN";
  # print "P7:Z:ND\t";
}

#########################################################Correct T5
$closest=6;
$secondclosest=6;
foreach  (@T5ref)
{
  my $dis=hd($T5,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  $T5=$closestBC;
  # print "T5:Z:$closestBC\t";
}else{
  $T5="NNNNNN";
  # print "T5:Z:ND\t";
}

#########################################################Correct P5
$closest=6;
$secondclosest=6;
foreach  (@P5ref)
{
  my $dis=hd($P5,$_);
        if($dis<$closest){
    $closest=$dis;
    $closestBC=$_;
  }elsif($dis>$closest && $dis<$secondclosest){
    $secondclosest=$dis;
  }
}
if($closest==0 || ($closest<=2 && ($secondclosest-$closest)>=2)){
  $P5=$closestBC;
  # print "P5:Z:$closestBC\n";
}else{
  $P5="NNNNNNNN"
  # print "P5:Z:ND\n";
}
if (defined $DecodeHASH{"$P7$P5"}){
my $group=$DecodeHASH{"$P7$P5"};

print {$FHR1{"$group"}} "$R1line1\n$R1line2\n$R1line3\n$R1line4\n";
print {$FHR2{"$group"}} "$R2line1\n$R2line2\n$R2line3\n$R2line4\n";
print {$FHbarcode{"$group"}} "$name[0]\t$T7$P7$T5$P5\n";
# close $curR1;
# close $curR2;
# print "$P7$P5\t$group\n";
}else{
  #print "$P7$P5\n";
  $failN++;
}
$total++;
}
my $summary;
open $summary, ">", "Demultiplexsummary" or die();
print $summary "In total $total reads, $failN failed to be demultiplexed\n";
