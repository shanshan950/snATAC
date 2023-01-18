#!/usr/bin/perl
use strict;
use warnings;

my $usage="Usage: ./AddTag.pl <tag hash file> <bam by std>\n".
          "This script add customized tag on to bam\n".
          "Output is sam file\n";


########################## Make tag hash #########################################################
my $tagfilename=shift@ARGV;
my $tagfile;
open $tagfile,$tagfilename or die($usage);
my %taghash;
while(<$tagfile>){
chomp $_;
my @curline=split(/\s/,$_);
$taghash{$curline[0]}="$curline[1]";
}

# for (keys %taghash) {
#   print "$_\t$taghash{$_}\n";
# }
#print "tagHash completed;\n";
#########################Read bam file and add tag end the end ######################################
while (<>) {
chomp $_;
my @curline=split(/\s/,$_);
my $readname="@".$curline[0];
my $curtag=$taghash{$readname};
print "$_\tTA:Z:$curtag\n";
}
