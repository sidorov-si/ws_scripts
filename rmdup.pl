#!/usr/bin/perl

use strict;
use warnings;
## don't use it on very large BAMs, it will explode in your face

my $bam = shift @ARGV; 
my $cut = shift @ARGV; 
my %count; ## hash of string counts 

open BAM,"samtools view $bam |" or die "$!";
die "ERROR: please specify the cutoff."  if ( !defined $cut); 

system "samtools view -H $bam";  
while(<BAM>){
  my $line=$_; 
  my @tt = split(/\t+/); 
  $count{$tt[9]}++; 
  print $line if ($count{$tt[9]}<=$cut); 
}

close BAM; 
