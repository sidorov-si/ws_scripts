#!/usr/bin/perl 

use strict; 
use warnings; 

## this one is to convert AraCyc file into some resemblance of systematic gene naming. 

my $gmt = shift @ARGV; 
my $ann = shift @ARGV; 

open GMT,"<",$gmt or die "$!"; 
open ANN,"<",$ann or die "$!"; 

my %ann = (); 
my %gmt = (); 

while (<ANN>) { 
  chomp; 
  my @tt = split /\t+/; 
  $ann{$tt[1]} = $tt[0]; 
} 


while (<GMT>) { 
  my ($c1,$c2,$c3)=(0,0,0); 
  my @unclaimed; 
  chomp; 
  my @tt = split /\t+/; 
  my $pid = shift @tt; 
  my $name = shift @tt;
  my $newline = ""; 
  foreach my $gene (@tt) {
    if ($gene =~ /(AT\wG\d{5})/) {
      $newline=$newline."\t".$1; 
      $c1++;
      $gmt{$1} = "TURD";  
    } elsif (defined $ann{$gene}) {
      $newline=$newline."\t".$ann{$gene}; 
      $c2++;
      $gmt{$ann{$gene}} = "BIRD"; 
    } else {
      $c3++;
      push @unclaimed,$gene;  
    }
  }
  print STDERR "Found $c1 LOCUS IDs, $c2 defined SYMBOLS, $c3 undef SYMBOLS: @unclaimed\n"; 
  print "$pid\t$name$newline\n"; 
} 

print STDERR "==========================================================================================\n"; 
printf STDERR "The processed GMT file has %d unique gene IDs.\n",scalar keys %gmt; 

close GMT; 
close ANN;
