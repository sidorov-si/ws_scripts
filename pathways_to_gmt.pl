#!/usr/bin/perl 

use strict; 
use warnings; 

## simple script to convert biocyc pathway file (3-column version - path_id, path_name, gene_id)
## into a legit GMT file. 

my $pw = shift @ARGV; 

open PW,"<",$pw or die "$!"; 
my $id=""; 
my $name=""; 
my $genes=""; 
my %gmt=(); 

my @lines=<PW>; 

for (my $i=0; $i < scalar @lines; $i++) {
  chomp $lines[$i];
  my @tt = split /\t+/,$lines[$i]; 
  if ($i==0) {
    $id=$tt[0];
    $name=$tt[1]; 
    $genes=$tt[2]; 
  } elsif ($tt[0] ne $id) {
    print "$id\t$name\t$genes\n"; 
    $id=$tt[0];
    $name=$tt[1];
    $genes=$tt[2];
  } elsif ($i == scalar @lines -1) {
    $genes=$genes."\t".$tt[2];
    print "$id\t$name\t$genes\n";
  } elsif ($tt[0] eq $id) {
    $genes=$genes."\t".$tt[2]; 
  } 
} 

close PW; 
