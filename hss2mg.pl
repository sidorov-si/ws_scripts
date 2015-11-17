#!/usr/bin/perl 

## special fur conversion of Hummin symbols to murine entrez IDs 

use strict; 
use warnings; 
use List::MoreUtils qw(uniq);

my $sym = shift @ARGV; 
my $ann = shift @ARGV; ## ID, refseq

open SYM,"<",$sym or die "$!"; 
open ANN,"<",$ann or die "$!"; 

my %g2s;

while (<ANN>) {
  chomp;
  my @tt=split /\t+/;
  $g2s{uc($tt[15])}=$tt[1];
}
 
printf STDERR "Found %d unique symbol ids.\n",scalar keys %g2s;  

while (<SYM>) {
  chomp; 
  print STDOUT "$g2s{$_}\n" if (defined $g2s{$_}); 
}

close ANN;  
close SYM;
