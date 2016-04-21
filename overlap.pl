#!/usr/bin/perl 

use strict; 
use warnings;

my $hs = shift @ARGV; 
my $mm = shift @ARGV; 

open HS,"<",$hs or die "$!"; 
open MM,"<",$mm or die "$!"; 

my %hs=(); 
my %mm=(); 

while (<HS>)  { 
  chomp; 
  my @tt = split /\t+/;
  my @kk = split /,/,$tt[1]; 
  $hs{$tt[0]}=\@kk; 
}

while (<MM>) { 
  chomp; 
  my @tt = split /\t+/;
  my @kk = split /,/,$tt[1];
  $mm{$tt[0]}=\@kk;
}

foreach my $k1 (keys %hs) {
  foreach my $k2 (keys %mm) {
    my @hs1 = @{$hs{$k1}};
    my @mm1 = @{$mm{$k2}};
    $_=uc for @mm1; 
    my (@union, @intersection, @difference);
    my %count = ();
    foreach my $element (@hs1, @mm1) { $count{$element}++ }
    foreach my $element (keys %count) {
      push @union, $element;
      push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
    }
    printf "%d\t%d\t%d\t%d\t%d\n",$k1,$k2,scalar @hs1,scalar @mm1,scalar @intersection; 
  }
} 

close HS; 
close MM;  
