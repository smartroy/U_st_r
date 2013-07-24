#!/usr/bin/perl -w
use diagnostics;
use strict;

my $input=shift @ARGV;
open IN_DATA,"<",$input;
my $count=0;
my @temp;
my $sum;

while(<IN_DATA>){
	chomp;

    @temp=split /\s+/, $_;
    $sum+=$temp[0];
    $count++;
}
print $input."   ".$sum/$count."\n";
