#!/usr/bin/perl -w
use strict;
use diagnostics;

my $input=shift @ARGV;
open IN_DATA,"<",$input;
open OUT_DATA,">",($input."_core");
my @temp;
my $i;
my $pwr;
while(<IN_DATA>){
	@temp=split /\s+/, $_;
	$pwr=0.0;
	for($i=0;$i<17;$i++){
		$pwr+=$temp[$i];
	}
	printf OUT_DATA ("%.6f\n",$pwr);
}