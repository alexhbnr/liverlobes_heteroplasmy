#!/usr/local/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Calculate the average coverage from the SSP files
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $dir=shift || ".";
$dir=~s/\/$//;

my %coverage;
my %calling;
print "sample\tcoverage\tfrac_larger_500x\n";
opendir DH,$dir||$!;
while(my $file=readdir(DH)){
	 if($file=~/ssp$/){
		my $cov;
		my $call;
		$file=~/^(\S+)\.ssp/;
		my $sample=$1;
		open IN,$dir."/".$file||$!;
		while(<IN>){
			my @array=split /\t/,$_;
			$cov+=$array[2];
			$call++ if ($array[2]>=500);
		}
		close IN;
		$coverage{$sample}=$cov/16569;
		$calling{$sample}=$call/16569;
		$coverage{$sample}=sprintf("%.0f",$coverage{$sample});
		print $sample, "\t", $coverage{$sample}, "\t", $calling{$sample}, "\n";
	}
}
closedir DH;
