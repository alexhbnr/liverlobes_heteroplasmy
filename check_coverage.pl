#!/usr/local/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Identify heteroplasmies that have < 20% or > 200% of the average
# coverage of an individual
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $file=shift||die "perl check_coverage.pl [full list file] [coverage file]\n";
my $cov=shift||"coverage.csv";

my %coverage;

open IN,$cov||$!;
while(<IN>){
	chomp;
	my @array=split /\t/,$_;
	$coverage{$array[0]}=$array[1];
}
close IN;

my $count=0;
open IN, $file||$!;
while(<IN>){
		if ($_=~/^#/){
			print $_;
			next;
		}
		my @array=split /\t/,$_;
		$array[0]=~/^(\S+)_(pois|emp)\.log/;
		my $sample=$1;
		if ($array[9]/$coverage{$sample}>0.2 and $array[9]/$coverage{$sample}<2){
			print $_;
		}
		else{
			$count++;
		}

}
close IN;

print "# $count are removed due to the lower coverage < 20%\n";
