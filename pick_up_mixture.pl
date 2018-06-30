#!/usr/local/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Remove samples which pass a threshold regarding the number of
# heteroplasmies that are caused with a high likelihood by a
# contamination from a different haplogroup
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $dir=shift||die "perl pick_up_mixture.pl [directory] [number of hg mutation allowed]\n";
my $no=shift||"5";
$dir=~s/\/$//;


print "#Minor allele at >=$no positions define one haplogroup\n";
opendir DH,$dir||$!;
while(my $file=readdir(DH)){
    next if ($_=~/-mixture.log$/); 
	my $input=$dir."/".$file;
	open IN,$input||$!;
	$file=~/([A-Za-z0-9_-]+)(_2_20_20)?(_(pois|emp))?\.(log|hg)/;
	my $name=$1;
	my $fre;
	my $no_count;
	while(<IN>){
		if($_=~/\t hg observed >1:\w+->(\d+)/){
			$no_count=$1;
		}
		if($_=~/>minor_frequency\t(0.\d+)/){
			$fre=$1;
		}
	}
	if($no_count>=$no){
		print $name,"\t","$fre\t","$no_count\n";
	}
	close IN;
}
closedir DH;
