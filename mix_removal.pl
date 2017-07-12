#!/usr/local/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Script to remove all heteroplasmies from the final list that haven't been
# flagged as quality check failed
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $input=shift||die "perl mix_heteroplasmy_remove.pl [heteroplasmy list] [sample to remove]\n";
my $file=shift||"mixture/test/mix.list";

my @list;
open IN,$file||die "can not find $file\n";
while(<IN>){
	next if ($_=~/^#/);
	next if ($_=~/^\s/);
	$_=~/^(\w+)_/;
	push @list,$1;
}
close IN;

my $count;
open IN, $input||$!;
while(my $l=<IN>){
	if ($l=~/^#/){
		print $l;
		next;
	}
	if(&check($l,[@list])){
		$count++;
		next;
	}
	else{
		print $l;
	}
}

print "# $count heteroplasmies are removed, because they are found in ","@list\t","\n";

sub check(){
	my $query=shift;
	my $array=shift;
	foreach my $int (@$array){
		return 1 if ($query=~/^$int\_/);
	}
	return 0;
}
