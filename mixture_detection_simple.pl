#!/usr/local/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Detect heteroplasmies that are caused with a high likelihood by a
# contamination from a different haplogroup
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;
use lib "/mnt/genotyping/sk_pipelines/source/heteroplasmy/mingkun_libs";
use phylotree;
use Getopt::Std;

use vars qw( $opt_m $opt_i);

&init();

my $ref_file="/mnt/genotyping/sk_pipelines/source/heteroplasmy/humanmt.fasta";
my $ref;
open IN,$ref_file||$!;
while(<IN>){
	next if($_=~/^>/);
	chomp;
    $_ =~ s/\s+$//g;	
	$ref.=$_;
}
my @ref_array=split //,$ref;

my $phylotree_o=phylotree->new($opt_m);
my $minor_fre;
my $count;
my $total;
open IN,$opt_i||die $!;
my $indi="tested";

my %hash;
while(<IN>){
	chomp;
    next if ($_=~/^#/);
	$total++;
	my @array_c=split /\t/,$_;
	
	my $second_c=$array_c[12];
	my $pos=$array_c[0];
	my $ref_allele=$ref_array[$pos-1];
	
	$second_c=$array_c[10] if ($second_c eq $ref_allele);
	
	my $c2=$ref_allele.$second_c;

	print $pos,"\t";
	
	if (defined $phylotree_o->{$pos}->{$c2} and  @{$phylotree_o->{$pos}->{$c2}}<=100){
		$minor_fre+=$array_c[9];
		$count++;
		print $c2,": ","@{$phylotree_o->{$pos}->{$c2}}","; ";
		push @{$hash{'minor'}},@{$phylotree_o->{$pos}->{$c2}};
	}
	elsif(defined $phylotree_o->{$pos}->{$c2} and  @{$phylotree_o->{$pos}->{$c2}}>100){
		$minor_fre+=$array_c[9];
		$count++;
		print $c2,": ",">100hg","; ";
		push @{$hash{'minor'}},@{$phylotree_o->{$pos}->{$c2}};
	}	
	else{
		print $c2,":none",";";
	}
	print "\n";
}

close IN;
print ">minor_frequency\t",$minor_fre/$count,"\t","$count/$total\n" if ($count>0 and $total>0);

foreach my $indi (keys %hash){
	print ">$indi","\n";
	my $classify=&class($hash{$indi});
	my @array_1;
	my @array_2;
	foreach my $int (sort {$classify->{$b}<=>$classify->{$a}} keys %$classify){
		if($classify->{$int}>1){
			push @array_2,$int."->".$classify->{$int};
		}
		else{
			push @array_1,$int;
		}
	}
	
	print "\t hg observed >1:","@array_2","\n";
	print "\t hg observed 1:","@array_1","\n";
}


close IND;

sub class(){
	my $array=shift;
	my %hash;
# 	my %result;
	foreach (@$array){
		$hash{$_}++;
	}
	return \%hash;
}

sub _min(){
	my $no1=shift;
	my $no2=shift;
	if($no1>$no2){
		return $no2;
	}
	else{
		return $no1;
	}
}



sub arrange(){
	my $array=shift;
	my %hash;
	my @result;
	for(my $i=0;$i<@$array;$i++){
		next unless ($array->[$i]);
		$hash{$array->[$i]}=$array->[$i+1];
		$i++;
	}
	foreach (sort {$hash{$b}<=>$hash{$a}} keys %hash){
		push @result,$_;	
	}
	return [@result];
}

sub init(){
	getopts( 'm:i:' );
	unless(defined ($opt_i)){
		die "perl mixture_detection_simple.pl -i -m
			-i input file[heteroplasmy file]
			-m phylotree file[made by Mark B]
			input file is ssp file, line:45 is the requirement for hetero calling\n"
	}
	unless(defined ($opt_m)){
		$opt_m="/mnt/genotyping/sk_pipelines/source/heteroplasmy/mtDNAtree.build13.REF.txt";
	}
}
