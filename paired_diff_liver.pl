#!/usr/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Calculate the distribution of the major allele differences between the
# paired liver samples
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $dir=shift||die "perl paired_diff_liver.pl [dir include all the ssp] [outdir] [distribution filename]\n";
my $outdir=shift||die "perl paired_diff_liver.pl [dir include all the ssp] [outdir] [distribution filename]\n";
my $distfile=shift||die "perl paired_diff_liver.pl [dir include all the ssp] [outdir] [distribution filename]\n";

my @samples;
my %hash;
my %hash1;

opendir DH,$dir||$!;
while(my $file=readdir(DH)){
	if ($file=~/^Li/) {
        $file=~/(\S+)\_(\S+)\.ssp$/;
        my $sample=$1;
        my $rep=$2;
		push @samples, $sample unless (grep {$sample eq $_} @samples);
		open IN,"$dir/$file"||$!;
		while(<IN>){
		    next if ($_=~/^pos/);
			my @array=split /\t/,$_;
			$hash{"$sample,$rep,$array[0]"}=$array[3] if ($array[2]);
			$hash1{"$sample,$rep,$array[0]"}=$_ if ($array[2]);
		}
		close IN;
	}
}
closedir DH;

open LOG,">$distfile";

foreach my $sample (@samples){
    open OUT,">$outdir/$sample.paired.diff"||$!;
	my $diff=0;
	my @array;
	for(my $pos=1;$pos<=16569;$pos++){
		next if (($pos>=16181 and $pos<=16195) || ($pos>=303 and $pos<=315));
                next if ($pos>512 && $pos<525);
                next if($pos>3104 && $pos<3110);
                next if($pos>8269 && $pos<8291);
                next if($pos>16566 or $pos<4);
                
                
		if(defined($hash{"$sample,1,$pos"}) and defined($hash{"$sample,2,$pos"})){
			if ($hash{"$sample,1,$pos"} ne $hash{"$sample,2,$pos"}){
				$diff++;
				push @array,$pos;
                print OUT $sample,"\t1\t", $hash1{"$sample,1,$pos"};
                print OUT $sample,"\t2\t", $hash1{"$sample,2,$pos"};
			}
			my %tmp1;
			my %tmp2;
			my @arr1=split /\t/,$hash1{"$sample,1,$pos"};
			my @arr2=split /\t/,$hash1{"$sample,2,$pos"};
			my @alleles;
			$tmp1{$arr1[3]}=$arr1[6]/$arr1[2];
			$tmp1{$arr1[4]}=$arr1[7]/$arr1[2];
			$tmp2{$arr2[3]}=$arr2[6]/$arr2[2];
			$tmp2{$arr2[4]}=$arr2[7]/$arr2[2];
 			my $frediff;
			if($arr1[3] eq $arr2[3]){
			      $frediff=abs($tmp1{$arr1[3]}-$tmp2{$arr2[3]});
			}
			elsif($arr1[3] eq $arr2[4]){
			      $frediff=abs($tmp1{$arr1[3]}-$tmp2{$arr2[4]});
			}
			elsif($arr1[4] eq $arr2[3]){
			      $frediff=abs($tmp1{$arr1[4]}-$tmp2{$arr2[3]});
			}
			elsif($arr1[4] eq $arr2[4]){
			      $frediff=abs($tmp1{$arr1[4]}-$tmp2{$arr2[4]});
			}
			else{
			      die qq(!!! $hash1{"$sample,1,$pos"}\n $hash1{"$sample,2,$pos"}\n);
			}
            print LOG $frediff,"\n";
		}
		
	}
    close OUT;
	print $sample,"\t",$diff,"\t";
	map {print $_,"\t"} @array if (@array);
	print "\n";
}
close LOG;
