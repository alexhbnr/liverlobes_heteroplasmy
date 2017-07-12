#!/usr/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Pairwise comparisons of all samples of a tissue regarding the overlap
# of major and minor alleles and the minor allele frequencies in order
# to identify cross-contamination
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $dir=shift||die "perl search_dataset_mixture.pl [dir has all ssp] [dir has all heteroplasmy list] [outdir]\n";
my $dir_hetero=shift||die "perl search_dataset_mixture.pl [dir has all ssp] [dir has all heteroplasmy list] [outdir]\n";
my $outdir=shift||die "perl search_dataset_mixture.pl [dir has all ssp] [dir has all heteroplasmy list] [outdir]\n";

my @samples;

my %major;
my %minor;
my %fre;
opendir DH,$dir ||$!;
while(my $file=readdir(DH)){
    if ($file=~/ssp$/){
	    my $sample;
	    my $tissue;
		if ($file=~/^Bl/) {
            $file=~/(\D+)\_(\d+)/;
            $sample=$2;
            $tissue=$1;
		} else {
            $file=~/(\D+)\_(\d+\_\d)/;
            $sample=$2;
            $tissue=$1;
		}
		push @samples, $sample unless (grep {$_ eq $sample} @samples);
		open IN, "$dir"."/"."$file"||$!;
		while(<IN>){
            next if ($_=~/^pos/);
			my @array=split /\t/,$_;
			next if (($array[0]>=16181 and $array[0]<=16195) || ($array[0]>=303 and $array[0]<=315));
			next if ($array[0]>=8269 and $array[0]<=8290);
            $major{"$sample,$tissue,$array[0]"}=$array[3];	
            $minor{"$sample,$tissue,$array[0]"}=$array[4];
            $fre{"$sample,$tissue,$array[0]"}=$array[7]/$array[2]||0;
		}
		close IN;
	}
}
closedir DH;

opendir DH,$dir_hetero ||$!;
while(my $file=readdir(DH)) {
    if($file=~/log$/) {
	    my $sample;
	    my $tissue;
		if ($file=~/^Bl/) {
            $file=~/(\D+)\_(\d+)/;
            $sample=$2;
            $tissue=$1;
		} else {
            $file=~/(\D+)\_(\d+\_\d)/;
            $sample=$2;
            $tissue=$1;
		}

        my @hetero;
        open IN,"$dir_hetero"."/"."$file"||$!;
        while(<IN>){
            my @array=split /\t/,$_;
            push @hetero,$array[1];
        }
        close IN;
	
        my $output="$outdir"."/"."$tissue"."_"."$sample".".mix";
        open OUT,">$output" ||$!;
        
        foreach my $int (@samples){
            next if ($sample eq $int);
            my $total_diff=0;
            my $found=0;
            my $significant=0;
            my $fre_sum=0;
            for(my $i=1;$i<=16569;$i++){
                next unless (defined $major{"$sample,$tissue,$i"});
                next unless (defined $major{"$int,$tissue,$i"});
                            
                if($major{"$sample,$tissue,$i"} ne $major{"$int,$tissue,$i"}){
                    $total_diff++;
                    if ($minor{"$sample,$tissue,$i"} eq $major{"$int,$tissue,$i"}){
                        $found++;
                        $fre_sum+=$fre{"$sample,$tissue,$i"};
                        $significant++ if(grep {$_ eq $i} @hetero);
                    }
                    
                }
            }
            my $frequency=($found==0)?0:$fre_sum/$found;
            if (($total_diff!=0) || ($found!=0) || ($significant!=0) || ($frequency!=0)) {
                print OUT $sample,"\t",$tissue,"\t",$#hetero+1,"\t","$int\t","$total_diff\t","$found\t","$significant\t",$frequency,"\n";
            }
        }
        close OUT;
    }
}
closedir DH;
