#!/usr/bin/perl
########################################################################
# Heteroplasmy liver lobes project:
# Identify heteroplasmies which might be contaminated by samples that
# were processed at the same time in the wet lab
#
# Copyright (c) 2017 Mingkun Li
########################################################################
use strict;

my $file=shift||die "perl mix_in_batch_parse.pl [mix file _01_70]\n";

my $prop=shift||0.8;
my $level=shift||0.01;
my $level2=shift||0.6;

my @flag;

my %hash;
open IN,$file||$!;
while(<IN>){
      chomp;
      my @array=split /\t/,$_;
      my @arr=split /_/,$array[0];
      my $rec=$arr[0];

      
      my @arr=split /_/,$array[3];
      my $don=$arr[0];
      
      my $fre;
      if ($array[4]){
	    $fre=$array[5]/$array[4]; 
      }
      else{
	    $fre=0;
      }
      $hash{"$rec,$don"}=$_;
      
      my $lev;
      if($array[2]){
	    $lev=$array[6]/$array[2];
      }
      else{
	    $lev=0;
      }
      
      if($fre>=$prop and $array[7]>=$level and $lev>=$level2){
	    push @flag,"$rec,$don" unless (grep {$_ eq "$rec,$don"} @flag);
      }
      
}
close IN;

foreach my $int (@flag){
      my @array=split /\,/,$int;
      print $hash{"$int"},"\n";
}
