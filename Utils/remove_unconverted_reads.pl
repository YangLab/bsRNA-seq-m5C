#!/usr/bin/env perl
#########################################################################
# File Name: remove_unconverted_reads.pl
# Author: Zhenghena
# mail: zhanghena@gmail.com
# Created Time: 23 Aug 2018 09:25:43 AM CST
#########################################################################

my $version=1.00;

use strict;
use Getopt::Long;
use warnings;
use File::Basename;
use Bio::DB::Sam;

sub usage{
    print <<"USAGE";
    Version $version
Usage:  perl $0 -i -t
options:
        -i              input bam file
        -t              input "watson" or "crick", no other options
        -n              max number of unconverted Cs in a read.(default <=3)
        -h help
        ex:
            perl filter_unconverted_reads.pl -i *.bam -t crick
USAGE
exit(1);
}
####################### get options to overlap default values;
my ($in,$t,$n,$help);
$n=3;
GetOptions(
            'i=s'           =>  \$in,
            't=s'           =>  \$t,
            'n=i'           =>  \$n,
            'h'             =>  \$help,
            )or usage();
usage() if $help or !$in or !$t;

open my $BAM,"samtools view -h $in |" or die;

#################################################################################################
goto L if $t !~ /watson/i;
my ($aa,$seq,$N);
while (<$BAM>){
    chomp $_;
	if(/^(\@)/){
    print $_."\n";
    }
    else{
        my @F=split( /\t/, $_);
        ($seq,$N)=c2t($F[9]);
		next if($N>$n);
		#print join("\t",$F[0],$_,$N)."\n";
		print join("\t",@F[0..$#F])."\n";
    }
}


close IN;
exit (0);

#################################################################################################
L:

while (<$BAM>){
	chomp $_;
	if(/^(\@)/){
	print $_."\n";
	}
     else{
        my @F=split( /\t/, $_);
        ($seq,$N)=g2a($F[9]);
		next if($N>$n);
		#print join("\t",$F[0],$N)."\n";
		print join("\t",@F[0..$#F])."\n";
    }


}
close IN;
#################################################################################################

sub c2t{
    my ($str) = @_ ;
    my (@arr, @brr, $seq); 
    @arr = split(//, $str);
    my $j=0;
    for (my $i=0; $i<=$#arr; $i++) {
        if ($arr[$i] eq 'C' or $arr[$i] eq 'c') {
            $j++;
            push @brr, 'T';
        }else{
            push @brr, $arr[$i];
        }
    }
    
    $seq = join("", @brr);
    return ($seq,$j);
}

sub g2a{
    my ($str) = @_ ;
    my (@arr, @brr, $seq);  
    @arr = split(//, $str);
    my $j=0;
    for (my $i=0; $i<=$#arr; $i++) {
        if ($arr[$i] eq 'G' or $arr[$i] eq 'g') {
            $j++;
            push @brr, 'A';
        }else{
            push @brr, $arr[$i];
        }
    }
  
    $seq = join("", @brr);
    return ($seq,$j);
}
