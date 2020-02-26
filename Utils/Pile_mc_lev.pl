#!/usr/bin/env perl
#########################################################################
# File Name: pile_mc_lev.pl
# Author: Zhenghena
# mail: zhanghena@gmail.com
# Created Time: 25 May 2018 09:25:43 AM CST
#########################################################################
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
#use Bio::DB::Sam;

my $version='1.00';
sub usage{
    print <<"USAGE";
    Version $version
Usage:  perl $0 -i -t -r --rlength
options:
        -i              input bam file
        -t              input "watson" or "crick", no other options
        -r              raw reference directory
        --region        region for samtools -r
        --list          list file for samtools -l
        --rlength       reads lengh of bam file
        --overhang      delete the bases on the ends(default >=1: all)
        --reads         more than n reads(default >=10)
        --variants      more than n variants(default >=2)
        --minBQ         minimum base quality to filter reads (default >=0)
        --cRatio        methylation ratio(default >=0.05)
        --depth         max per-BAM depth to avoid excessive memory usage (default 80000)
        --all           No cutoff(overhang=1,reads=1,cRatio=0,var=0), 
                        suppress all other parameters.

        -h help
        ex:
            perl rRNAs_pile_mc_lev.pl -i TMEM160.bam -t crick -r /picb/rnomics3/Zhanghena/project/m5C/annotation/combine_genome_ERCC_Rluc/Combined_hg38_ERCC_Rluc.fa --rlength 101
USAGE
exit(1);
}
####################### get options to overlap default values;
my ($in,$t,$ref,$region,$list,$rlength,$over,$reads,$variants,$minBQ,$cRatio,$depth,$all,$help);
$over=1; $variants=2;$minBQ=0; $cRatio=0.05;$reads=1;$depth=80000;
GetOptions(
            'i=s'           =>  \$in,
            't=s'           =>  \$t,
            'r=s'           =>  \$ref,
            'region=s'      =>  \$region,
            'list=s'        =>  \$list,
            'rlength=s'     =>  \$rlength,
            'overhang=i'    =>  \$over,
            'reads=i'       =>  \$reads,
            'variants=i'    =>  \$variants,
            'minBQ=i'       =>  \$minBQ,
            'cRatio=f'      =>  \$cRatio,
            'depth=i'       =>  \$depth,
            'all'           =>  \$all,
            'h'             =>  \$help,
            )or usage();
usage() if $help or !$in or !$t or !$ref or !$rlength;
if($all){
    $over=1; $variants=0; $cRatio=0;$reads=1;
}
my $paras="";
$paras.="-r $region " if $region;
$paras.="-l $list " if $list;

open my $pile,"samtools mpileup -f $ref $in -B -OI --output-QNAME -Q 0 -d $depth $paras|" or die;

##paired reads 会自动检测，overlap的部分的base quality为0
#open my $pile,"sambamba mpileup -t 16 $in --tmpdir=Temp --samtools -f $ref -B -OI --output-QNAME -Q 0 -d $depth $paras |" or die;
##samtools 1.6
#################################################################################################
#ÕýÁ´´¦Àí
goto L if $t !~ /watson/i;
my %idx=(A=>7,C=>9,G=>11,N=>13,T=>15);

#rDNA    9       G       5       ^].^].^].^].^]. DDDDD   
#rDNA    10      C       7       TTTTT^]t^]t     DDDBDII 
my $aa;
while ($aa=<$pile>){
    chomp $aa;
    my @F=split( /\t/, $aa);
    $F[2] = uc($F[2]);
    next unless $F[2] =~ /^[Cc]$/ and $F[3] > 0;
    $F[4]=~s/(\^.)+|(\$)+//g;   ### remove starts, qualities, and ends;
    my $temp=$F[4];             ### remove insertions and deletions.    ### updated by sszhu1007@gmail.com 13.06.13 18:49:05 ###
    while($temp=~/(\+|-)(\d+)/g){
        my $ss=$1;
        $ss="\\$ss" if $ss eq "+";
        $F[4]=~s/$ss$2 @{["." x $2]}//x;
        }
    my @a=split //,$F[4];
    ##########################################
    my @positions=map{chomp;$_} split /,/,$F[6];
    my @scores=split //,$F[5];
    my @p_sites=grep {$a[$_]!~/[><\*]/ and $positions[$_]>=$over and $positions[$_]<=$rlength-$over+1 } 0..$#positions; 
    #my @p_sites=grep {$a[$_]!~/[><\*]/} 0..$#positions; 
        my $ct_total=@p_sites;     #### contain the minBQ filtered bases;
        if($minBQ){
            @p_sites=grep {ord($scores[$_])-33>=$minBQ} @p_sites;
        }
    my $ct=@p_sites; 
    next if $ct <1;
    next unless !$reads or ($reads and $ct_total>=$reads);           #### 1 reads
   
    my @a_new=@a[@p_sites];
    my @positions_new=@positions[@p_sites];
    my @scores_new=@scores[@p_sites];

    my (%ct,$vct,@locs,@ids); $vct=0;
    map{$ct{$_}=0} qw(A C G N T);
    map{$ct{$a_new[$_]=~/[,.]/?uc($F[2]):($vct++,(push @locs,$_),uc $a_new[$_])}++} 0..$#a_new;
    next unless !$variants or ( $variants and $ct{"C"}>=$variants);      #### number of methylated Cs
    next unless !$cRatio or ($cRatio and $ct{"C"}/($ct{"C"}+$ct{"T"}) > $cRatio); ## methylation level
    next if $ct{"C"}+$ct{"T"}<1; ### 
    print join("\t","$F[0]:$F[1]","C","+",$ct,$ct{"C"},$ct{"T"},$ct{"C"}/($ct{"C"}+$ct{"T"}),($ct{"C"}+$ct{"T"})/$ct,(map {$_?$_:0} map { $ct{$_}} qw(A C G N T)) ,$ct_total)."\n";  ### updated by zhanghena@gmail.com 18.05.07 ###
}
close IN;
exit (0);
###################################################################################
L:

#rDNA    9       G       5       ^].^].^].^].^]. DDDDD   
#rDNA    10      C       7       TTTTT^]t^]t     DDDBDII 

while ($aa=<$pile>){
    chomp $aa;
    my @F=split( /\t/, $aa);
    $F[2] = uc($F[2]);
    next unless $F[2] =~ /^[Gg]$/ and $F[3] > 0;
    $F[4]=~s/(\^.)+|(\$)+//g;   ### remove starts, qualities, and ends;
    my $temp=$F[4];             ### remove insertions and deletions.    ### updated by sszhu1007@gmail.com 13.06.13 18:49:05 ###
    while($temp=~/(\+|-)(\d+)/g){
        my $ss=$1;
        $ss="\\$ss" if $ss eq "+";
        $F[4]=~s/$ss$2 @{["." x $2]}//x;
        }
    my @a=split //,$F[4];
    my @positions=map{chomp;$_} split /,/,$F[6];
    my @scores=split //,$F[5];
    my @p_sites=grep {$a[$_]!~/[><\*]/ and $positions[$_]>=$over and $positions[$_]<=$rlength-$over+1 } 0..$#positions; 
        my $ct_total=@p_sites;     #### contain the minBQ filtered bases;
        if($minBQ){
            @p_sites=grep {ord($scores[$_])-33>=$minBQ} @p_sites;
        }
    my $ct=@p_sites; 
    
    next if $ct <1;
    next unless !$reads or ($reads and $ct_total>=$reads);           #### 1 reads
   
    my @a_new=@a[@p_sites];
    my @positions_new=@positions[@p_sites];
    my @scores_new=@scores[@p_sites];

    my (%ct,$vct,@locs,@ids); $vct=0;
    map{$ct{$_}=0} qw(A C G N T);
    map{$ct{$a_new[$_]=~/[,.]/?uc($F[2]):($vct++,(push @locs,$_),uc $a_new[$_])}++} 0..$#a_new;
    next unless !$variants or ( $variants and $ct{"G"}>=$variants);      #### number of methylated Cs
    next if $ct{"G"}+$ct{"A"}<1; ### 
    next unless !$cRatio or ($cRatio and $ct{"G"}/($ct{"G"}+$ct{"A"}) > $cRatio); ## methylation level    
    print join("\t","$F[0]:$F[1]","C","-",$ct,$ct{"G"},$ct{"A"},$ct{"G"}/($ct{"G"}+$ct{"A"}),($ct{"G"}+$ct{"A"})/$ct,(map {$_?$_:0} map { $ct{$_}} qw(A C G N T)) ,$ct_total)."\n";  ### updated by zhanghena@gmail.com 18.05.07 ###
}
close IN;

#*****************************************************************************************************

sub watson_mC_type{
	my @base = @_ ;
	
	if ($base[0] =~ /[Gg]/) {
		return "CpG";
	} elsif ($base[0] =~ /[AaCcTt]/) {
		if ($base[1] =~ /[Gg]/) {
			return "CHG";
		} elsif ($base[1] =~ /[AaCcTt]/) {
			return "CHH";
		}
	}

	return "NA";
}

sub crick_mC_type{
	my @base = @_ ;
	if ($base[0] =~ /[Cc]/) {
		return "CpG";
	} elsif ($base[0] =~ /[AaGgTt]/) {
		if ($base[1] =~ /[Cc]/) {
			return "CHG";
		} elsif ($base[1] =~ /[AaGgTt]/) {
			return "CHH";
		}
	}

	return "NA";
}
