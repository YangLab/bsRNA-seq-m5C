#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;
use Getopt::Long;

#! usr/bin/perl -w
use strict;
system "echo `date` dir: `pwd` command:'$0 @ARGV' >>$0.arg";

# Store command line arguements
 my $genomeDir = '';
 my $genePred = '';	
 GetOptions ('genomeDir=s' => \$genomeDir, 'genePred=s' => \$genePred); #maybe add option for single multifasta file

# Check command line args are minimally valid
&checkArgs($genomeDir, $genePred) || die "Check command line arguments: $!";

# Annotate every genic nucleotide
open (ANNOT, "$genePred")  || die "Cant open genePred file: $!";

my $counter = 0;
my $header = <ANNOT>;
my $chrom_seq;
my $last_chr = '';

while (my $line = <ANNOT>){
	chomp ($line);

	# Split line and assign variables
	my @temp = split(/\t/, $line);
	my ($transcript_ID, $chr, $strand) = @temp[1..3];
	my ($txStart, $txEnd, $cdsStart, $cdsEnd) = @temp[4..7];
	my $exonCount = $temp[8];
	my @exonFrame = split(/,/, $temp[-1]);
	my @exon_starts = split (/,/, $temp[9]);
	my @exon_ends = split (/,/, $temp[10]);
	my $gene_name = $temp[12];
##################################################added intron information by zhn  
    my $intronCount = $temp[8] - 1;
    my @intron_starts = split (/,/, $temp[10]);
    splice @intron_starts, -1, 1;
    my @intron_ends = split (/,/, $temp[9]);
    splice @intron_ends, 0, 1;
    my $pre_size =  abs($txEnd - $txStart);
###################################################
    
	# Load genomic sequence of chromosome (if necessary)
	if ($chr ne $last_chr){
		$chrom_seq = &loadGenome($genomeDir, $chr);
		$last_chr = $chr;
	}

	# Determine size of transcript
	my $gene_size = 0;
	my $num_exons = scalar (@exon_starts);
	for (my $i = 0; $i < $num_exons; $i++){
		$gene_size += abs ($exon_ends[$i] - $exon_starts[$i]);
	}
######################################################## Determine size of CDS by zhn
    my $cds_size = 0;
    my $UTR5_size = 0;
    my $UTR3_size = 0;
    my $utr5_st   = "NA";
    my $utr5_end  = "NA";
    my $cds_st    = "NA";
    my $cds_end   = "NA";
    my $utr3_st   = "NA";
    my $utr3_end  = "NA";
    if( abs($cdsStart - $cdsEnd) > 0){
        for (my $i = 0; $i < $num_exons; $i++){
            if ($exon_starts[$i] <= $cdsStart && $exon_ends[$i] >= $cdsStart && $exon_ends[$i] <= $cdsEnd) {
                # 1. #####|###   |
                $cds_size += abs ($exon_ends[$i] - $cdsStart);
                if ($strand eq "+"){
                    $UTR5_size += abs ($cdsStart - $exon_starts[$i]);
                }
                else{
                    $UTR3_size += abs ($cdsStart - $exon_starts[$i]);
                }
            }
            elsif ($exon_starts[$i] <= $cdsEnd && $exon_ends[$i] >= $cdsEnd && $exon_starts[$i] >= $cdsStart){ 
                #2. |    ####|####
                $cds_size += abs ($cdsEnd - $exon_starts[$i]);
                if ($strand eq "+"){
                    $UTR3_size += abs ($exon_ends[$i] - $cdsEnd);
                }
                else{
                    $UTR5_size += abs ($exon_ends[$i] - $cdsEnd);
                }
            }
            elsif ($exon_starts[$i] >= $cdsStart && $exon_ends[$i] <= $cdsEnd ){ #3. |  ######  |
                $cds_size += abs ($exon_ends[$i] - $exon_starts[$i]);
            }
            elsif ($exon_starts[$i] <= $cdsStart && $exon_ends[$i] <= $cdsStart ){ #4. ######  |        |
                if ($strand eq "+"){
                    $UTR5_size += abs ($exon_ends[$i] - $exon_starts[$i]);
                }
                else{
                    $UTR3_size += abs ($exon_ends[$i] - $exon_starts[$i]);
                }
            }
            elsif ($exon_starts[$i] >= $cdsEnd && $exon_ends[$i] >= $cdsEnd ){ #5. |     |    #######
                if ($strand eq "+"){
                    $UTR3_size += abs ($exon_ends[$i] - $exon_starts[$i]);
                }
                else{
                    $UTR5_size += abs ($exon_ends[$i] - $exon_starts[$i]);
                }
            }
            elsif ($exon_starts[$i] <= $cdsStart && $exon_ends[$i] >= $cdsEnd ){ #6. ####|#####|#### 
                $cds_size += abs ($cdsEnd - $cdsStart);
                if ($strand eq "+"){
                    $UTR5_size += abs ($cdsStart - $exon_starts[$i]);
                    $UTR3_size += abs ($exon_ends[$i] - $cdsEnd);
                }
                else{
                    $UTR5_size += abs ($exon_ends[$i] - $cdsEnd);
                    $UTR3_size += abs ($cdsStart - $exon_starts[$i]);
                }
            }
        }
    $utr5_st   = 1;
    $utr5_end  = $utr5_st + $UTR5_size - 1;
    $cds_st    = $utr5_end + 1;
    $cds_end   = $cds_st + $cds_size - 1;
    $utr3_st   = $cds_end + 1;
    $utr3_end  = $gene_size;
    #print $gene_name, "\t", $strand,"\t", $gene_size,"\t", $UTR5_size, "\t",$cds_size,"\t", $UTR3_size, "\t";
    }    

########################################################
	# Conversion from 0-based to 1-based coordinates
	$txStart = $txStart + 1;
	$cdsStart = $cdsStart + 1;
	foreach my $x (@exon_starts){
		$x = $x + 1
	}
    foreach my $x (@intron_starts){
		$x = $x + 1
	} ###zhn

	# Label the exonic loci
	my $nucl_j = -1;
	my $count = 0;
    my $cds_count = 0;
	my @temp_out = ();
########################################################### added loci in pre_mRNA and splice information by zhn
	for (my $i = 0; $i < $exonCount; $i++){ # Go through each exon delimiter
		foreach my $coord (($exon_starts[$i])..($exon_ends[$i])){ # Go through each position of the exon
			$count++;
			my $mrna_pos = &determine_mrna_pos($strand, $gene_size, $count);
            my $pre_pos = &determine_pre_pos($strand, $coord, $txStart, $txEnd,); ###zhn
			my $feature = &utrs_or_cds($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd);
            my $cds_pos = "NA";
            if ($feature eq "cds"){
                    $cds_count++;
                    $cds_pos = &determine_cds_pos($strand, $cds_size, $cds_count);                
            }
			my $nucleotide = uc(substr ($chrom_seq, $coord-1, 1));
            my $rel_pos = "NA";
            $rel_pos = &determine_rel_pos($mrna_pos, $feature, $UTR5_size, $cds_size, $UTR3_size);
			#my $codon_pos = &determine_codon_pos($feature, $nucl_j, $strand); 
			if (($nucleotide eq "C" && $strand eq "+") || ($nucleotide eq "G" && $strand eq "-")){
				print $chr, "\t", $coord-1, "\t", $coord, "\t", join('|', ($transcript_ID, $gene_name, $feature, $mrna_pos, $pre_pos,$cds_pos,join("-",$exon_starts[$i],$exon_ends[$i]))), "\t";
				print $nucleotide, "\t", $strand,"\t", join(':',$chr,$coord,"C",$strand),"\t",$gene_size,"\t", $UTR5_size, "\t",$cds_size,"\t", $UTR3_size, "\t";
            	print join("\t",$utr5_st,$utr5_end,$cds_st,$cds_end,$utr3_st,$utr3_end, $rel_pos),"\t";
            
            	my $dis_donor = "NA";
            	my $dis_acceptor = "NA";
            	if($i == 0){
                	$dis_donor = $coord - $exon_ends[$i]; 
            	}
            	elsif ($i == $exonCount - 1){
                	$dis_acceptor = $coord - $exon_starts[$i];
            	}
            	else{
                	$dis_acceptor  = $coord - $exon_starts[$i];
             	   $dis_donor     = $coord - $exon_ends[$i];
            	}
           		if ($strand eq "+ ") {
            		print $dis_donor,"\t",$dis_acceptor,"\n";
          		}
          		else{
              		print $dis_acceptor,"\t",$dis_donor,"\n";
         		}
			}
		}
	}
######################################################### added intronic and loci in pre_mRNA and splice information by zhn
    for (my $i = 0; $i < $intronCount; $i++){ # Go through each exon delimiter
		foreach my $coord (($intron_starts[$i])..($intron_ends[$i])){ # Go through each position of the exon
			#$count++;
			my $mrna_pos = "NA";
            my $pre_pos = &determine_pre_pos($strand, $coord, $txStart, $txEnd,);
			my $feature = &intronic($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd);
			my $nucleotide = uc(substr ($chrom_seq, $coord-1, 1));
			#my $codon_pos = &determine_codon_pos($feature, $nucl_j, $strand);
			if (($nucleotide eq "C" && $strand eq "+") || ($nucleotide eq "G" && $strand eq "-")){
				print $chr, "\t", $coord-1, "\t", $coord, "\t", join('|', ($transcript_ID, $gene_name, $feature, $mrna_pos, $pre_pos,"NA",join("-",$intron_starts[$i],$intron_ends[$i]))), "\t";
				print $nucleotide, "\t", $strand,"\t",join(':',$chr,$coord,"C",$strand),"\t", $gene_size,"\t", $UTR5_size, "\t",$cds_size,"\t", $UTR3_size, "\t"; 
            	print join("\t",$utr5_st,$utr5_end,$cds_st,$cds_end,$utr3_st,$utr3_end, "NA"),"\t";
            
            	my $dis_acceptor  = $coord - $intron_ends[$i];
            	my $dis_donor     = $coord - $intron_starts[$i];
            	if ($strand eq "+ ") {
                	print $dis_donor,"\t",$dis_acceptor,"\n";
            	}
            	else{
                	print $dis_acceptor,"\t",$dis_donor,"\n";
            	}
			}
        }
    }
}

######################### sub routunes ##########################
sub determine_codon_pos(){
	my ($feature, $nucl_j, $strand) = @_;
	my $codon_pos = 'NA';

	if ($feature eq 'cds'){
		$nucl_j++;
		$codon_pos = $nucl_j % 3;

		# correct $codon_pos for strandedness
		if ($strand eq '-'){
			if ($codon_pos == 0){ $codon_pos = 2; }
			elsif ($codon_pos == 2){ $codon_pos = 0; }
		}
	}
	return $codon_pos;
}

sub determine_mrna_pos(){
	my ($strand, $gene_size_ref, $count) = @_;
	my $mrna_pos = $count;
	if ($strand eq '-'){
		$mrna_pos = $gene_size_ref - $count + 1 ;;
	}  
	return $mrna_pos;
}

sub determine_cds_pos(){
	my ($strand, $cds_size, $cds_count) = @_;
	my $cds_pos = $cds_count;
	if ($strand eq '-'){
		$cds_pos = $cds_size - $cds_count;
	}  
	return $cds_pos;
}

sub determine_pre_pos(){
	my ($strand, $coord,$txStart, $txEnd) = @_;
	my $pre_pos = $coord - $txStart + 1;
	if ($strand eq '-'){
		$pre_pos = $txEnd - $coord + 1;
	}  
	return $pre_pos;
}

sub determine_rel_pos(){
    my ($mrna_pos, $feature, $UTR5_size, $cds_size, $UTR3_size) = @_;
    my $rel_pos = 0;
    if ($feature eq "5utr"){
        $rel_pos = $mrna_pos/$UTR5_size;
    }
    if ($feature eq "cds"){
        $rel_pos = 1 + ($mrna_pos - $UTR5_size)/$cds_size;
    }
    if ($feature eq "3utr"){
        $rel_pos = 2 + ($mrna_pos - $UTR5_size - $cds_size)/$UTR3_size;
    }
	return $rel_pos;
}



sub checkArgs(){
	my ($genDir, $genPred) = @_;
	my $pass = 1;
	if (! -d $genDir){
		print "No such directory exists: $genDir \n";
		$pass = 0;
	}
	if (! -f $genPred){
		print "No such file exists: $genDir \n";
		$pass = 0;
	}
	return ($pass);
}

sub loadGenome(){
	my ($genDir, $chr) = @_;
	my $file = $genDir.'/'.$chr.'.fa';
	my $seqio_obj = Bio::SeqIO->new(-file => $file);
	my $seq_obj = $seqio_obj->next_seq;
	return ($seq_obj->seq);
}

sub intronic (){
	my ($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd) = @_;
	my $feature = 'na';
	#my ($utr5, $cds, $utr3) = (0,0,0);

	if( abs($cdsStart - $cdsEnd) <= 1){ ## really if the cdsStart = cdsEnd.  Need to specify this way due to 0 => 1 based conversion.
		$feature = 'ncRNA_intronic';
    }
    else{$feature = 'intronic';}
	return ($feature);
}


sub utrs_or_cds (){
	my ($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd) = @_;
	my $feature = 'na';
	#my ($utr5, $cds, $utr3) = (0,0,0);

	if( abs($cdsStart - $cdsEnd) <= 1){ ## really if the cdsStart = cdsEnd.  Need to specify this way due to 0 => 1 based conversion.
		$feature = 'ncRNA_exonic';
	}
	elsif ($coord >= $txStart && $coord < $cdsStart){
		if ($strand eq '+'){$feature = '5utr';}
		else {$feature = '3utr';}
	}
	elsif ($coord > $cdsEnd && $coord <= $txEnd){
		if ($strand eq '+'){$feature = '3utr';}
		else {$feature = '5utr';}
	}
	else{$feature = 'cds';}

	return ($feature);
}
