#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use Bio::Cigar;  # Cpan this if you want to run this script !!
use List::Util qw(max);

use vars qw/ %opt /;

# this program is subsequently packed via App::FatPacker
# see http://perltricks.com/article/58/2014/1/5/The-easy-way-to-build-stand-alone-Perl-apps/
# install fatpack via cpan App::FatPacker
#
# fatpack file hapture.pl > hapture
# either perl hapture or -run as executable by chmod 755


sub init(){
	getopts( "hv:s:i:g:", \%opt ) or usage();
	usage() if $opt{h};
	print STDERR "Require path specification for VCF file: -v\n" and exit if not defined $opt{v};
	print STDERR "Require path specification for SAM file: -s\n" and exit if not defined $opt{s};
	print STDERR "Require individual ID: -i\n" and exit if not defined $opt{i};
	print STDERR "Require group ID: -g \n" and exit if not defined $opt{g};
}


sub usage(){
print STDERR << "EOF";
    This program gathers variant haplotype sites from SAM alignment file, and reports a summary file for those variant sites

    usage: $0 [-h] -v vcf_file -s sam_file -i int

     -h        : this (help) message
     -v file   : variant caller file - VCF format (!! assumed the position is sorted)
     -s file   : sequence alignment file - SAM format
     -i int    : individual ID (integer value or unbroken string)
     -g str    : group ID (unbroken string)

    example: $0 -v s1.vcf -s s1.sam -i 0 -g sebastes

EOF
        exit;
}

#----- get user's input parameters ------

init();

#---------------------------------------


#----- read vcf file -------------

# Objective: keeps variants' info into memory so that
# I can tell whether there are any variant sites for
# the alignment read entry from the SAM file

my $vcf; # a hash reference that keeps track of essential vcf info: reference variant, derived variants, pos
my $hap;

open VCF, $opt{v};
while (<VCF>) {
	next if /^#/;
	my @line = split "\t";
	# for now, we only concern SNP site, (exclude indel and any other complex events)

	# reference allele
	next if length($line[3]) > 1;
	# derived allele
	my $max_len = 0;
	my @snp = split ",", $line[4];
	for my $deriv (@snp) {
		$max_len = length($deriv) if length($deriv) > $max_len;
	}
	next if $max_len > 1;

	#print $line[0], "\t", $line[3], "\t", $line[4], "\n";
	push @{$vcf->{$line[0]}}, $line[1];
	push @{$vcf->{"ref_".$line[0]}}, $line[3];
	push @{$vcf->{"der_".$line[0]}}, $snp[0];
}
close VCF;

#----------------------------------------


open SAM, $opt{s};
while(<SAM>) {
	next if /^\@/;
	my @lines = split "\t";
	my $id = $lines[2];
    next if $lines[1] >= 256; # skip entries that are secondary alignment or to multiple sites
	my $st_qpos = $lines[3]; # starting query position
	#my $mapq = $lines[4]; # mapping quality score
	#skip if the alignment id is not found in the vcf hash ref
	next if not defined $vcf->{$id};

  next if $lines[5] eq "*";
	my $cigar = Bio::Cigar->new($lines[5]);
	my @qseq = split "", $lines[9];
	my @qseq_qual = split "", $lines[10];

	next if $#qseq < 1;
	#print "Query length is ", $cigar->query_length, "\n";
	#print "Reference length is ", $cigar->reference_length, "\n";
	my $hapRead={};
	$hapRead->{"seq"} = "";
	my $ct=0;
	for my $rpos (@{$vcf->{$id}}) {
		my $rpos_adj = $rpos - $st_qpos +1;
		$ct++;

		if ($cigar->reference_length < $rpos_adj || $rpos_adj < 1) {
			$hapRead->{"seq"} .= "N"; # marked as unknown
			push @{$hapRead->{"qual"}}, "_";
			next;
		}
		my ($qpos, $op) = $cigar->rpos_to_qpos($rpos_adj);

		if (not defined $qpos) {
			$hapRead->{"seq"} .= "X"; # marked as deletion site
			push @{$hapRead->{"qual"}}, "_";
 		}
		else {
			$hapRead->{"seq"} .= $qseq[$qpos-1];
			push @{$hapRead->{"qual"}}, $qseq_qual[$qpos-1];
		}

		#print $qpos-1, "\t", $#qseq, "\t", $lines[0], "\t", $id, "\n";
		#print join "\t", $id, $rpos, $qseq[$qpos-1], $qseq_qual[$qpos-1], "\n" if $qpos != -1;
	}

	$hap->{$id}->{$hapRead->{"seq"}}->{"ct"}++;
	for my $i (0..$#{$vcf->{$id}}) {

	  if(! defined ${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i]) {
	 	${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i] = 0;
	 	${$hap->{$id}->{$hapRead->{"seq"}}->{"sC"}}[$i] = 0;
	 	}

		my $q = 10**(-(ord(${$hapRead->{"qual"}}[$i])-33)/10);
		#${$hap->{$id}->{$hapRead->{"seq"}}->{"logC"}}[$i]+= log(1-$q) ;
	 	#${$hap->{$id}->{$hapRead->{"seq"}}->{"logW"}}[$i]+= log($q);
	 	${$hap->{$id}->{$hapRead->{"seq"}}->{"sC"}}[$i]+= 1-$q; # collecting the sum of prob phred site score for future use
	 	${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i] = max(1-$q, ${$hap->{$id}->{$hapRead->{"seq"}}->{"maxC"}}[$i]);

	}
	#${$hap->{$id}->{$hapRead->{"seq"}}->{"mapq"}->{$mapq}}++; # collecting the MapQ alignment score

}


#--- output a haplotype summary file -----

for my $id (keys %{$hap}){
	for my $h (keys %{$hap->{$id}}){

		print join "\t", $opt{g}, # group label
				 $opt{i}, # individual id label
				$id, #locus id
				$h, # haplotype sequence
				$hap->{$id}->{$h}->{"ct"},  # number of occurence observed for this haplotype or read depth
				#(join ",", @{$hap->{$id}->{$h}->{"logC"}}), # log scale phred stat for being a correct base
				#(join ",", @{$hap->{$id}->{$h}->{"logW"}}), # log scale phred stat for being a miscalled base
				#(join ",", @{$vcf->{$id}}); # variant position
				(join ",", @{$hap->{$id}->{$h}->{"sC"}}), # sum of prob for calling correct base
				(join ",", @{$hap->{$id}->{$h}->{"maxC"}})."\n";# max phred score

		#print "\t";

    #print join ",", map {$_, ${$hap->{$id}->{$h}->{"mapq"}->{$_}}} (keys %{$hap->{$id}->{$h}->{"mapq"}});
	  #print "\n";
	}
}

#----------------------------------------
