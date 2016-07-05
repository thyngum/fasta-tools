#!/usr/bin/perl

# Translates the sequences in a multi-FASTA sequence file. Designed for a file with  
# coding (spliced) sequences of genes (CDSs).

#   fasta-translate.pl -s -in <filename> [-out <output>]

# 	-in 	Input multi-FASTA file
#	-out 	Ouput filename (outputs to SDTOUT by default)

# Script written by Alejandro Llanes (thyngum@gmail.com)

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output );

die "Usage: fasta-translate.pl -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my $seqio_in = Bio::SeqIO->new(-file => $filename, 
                               -format => 'fasta' ) or die $!;

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
                                 -format => 'fasta');	
}

$count = 0;
while ( $seq = $seqio_in->next_seq ) {
	$count++;

	die "This script can only be used with DNA sequences!\n" unless ( $seq->alphabet() eq 'dna' );

	$tmp_seq = $seq->translate();
	$seq->seq($tmp_seq->seq());

	if ( $output ) {
		$seqio_out->write_seq($seq);	
	}
	else {
		print ">" . $seq->display_id() . " " . $seq->desc() . "\n";
		print sblock($seq->seq());
	}	
}	
