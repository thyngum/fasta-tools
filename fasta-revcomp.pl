#!/usr/bin/perl

# Converts each sequence in a multi-FASTA file to its reverse, complement or 
# reverse-complement counterpart.

# 	fasta-revcomp.pl [-r] [-c] -in <filename> [-out <output>]

# 		-in 	Input multi-FASTA file
#		-out 	Ouput filename (outputs to SDTOUT by default)
# 		-r   	Reverse the sequence
# 		-c   	Convert sequence to its to complement

# Script written by Alejandro Llanes (thyngum@gmail.com)

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use utils;

GetOptions ( 'in=s' => \$input, 'r' => \$reverse, 'c' => \$complement, 'out=s' => \$output );

die "Usage: fasta-revcomp.pl [-r] [-c] -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

if ( ( ! $reverse ) && ( ! $complement ) ) {
	die "Nothing to do, use -r to reverse the sequence and/or -c to get its complement!\n";
}

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

	$tmp = $seq->seq();
	if ( $complement ) {
		$tmp = complement $tmp;
	}
	if ( $reverse ) {
		$tmp = reverse $tmp;
	}
	$seq->seq($tmp);

	if ( $output ) {
		$seqio_out->write_seq($seq);	
	}
	else {
		print ">" . $seq->display_id() . " " . $seq->desc() . "\n";
		print sblock($seq->seq());
	}	
}

