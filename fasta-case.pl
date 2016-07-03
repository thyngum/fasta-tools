#!/usr/bin/perl

# Change the case of sequences in a multi-FASTA file.

#   fasta-case.pl <-u/l> -in <filename> [-out <output>]

# 	-in 	Input multi-FASTA file
#	-out 	Ouput filename (outputs to SDTOUT by default)
# 	-u/l	Change to uppercase/lowercase. 

# Script written by Alejandro Llanes (thyngum@gmail.com)

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use utils;

GetOptions ( 'in=s' => \$input, 'u' => \$uppercase, 'l' => \$lowercase, 'out=s' => \$output );

die "Usage: fasta-case.pl <-u/l> [-out <output>] -in <filename>\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

if ( ( ! $uppercase ) && ( ! $lowercase ) ) {
	die "Nothing to do, use -u for uppercase or -l for lowercase!\n";
}

my $seqio_in = Bio::SeqIO->new(-file => $filename, 
                                   -format => 'fasta' ) or die $!;

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
                                 -format => 'fasta');	
}

while ( my $seq = $seqio_in->next_seq ) {
	if ( $uppercase ) {
		$tmp = uc($seq->seq);
	}
	else {
		$tmp = lc($seq->seq);
	}
	$seq->seq($tmp) ;
	
	if ( $output ) {
		$seqio_out->write_seq($seq);	
	}
	else {
		print ">" . $seq->display_id() . " " . $seq->desc() . "\n";
		print sblock($seq->seq());
	}	
}
