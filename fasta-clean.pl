#!/usr/bin/perl

# Clean the headers of the sequences in a multi-FASTA file leaving only primary IDs.

#   fasta-clean.pl -in <filename> [-out <output>] 

# 	-in 	Input multi-FASTA file
#	-out 	Ouput filename (outputs to SDTOUT by default)

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output );

die "Usage: fasta-clean.pl -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my $seqio_in = Bio::SeqIO->new(-file => $filename, 
                                   -format => 'fasta' ) or die $!;

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
                                 -format => 'fasta');	
}

while ( my $seq = $seqio_in->next_seq ) {

	my $tmp = $seq->display_id();
	$seq->desc("");
	$seq->display_id($tmp);

	if ( $output ) {
		$seqio_out->write_seq($seq);	
	}
	else {
		print ">" . $seq->display_id() . "\n";
		print sblock($seq->seq());
	}		
}
