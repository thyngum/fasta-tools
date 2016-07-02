#!/usr/bin/perl

# Extracts the sequences from a multi-FASTA file and saves them to individual FASTA 
# files. Output files are named after the corresponding primary IDs.

# 	fasta-split.pl -in <filename>

# 		-in 	Input multi-FASTA file

# Script written by Alejandro Llanes (thyngum@gmail.com)

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'in=s' => \$input );

die "Usage: fasta-split.pl -in <filename>\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my ($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);

my $seqio_input = Bio::SeqIO->new(-file => $filename, 
                                  -format => 'fasta' ) or die $!;

my $count = 0;
while ( my $seq = $seqio_input->next_seq ) {
	my $output = $path . $seq->display_id() . ".fasta";
	print STDERR "$output\n";
	my $seqio_output = Bio::SeqIO->new(-file => ">$output" ,
                                    -format => 'fasta');
	$seqio_output->write_seq($seq);

	$count++;
}

print STDERR "\nOriginal file \'$filename\' split into $count output file(s).\n";
