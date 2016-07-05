#!/usr/bin/perl

# Converts a FASTQ file to an equivalent pair of FASTA/QUAL files with the same 
# basename of input file.

#   fastq-to-fasta-qual [-f <format>] -in <filename>

# 	-in	Input FASTQ file
# 	-f 	Input file format: fastq (default), fastq-sanger, fastq-illumina or fastq-solexa.

use Bio::SeqIO;
use Bio::Seq::Quality;
use File::Basename;
use File::Spec;
use Getopt::Long;
 
GetOptions ( 'in=s' => \$input, 'f=s' => \$format );

die "Usage: fastq-to-fasta-qual.pl [-f <format>] -in <filename>\n" if ( ! $input );

$format = 'fastq' if ( ! $format );

$filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );
 
($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
$fasta_file = $path . $name . ".fasta";
$qual_file = $path . $name . ".qual";
 
$in_seq_obj = Bio::SeqIO->new(-file   => $filename,
                              -format => $format, );

$out_fasta_obj = Bio::SeqIO->new(-file   => ">$fasta_file",
                                 -format => 'fasta', );

$out_qual_obj = Bio::SeqIO->new(-file   => ">$qual_file",
                                -format => 'qual', );

while ( $seq = $in_seq_obj->next_seq ) {
	$out_fasta_obj->write_seq($seq);
	$out_qual_obj->write_seq($seq);
	$count++;
}
