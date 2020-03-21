#!/usr/bin/perl
 
# Convert a pair of FASTA/QUAL files with the same basename to an
# FASTQ file. Output to STDOUT.

# fasta-qual-to-fastq [-f <format>] -in <filename>

# -f <format>	Output file format can be fastq (default), fastq-sanger, 
#		fastq-illumina or fastq-solexa.

# <filename> is the name of the FASTA file, including the suffix,
# i.e. if the FASTA file is 'file.fna', the QUAL file should be
# 'file.qual'. 

# TODO: print to output file. 
 
use Bio::SeqIO;
use Bio::Seq::Quality;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'in=s' => \$input, 'f=s' => \$format );

die "Usage: fasta-qual-to-fastq.pl -in <filename>\n\n",
    "<filename> is the name of the input FASTA file, including the\n",
    "suffix. There should be a QUAL file with the same basename.\n\n" if ( ! $input );

$format = 'fastq' if ( ! $format );

$filename = File::Spec->rel2abs($input);
die "FASTA file \'$filename\' not found.\n" if ( ! -e $filename );
 
($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
$qualfile = $path . $name . ".qual";
die "QUAL file \'$qualfile\' not found.\n" if ( ! -e $qualfile );

$in_seq_obj = Bio::SeqIO->new( -file   => $filename,
		                       -format => 'fasta', );
 
$in_qual_obj = Bio::SeqIO->new( -file   => $qualfile,
		                        -format => 'qual', );
 
$out_fastq_obj = Bio::SeqIO->new( -format => $format );
 
while (1) {
	# Create objects for both a seq and its associated qual
	$seq_obj  = $in_seq_obj->next_seq || last;
	$qual_obj = $in_qual_obj->next_seq;
 
	die $! unless $seq_obj->id eq $qual_obj->id;
 
	# Use seq and qual object methods feed info for new BSQ object.
	$bsq_obj = Bio::Seq::Quality->new( -id   => $seq_obj->id,
	  	                               -seq  => $seq_obj->seq,
			                           -qual => $qual_obj->qual,
			                         );
	 
	# Print the output.
	$out_fastq_obj->write_fastq($bsq_obj);
}

