#!/usr/bin/perl

# Converts a FASTQ file to FASTA format, discarding quality info. 

#   fastq-to-fasta.pl -in <filename> [-out <output>]

# 	-in 	Input FASTQ file
#	-out 	Ouput filename (outputs to SDTOUT by default)

# Script written by Alejandro Llanes (thyngum@gmail.com)

use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output );

die "Usage: fastq-to-fasta.pl -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	open OUT, ">$output_filename";
}

open IN, $filename or die "File '$filename' not found.\n";

my $line_count = 0;
while ( my $line = <IN> )  {  
	my $desc = $line;
	my $seq = <IN>;
	my $plus = <IN>;
	my $qual = <IN>;

	$line_count = $line_count + 4;

	die "Incorrect FASTQ format around line $line_count\n" unless ( $desc =~ /^\@/ and $plus =~ /^\+/ ); 

	chomp $desc;
	chomp $seq;
	chomp $plus;
	chomp $qual;

	( my $id ) = split /\s/, $desc;
	$id = substr $id, 1;

	if ( $output ) {
		print OUT ">$id\n";
		print OUT sblock($seq);	
	}
	else {
		print ">$id\n";
		print sblock($seq);
	}

}	
close IN;

if ( $output ) {
	close OUT;
}

$count = $line_count/4;

unless ( $output ) {
	print "\n";
}

print STDERR "Found $count sequences in input file.\n";
