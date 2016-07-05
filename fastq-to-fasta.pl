#!/usr/bin/perl

# Converts a FASTQ file to FASTA, discarding quality info. Intended to be used
# as a more efficient alternative to fastq-to-fasta-qual.pl for large FASTQ files.

#   fastq-to-fasta.pl -in <filename> [-out <output>]

# 	-in 	Input FASTQ file
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

die "Usage: fastq-to-fasta.pl -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

open IN, $filename or die("Could not open input file.\n");

my $line_count = 0;
my $count = 0;
while ( my $line = <IN> )  {  
	$line_count++;
	my $start = substr $line, 0, 1;
	my $dummy = substr $line, 0, 2;		
	if ( ( $start eq '@' ) and ( $dummy ne '@@' ) ) {
		my $desc = $line;
		my $seq = <IN>;
		my $plus = <IN>;
		my $qual = <IN>;
		chomp $desc;
		chomp $seq;
		chomp $plus;
		chomp $qual;
		$line_count = $line_count + 3;
		if ( $plus eq '+' ) {
			$count++;
			print ">" . substr($desc, 1) . "\n$seq\n";
		}
		else {
			die "Incorrect FASTQ format around line " . ( $line_count - 1 ) . ".\n" ;	
		}
	}	
}	
close IN;

print STDERR "Found $count reads in input file.\n";
