#!/usr/bin/perl

# Fix a pair of FASTQ files containing paired-end Illumina reads. Read the input prefix_1 and prefix_2
# FASTQ files and throw out orphaned reads. Save the corrected files with the prefix passed via the
# -out option, or append "_fixed" to the original prefix if no output is specified 

#   fastq-filter.pl -f <format> -size <n> -fq1 <filename> -fq2 <filename> [-out <prefix>]

# 	-f 	Format can be 'fasta' or 'fastq' (default is fasta)
# 	-size  	Filter by size, extract reads with a length equal or larger than <n>
# 	-gc 	Filter by GC-content, extract reads with a GC percent in the given range, e.g. 20..40

# WARNING! This script is NOT optimized for memory use! 

# It is designed to work 
# only with reads with multiplexed protocols (named xxx 1:0:... and xxx 2:0:...)

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output );

die "fastq-filter.pl -in <prefix> [-out <prefix>]\n" if ( ! $input );



$fastq1 = $input . "_1.fastq";
$fastq2 = $input . "_2.fastq";

$file_fastq1 = File::Spec->rel2abs($fastq1);
die "File \'$file_fastq1\' not found.\n" if ( ! -e $file_fastq1 );

$file_fastq2 = File::Spec->rel2abs($fastq2);
die "File \'$file_fastq2\' not found.\n" if ( ! -e $file_fastq2 );

if ( $output ) {
	$out_fastq1 = $output . "_1.fastq";
	$out_fastq2 = $output . "_2.fastq";
}
else {
	$out_fastq1 = $input . "_fixed_1.fastq";
	$out_fastq2 = $input . "_fixed_2.fastq";	
}

# Process _1.fastq
open FASTQ1, $file_fastq1 or die("Could not open input file.\n");

%descs1 = ();
%seqs1 = ();
%quals1 = ();
$count1 = 0;

$line_count = 0;
while ( $line = <FASTQ1> )  {  
	$line_count++;

	$start = substr $line, 0, 1;
	$dummy = substr $line, 0, 2;		
	if ( ( $start eq '@' ) && ( $dummy ne '@@' ) ) {
		$desc = $line;
		$seq = <FASTQ1>;
		$plus = <FASTQ1>;
		$qual = <FASTQ1>;
		chomp $desc;
		chomp $seq;
		chomp $plus;
		chomp $qual;
		$line_count = $line_count + 3;
		if ( $plus eq '+' ) {
			$count1++;
		}
		else {
			die "Incorrect format in $fastq1 around line " . ($line_count - 1) . ".\n" ;	
		}
				
		@split = split(" ", $desc);
		$descs1{$split[0]} = $desc; 
		$seqs1{$split[0]} = $seq;
		$quals1{$split[0]} = $qual;	
	}	
}
close FASTQ1;
print "$fastq1 has $count1 reads\n"; #

# Process _2.fastq
open FASTQ2, $file_fastq2 or die("Could not open input file.\n");

open FASTQ1_OUTPUT, ">$out_fastq1" or die("Could not create output file.\n");
open FASTQ2_OUTPUT, ">$out_fastq2" or die("Could not create output file.\n");

$count2 = 0;
$count_ok = 0;

$line_count = 0;
while ( $line = <FASTQ2> )  {  
	$line_count++;

	$start = substr $line, 0, 1;
	$dummy = substr $line, 0, 2;		
	if ( ( $start eq '@' ) && ( $dummy ne '@@' ) ) {
		$desc = $line;
		$seq = <FASTQ2>;
		$plus = <FASTQ2>;
		$qual = <FASTQ2>;
		chomp $desc;
		chomp $seq;
		chomp $plus;
		chomp $qual;
		$line_count = $line_count + 3;
		if ( $plus eq '+' ) {
			$count2++;
		}
		else {
			die "Incorrect format in $fastq2 around line " . ($line_count - 1) . ".\n" ;	
		}
				
		@split = split(" ", $desc);
		$current = $split[0];
		if ( $seqs1{$current} ) {
			$count_ok++;
			print FASTQ1_OUTPUT $descs1{$current} . "\n" . $seqs1{$current} . "\n$plus\n" . $quals1{$current} . "\n";
			print FASTQ2_OUTPUT "$desc\n$seq\n$plus\n$qual\n";
		}
	}	
}
close FASTQ2;
print "$fastq2 has $count2 reads\n"; #

close FASTQ1_OUTPUT;
close FASTQ2_OUTPUT;

print "\n$count_ok reads written to output files.\n";

