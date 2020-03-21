#!/usr/bin/perl

# Filters sequences from a FASTA/Q file according to specific set of criteria. Save  
# the filtered sequences to the file specified via the -out option, or append an  
# underscore (_) to the original filename if no output is specified.

#  fasta-filter.pl [-f <format>] [-n <size>] [-gc <range>] -in <filename> [-out <output>]

# 	-f <format>		Format can be 'fasta' or 'fastq' (default is fasta)
# 	-n <size>  		Filter by size, extract all the sequences with a length
#             		greater than <size>
# 	[-gc <range>]	Filter by gc-content, extract all the sequence with a gc percent 
#	             	in the given range, e.g. 20..40

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'in=s' => \$input, 'f=s' => \$format, "n=i" => \$length, 'gc=s' => \$gc_range, 'out=s' => \$output );

die "Usage: fasta-filter [-f <format>] [-n <size>] [-gc <range>] -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my $length = 0 if ( ! $length );
my $min_gc = 0;
my $max_gc = 100;

if ( $gc_range ) {
	my @split = split(/\.\./, $gc_range);
	$min_gc = $split[0];
	$max_gc = $split[1];
	
	die "Incorrect gc range!\n" if ( $min_gc > $max_gc or $min_gc < 0 or $max_gc < 0 or $min_gc > 100 or $max_gc > 100 );
}


if ( $length == 0 and ( $min_gc == 0 and $max_gc == 100 ) ) {
	die "Nothing to do!\n"; 
}
else {
	print "Extracting sequences larger than $length bp with a GC-content between $min_gc and $max_gc% ...\n";
}

if ( $format ) {
	if ( ( $format ne 'fasta' ) and ( $format ne 'fastq' ) ) {
		die "Can't work with the \'$format\' format.\n";
	}
}
else {
	$format = "fasta";
}

($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
if ( ! $output ) {
	$output = $path . $name . "_" . $suffix;			
}
else {		
	$output = File::Spec->rel2abs($output);
}

$count = 0;
$count_ok = 0;

# Process a FASTA file
if ( $format eq "fasta" ) {
	$seqio_object = Bio::SeqIO->new(-file => $filename, 
                                    -format => 'fasta' ) or die $!;
	$seqio_output = Bio::SeqIO->new(-file => ">$output" ,
                                    -format => 'fasta');

	while ( $seq = $seqio_object->next_seq ) {
		$count++;
	
		if ( $seq->length() > $length ) {
			$gc_content = $seq->seq =~ tr/gcGC//;
			$gc_percent = $gc_content / $seq->length() * 100; 

			if ( $gc_percent >= $min_gc && $gc_percent <= $max_gc ) {		
				$count_ok++;
				# print STDERR "$count_ok: " . $seq->display_id() . " (" . $seq->length() . ")\n";		
				$seqio_output->write_seq($seq);
			}		
		}
	}
}

# Process a FASTQ file
else {
	open FILE, $filename or die("Could not open input file.\n");
	open OUTPUT, ">$output" or die("Could not create output file.\n");

	$line_count = 0;
	while (<FILE>)  {  
		$line = $_;
		$line_count++;
		$start = substr $line, 0, 1;
		$dummy = substr $line, 0, 2;		
		if ( ( $start eq '@' ) && ( $dummy ne '@@' ) ) {
			$desc = $line;
			$seq = <FILE>;
			$plus = <FILE>;
			$qual = <FILE>;
			chomp $desc;
			chomp $seq;
			chomp $plus;
			chomp $qual;
			$line_count = $line_count + 3;
			if ( $plus eq '+' ) {
				$count++;
			}
			else {
				die "Incorrect FASTQ format around line " . ($line_count - 1) . ".\n" ;	
			}
					
			if ( length($seq) >= $length ) {
			
				$gc_content = $seq =~ tr/gcGC//;
				$gc_percent = $gc_content / length($seq) * 100; 
			
				if ( $gc_percent >= $min_gc && $gc_percent <= $max_gc ) {	
					$count_ok++;
					print OUTPUT "$desc\n$seq\n$plus\n$qual\n";
				}
			}			
				
		}	
	}	
	close FILE;
	close OUTPUT;
	
}

print STDERR "\n$count_ok of $count sequences written to $output\n";

