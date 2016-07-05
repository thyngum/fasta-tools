#!/usr/bin/perl

# Cluster identical sequences in a multi-FASTA file.

#   fasta-consolidate.pl -in <filename> [-out <output>] 

# 	-in 	Input multi-FASTA file
#	-out 	Ouput filename (outputs to SDTOUT by default)

# WARNING: Probably not suitable for very large sequences or datasets!

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output );

die "Usage: fasta-consolidate.pl -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my $seqio_in = Bio::SeqIO->new(-file => $filename, 
                                   -format => 'fasta' ) or die $!;

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
                                 -format => 'fasta');	
}

my %seqs = ();
my %labels = ();
my $count = 0;

# Hash sequences from input file by using the sequence itself as a key
while ( my $seq = $seqio_in->next_seq ) {
	$count++;
	$seqs{uc $seq->seq()}++;
	$labels{uc $seq->seq()} = $labels{$seq->seq()} . $seq->display_id() . " ";
}
# Count the number of elements after hashing
my $count_consolidated = keys %seqs;

# Write output
if ( $count_consolidated != $count ) {
	open (OUTPUT, ">$output") if ( $output );

	my $cluster = 0;
	my $table = "Cluster\tNum\tLabels\n";
	foreach my $key ( keys %seqs ) {
		$cluster++;
		if ( $seqs{$key} > 1 ) {
			if ( $output ) {			
				print OUTPUT ">" . $labels{$key} . "\n" . sblock($key);
			}
			else {
				print ">" . $labels{$key} . "\n" . sblock($key);
			}	
			$tmp = $labels{$key};		
			$tmp =~ s/\s/\t/g;
			$table .= $cluster . "\t" .  $seqs{$key} . "\t" .  $tmp . "\n";
		} 
		else {
			if ( $output ) {
				print OUTPUT ">" . $labels{$key} . "\n" . sblock($key);
			}
			else {
				print ">" . $labels{$key} . "\n" . sblock($key);
			}	
			$table .= $cluster . "\t" . $seqs{$key} . "\t" . $labels{$key} . "\n";
		}
	}

	close OUTPUT if ( $output );

	print STDERR $table;
}
else {
	print STDERR "Found $count unique sequences in input file, no need to consolidate!\n";
}
