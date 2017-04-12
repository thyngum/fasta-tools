#!/usr/bin/perl

# Compares two multi-FASTA files taking into account the number of sequences and 
# the sequences per se. 

#   fasta-diff.pl [-f <format>] -sub <filename> -qry <filename>

#      -f     Input format: fasta, fastq, genbank, etc. (guessed if omitted).
#      -d     Detailed output (if sequences differ, includes them in output).    

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'sub=s' => \$sub, 'qry=s' => \$qry, 'f=s' => \$format, 'd' => \$detailed );

die "Usage: fasta-diff.pl [-f <format>] -sub <filename> -qry <filename>\n" if ( ! $sub or ! $qry );

my $sub_file = File::Spec->rel2abs($sub);
die "File \'$sub_file\' not found.\n" if ( ! -e $sub_file );

my $qry_file = File::Spec->rel2abs($qry);
die "File \'$qry_file\' not found.\n" if ( ! -e $qry_file );

if ( ! $format ) {
	$seqio_sub = new Bio::SeqIO(-file   => $sub_file) or die $!;
	$seqio_qry = new Bio::SeqIO(-file   => $qry_file) or die $!;
}
else {
	$seqio_sub = new Bio::SeqIO(-format => $format,
	                           -file   => $sub_file) or die $!;
	$seqio_qry = new Bio::SeqIO(-format => $format,
	                           -file   => $qry_file) or die $!;	
}

# Process subject file
my $len;
my $count1 = 0;
my $total1 = 0;
my $max1;
my $min1;
my %seqs1 = ();
while ( my $seq = $seqio_sub->next_seq ) {
	$count1++;

	$len = $seq->length();
	$total1 += $len;
	$min1 = $len if ( !$min1 || $len < $min1 );
	$max1 = $len if ( !$max1 || $len > $max1 );

	$seqs1{$seq->display_id} = uc $seq->seq;
}
my $average1 = $total1 / $count1;

# Process query file
my $count2 = 0;
my $max2;
my $min2;
my %seqs2 = ();
while ( $seq = $seqio_qry->next_seq ) {
	$count2++;

	$len = $seq->length();
	$total2 += $len;
	$min2 = $len if ( !$min2 || $len < $min2 );
	$max2 = $len if ( !$max2 || $len > $max2 );	

	$seqs2{$seq->display_id} = uc $seq->seq;
}
my $average2 = $total2 / $count2;

# Compute global differences
my $count_diff = $count1 - $count2;
my $max_diff = $max1 - $max2;
my $min_diff = $min1 - $min2;
my $total_diff = $total1 - $total2;
my $average_diff = $average1 - $average2;

# Find differences in individual sequences
my $diff_str = "";
foreach my $key ( sort keys %seqs1 ) {
	if ( $seqs2{$key} ) {

		my $seq_diff = "";
		my $length_diff = length($seqs1{$key}) - length($seqs2{$key}); 

		if ( $length_diff ) {
			$seq_diff = $length_diff;
		}
		else {			
			$seq_diff = "(same length but sequences differ)" if ( $seqs1{$key} ne $seqs2{$key} );
		}
		 
		if ( $seq_diff ) {
			$diff_str = $diff_str . "$key\t$key\t$seq_diff\n";

			if ( $detailed ) {
				$diff_str = $diff_str . sblock($seqs1{$key}) . "---\n" . sblock($seqs2{$key});
			}
		}
	}
	else {
		$diff_str = $diff_str . "$key\t-\t-\n";
	}
}
foreach $key ( sort keys %seqs2 ) {
	unless ( $seqs1{$key} ) {
		$diff_str = $diff_str . "-\t$key\t-\n";
	}
}

# Print differences
if ( $count_diff != 0 ) {
	print STDERR "Count\t$count1\t$count2\t$count_diff\n";
}
if ( $max_diff != 0 ) {
	print STDERR "Maximum\t$max1\t$max2\t$max_diff\n";
}
if ( $min_diff != 0 ) {
	print STDERR "Minimum\t$min1\t$min2\t$min_diff\n";
}
if ( $average_diff != 0 ) {
	print STDERR "Average\t" . sprintf("%.0f", $average1) . "\t" . sprintf("%.0f", $average2) . "\t" . sprintf("%.4f", $average_diff) . "\n";
}
if ( $total_diff != 0 ) {
	print STDERR "Total\t$total1\t$total2\t$total_diff\n";
}

if ( $diff_str ne "" ) {
	print $diff_str;
}