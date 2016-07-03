#!/usr/bin/perl

# Compare two FASTA files according to the number of sequences and their basecounts.
# Prints nothing if no difference is found.

# 	fasta-diff.pl -sub <filename> -qry <filename>

#   	-qry	Query file
#   	-sub	Reference sequence file.

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'qry=s' => \$query, 'sub=s' => \$ref );

die "Usage: fasta-diff.pl -sub <filename> -qry <filename>.\n" if ( ! $query || ! $ref );

$query_file = File::Spec->rel2abs($query);
die "Query file not found.\n" if ( ! -e $query_file );

$ref_file = File::Spec->rel2abs($ref);
die "Reference file not found.\n" if ( ! -e $ref_file );

# Process reference sequence
$seqio_object = Bio::SeqIO->new(-file => $ref_file, 
								-format => $format ) or die $!;

$ref_count = 0;
$ref_size = 0;
$ref_A_count = 0;
$ref_T_count = 0;
$ref_C_count = 0;
$ref_G_count = 0;
$ref_N_count = 0;
$ref_a_count = 0;
$ref_t_count = 0;
$ref_c_count = 0;
$ref_g_count = 0;
$ref_n_count = 0;
$ref_gap_count = 0;

while ( $seq = $seqio_object->next_seq ) {
	$ref_count++;
	$len = $seq->length;
	$ref_size = $ref_size + $len;
	$ref_min = $len if ( !$ref_min || $len < $ref_min );
	$ref_max = $len if ( !$ref_max || $len > $ref_max );

	$ref_a_count += $seq->seq =~ tr/a//;
	$ref_t_count += $seq->seq =~ tr/t//;
	$ref_c_count += $seq->seq =~ tr/c//;
	$ref_g_count += $seq->seq =~ tr/g//;
	$ref_n_count += $seq->seq =~ tr/n//;

	$UCseq = uc($seq->seq);
	$ref_A_count += $UCseq =~ tr/A//;
	$ref_T_count += $UCseq =~ tr/T//;
	$ref_C_count += $UCseq =~ tr/C//;
	$ref_G_count += $UCseq =~ tr/G//;
	$ref_N_count += $UCseq =~ tr/N//;
	$ref_gap_count += () = $UCseq =~ /NN+/g;
}
$ref_average = sprintf("%.0f", $ref_size/$ref_count);
$ref_other = ($ref_size - ($ref_A_count + $ref_T_count + $ref_C_count + $ref_G_count + $ref_N_count));
$ref_lowq = $ref_a_count + $ref_t_count + $ref_c_count + $ref_g_count + $ref_n_count;
##

# Process query sequence
$seqio_object = Bio::SeqIO->new(-file => $query_file, 
								-format => $format ) or die $!;

$qry_count = 0;
$qry_size = 0;
$qry_A_count = 0;
$qry_T_count = 0;
$qry_C_count = 0;
$qry_G_count = 0;
$qry_N_count = 0;
$qry_a_count = 0;
$qry_t_count = 0;
$qry_c_count = 0;
$qry_g_count = 0;
$qry_n_count = 0;
$qry_gap_count = 0;

while ( $seq = $seqio_object->next_seq ) {
	$qry_count++;
	$len = $seq->length;
	$qry_size = $qry_size + $len;
	$qry_min = $len if ( !$qry_min || $len < $qry_min );
	$qry_max = $len if ( !$qry_max || $len > $qry_max );

	$qry_a_count += $seq->seq =~ tr/a//;
	$qry_t_count += $seq->seq =~ tr/t//;
	$qry_c_count += $seq->seq =~ tr/c//;
	$qry_g_count += $seq->seq =~ tr/g//;
	$qry_n_count += $seq->seq =~ tr/n//;

	$UCseq = uc($seq->seq);
	$qry_A_count += $UCseq =~ tr/A//;
	$qry_T_count += $UCseq =~ tr/T//;
	$qry_C_count += $UCseq =~ tr/C//;
	$qry_G_count += $UCseq =~ tr/G//;
	$qry_N_count += $UCseq =~ tr/N//;
	$qry_gap_count += () = $UCseq =~ /NN+/g;
}
$qry_average = sprintf("%.0f", $qry_size/$qry_count);
$qry_other = ($qry_size - ($qry_A_count + $qry_T_count + $qry_C_count + $qry_G_count + $qry_N_count));
$qry_lowq = $qry_a_count + $qry_t_count + $qry_c_count + $qry_g_count + $qry_n_count;
##
	
print STDERR "Count\t$qry_count\t$ref_count\t" . ($qry_count - $ref_count) .  "\n" if ( $qry_count != $ref_count );
print STDERR "Maximum\t$qry_max\t$ref_max\t" . ($qry_max - $ref_max) .  "\n" if ( $qry_max != $ref_max );
print STDERR "Minimum\t$qry_min\t$ref_min\t" . ($qry_min - $ref_min) .  "\n" if ( $qry_min != $ref_min );
print STDERR "Average\t$qry_average\t$ref_average\t" . ($qry_average - $ref_average) .  "\n" if ( $qry_average != $ref_average );
print STDERR "Total\t$qry_size\t$ref_size\t" . ($qry_size - $ref_size) .  "\n" if ( $qry_size != $ref_size );
print STDERR "  A\t$qry_A_count\t$ref_A_count\t" . ($qry_A_count - $ref_A_count) . "\n" if ( $qry_A_count != $ref_A_count );
print STDERR "  T\t$qry_T_count\t$ref_T_count\t" . ($qry_T_count - $ref_T_count) . "\n" if ( $qry_T_count != $ref_T_count );
print STDERR "  G\t$qry_G_count\t$ref_G_count\t" . ($qry_G_count - $ref_G_count) . "\n" if ( $qry_G_count != $ref_G_count );
print STDERR "  C\t$qry_C_count\t$ref_C_count\t" . ($qry_C_count - $ref_C_count) . "\n" if ( $qry_C_count != $ref_C_count );
print STDERR "  N (Gaps)\t$qry_N_count ($qry_gap_count)\t$ref_N_count ($ref_gap_count)\t" . ($qry_N_count - $ref_N_count) .  " (" . ($qry_gap_count - $ref_gap_count) . ")\n" if ( $qry_N_count != $ref_N_count );
print STDERR "  Other\t$qry_other\t$ref_other\t" . ($qry_other - $ref_other) . "\n" if ( $qry_other != $ref_other );
print STDERR "  L/Q\t$qry_lowq\t$ref_lowq\t" . ($qry_lowq - $ref_lowq) . "\n" if ( $qry_lowq != $ref_lowq );


