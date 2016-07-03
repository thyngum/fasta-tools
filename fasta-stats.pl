#!/usr/bin/perl

# Prints basic stats of a multi-FASTA file. It can also read any other sequence file
# format supported by BioPerl. 

# 	fasta-stats.pl [-f <format>] [-d] [-p] -in <filename>

#   	-f   	Input format: fasta, fastq, genbank, etc. (guessed if omitted)
# 		-c      Clean output (minimum descriptions, suitable for redirection)
#   	-d		Detailed stats (including base counts)
#   	-p		Also prints a tab-separated list with description and length of 
# 				individual sequences to STDOUT.

use Bio::Seq;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use utils;

GetOptions ( 'in=s' => \$input, 'f=s' => \$format, 'd' => \$detailed, 'p' => \$print );

die "Usage: fasta-stats.pl [-f <format>] [-d] [-p] -in <filename>\n" if ( ! $input );

$filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

if ( ! $format ) {
	$seqio_in = new Bio::SeqIO(-file   => $filename) or die $!;
}
else {
	die "Use fastq-stats.pl for FASTQ files!\n" if ( $format eq 'fastq' );

	$seqio_in = new Bio::SeqIO(-format => $format,
	                           -file   => $filename) or die $!;
}

my $type; 
my $len;
my $min;
my $max;
my $average;
my $count = 0;
my $total_size = 0;
my ( $A_count, $C_count, $G_count, $T_count ); 
my ( $a_count, $c_count, $g_count, $t_count );
my $other_count;
my $other_percent;
my $lc;
my $lc_percent;
my $gap_count;

while ( my $seq = $seqio_in->next_seq ) {
	$count++;

	if ( $count == 1 ) {
		$type = $seq->alphabet();
		die "RNA sequences not supported, please convert to DNA!\n" if ( $type eq 'rna' );
		print STDERR "Alphabet is 'protein' ...\n\n" if ( $type eq 'protein');
	}

	$len = $seq->length();
	$total_size += $len;
	$min = $len if ( !$min || $len < $min );
	$max = $len if ( !$max || $len > $max );

	if ( $type eq 'dna' ) {
		my ($A, $C, $G, $T, $a, $c, $g, $t, $other) = basecounts($seq->seq(), 1);

		$A_count += $A; $C_count += $C; $G_count += $G; $T_count += $T;
		$a_count += $a; $c_count += $c; $g_count += $g; $t_count += $t;
		$other_count += $other;

		$uc_seq = uc $seq->seq();
		$gap_count += () = $uc_seq =~ /N+/g;
	}

	if ( $print ) {		
		print $seq->display_id() . "\t" . $seq->length() . "\t";
		if ( $type eq 'dna' ) {
			my $current_gaps = () = $uc_seq =~ /N+/g;
			print $current_gaps . "\t";
		}
		print $seq->desc() . "\n";
	}
}
print "\n" if ( $print );

$average = $total_size / $count;
if ( $type eq 'dna' ) {
	$lc = $a_count + $t_count + $c_count + $g_count;
	$lc_percent = $lc / $total_size * 100;
	$A_count += $a_count; $C_count += $c_count; $G_count += $g_count; $T_count += $t_count;
	$other_percent = $other_count / $total_size * 100;
	$gc_content = ( $G_count + $C_count ) / ( $total_size - $other_count ) * 100;
}

unless ( $detailed ) {
	print STDERR "Count\t",  
	  "Total\t",
	  "Min\t",
	  "Average\t",
	  "Largest";

	if ( $type eq 'dna' ) {
		print STDERR "\tN\t",
		  "N%\t",
		  "Gaps\t",
		  "GC%";
	}
	print "\n";	

	print "$count\t",
	  "$total_size\t",
	  "$min\t",
	  sprintf("%.0f", $average), "\t",
	  "$max";

	if ( $type eq 'dna' ) {
		print "\t$other_count\t",
		  sprintf("%.2f", $other_percent), "\t",
		  "$gap_count\t",
		  sprintf("%.2f", $gc_content), "\t";
	}
	print "\n";	
}
else {
	print STDERR "Count\t$count\n",
	  "Maximum\t$max\n",
	  "Minimum\t$min\n",
	  "Average\t", sprintf("%.0f", $average), "\n",
	  "Total\t$total_size\n";

	if ( $type eq 'dna' ) {
		print "GC%\t", sprintf("%.2f", $gc_content), "\n",
		  "  A\t$A_count\t", sprintf("%.2f", ($A_count/$total_size)*100), "\n",
		  "  C\t$C_count\t", sprintf("%.2f", ($C_count/$total_size)*100), "\n",
		  "  G\t$G_count\t", sprintf("%.2f", ($G_count/$total_size)*100), "\n",
		  "  T\t$T_count\t", sprintf("%.2f", ($T_count/$total_size)*100), "\n",
		  "  Other\t$other_count\t", sprintf("%.2f", $other_percent), "\t($gap_count gaps)\n",
		  "  L/Q\t$lc\t", sprintf("%.2f", $lc_percent), "\n";
	}
}
