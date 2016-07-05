#!/usr/bin/perl

# Concatenates all the sequences in a multi-FASTA file. 

#   fasta-union.pl [-gap <n>] [-t] -in <filename> [-out <output>]

# 	-in 	Input multi-FASTA file
#	-out 	Ouput filename (outputs to SDTOUT by default)
# 	-gap	Add a gap of <n> Ns or Xs between each individual sequence.
# 	-t  	Write a table file (.tab) in EMBL format with the coordinates of
# 	    	of the original sequences (only for nucleotide sequences).

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output, 't' => \$tab, 'gap=i' => \$gap );

die "Usage: fasta-union.pl [-gap <n>] [-t] -in <filename> [-out <output>]\n" if ( ! $input );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );
my ($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);

my $seqio_object = Bio::SeqIO->new(-file => $filename, 
                                -format => 'fasta' ) or die $!;

if ( $output ) {
	my $output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
	                             -format => 'fasta');
}
if ( $tab ) {
	$tab_file = $path . $name . ".tab";
}
                    

my $count = 0;
my $union_seq = "";
my $union_length = 0;
my $tab_str = "";
my $type;
my $protein_gap = "X" x int($gap);
my $nucl_gap = "N" x int($gap);

while ( my $seq = $seqio_object->next_seq ) {
	$count++;

	if ( $count == 1 ) {
		$type = $seq->alphabet();
		die "The -t flag cannot be used with protein sequences!\n" if ( $type eq 'protein' and $tab );
	}
	
	if ( $gap ) {
		unless ( $count == 1 ) {
			if ( $type eq 'protein' ) {
				$union_seq = $union_seq . $protein_gap . $seq->seq();
			}
			else {
				$union_seq = $union_seq . $nucl_gap . $seq->seq();
			}			
			$union_length += $gap;
		}
		else {
			$union_seq .= $seq->seq();
		}
	}
	else {
		$union_seq .= $seq->seq();		
	}

	if ( $tab ) {
		$start = $union_length + 1;
		$end = $union_length + $seq->length();	
		$tab_str = $tab_str . "FT   misc_feature    $start..$end\n";
		$tab_str = $tab_str . "FT                   /note=\"" . $seq->display_id() . "\"\n";	
	}

	$union_length += $seq->length();
}	

if ( $output ) {
	$seq_union = Bio::Seq->new(-display_id => $name, 
	                           -seq => $union_seq ) or die $!;
	$seqio_out->write_seq($seq_union);	
}
else {
	print ">$name\n";
	print sblock($union_seq);	
}

if ( $tab ) {
	open TAB, ">$tab_file";
	print TAB "ID   $name; SV 1; linear; unassigned DNA; STD; UNC; $union_length BP.\n",
	  "XX\n",
	  "AC   $name\n",
	  "XX\n",
	  "XX\n",
	  "XX\n",
	  "FH   Key             Location/Qualifiers\n",
	  $tab_str;
	close TAB
}
