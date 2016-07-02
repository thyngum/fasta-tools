#!/usr/bin/perl

# Concatenates all the sequences in a multi-FASTA file. 

# 	fasta-union.pl [-aa] [-gap <n>] [-t] -in <filename> [-out <output>]

# 		-in 	Input multi-FASTA file
#		-out 	Ouput filename (outputs to SDTOUT by default)
#		-aa	Input sequences are amino acid sequences (default is nucleotide)
# 		-gap	Add a gap of <n> Ns or Xs between each individual sequence.
# 		-t	Write a table file (.tab) in EMBL format with the coordinates of
# 	  		of the original sequences (only for nucleotide sequences).

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

GetOptions ( 'in=s' => \$input, 'out=s' => \$output, 't' => \$tab, 'aa' => \$aa, 'gap=i' => \$gap );

die "Usage: fasta-union.pl [-aa] [-gap <n>] [-t] -in <filename> [-out <output>]\n" if ( ! $input );

die "Flags -aa and -t cannot be used together!\n" if ( $tab && $aa );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );
(my $name, my $path, my $suffix) = fileparse($filename, qr/\.[^.]*/);

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
if ( $gap ) {
	if ( $aa ) {
		$stretch = "X" x int($gap);
	}
	else {
		$stretch = "N" x int($gap);
	}	
}                         

my $count = 0;
my $base_count = 0;
my $union = "";
my $first = 1;
my $tab_str = "";
while ( my $seq = $seqio_object->next_seq ) {
	if ( $gap ) {
		if ( $first ) {
			$union = $union . $seq->seq();
			$first = 0;
		}
		else {
			$union = $union . $stretch . $seq->seq();
			$base_count = $base_count + $gap;
		}
	}
	else { 
		$union = $union . $seq->seq();		
	}
	my $length = $seq->length;	

	if ( $tab ) {
		$start = $base_count + 1;
		$end = $base_count + $length;	
		$tab_str = $tab_str . "FT   misc_feature    $start..$end\n";
		$tab_str = $tab_str . "FT                   /note=\"" . $seq->display_id() . "\"\n";	
		$base_count = $base_count + $length;	
	}

	$count++;
}	

if ( $output ) {
	$seq_union = Bio::Seq->new(-display_id => $name, 
	                           -seq => $union ) or die $!;
	$seqio_out->write_seq($seq_union);	
}
else {
	print ">$name\n";
	print $union . "\n";	
}

if ( $tab ) {
	open TAB, ">$tab_file";

	# Print the EMBL header here
	print TAB "ID   $name; SV 1; linear; unassigned DNA; STD; UNC; $base_count BP.\n";
	print TAB "XX\n";
	print TAB "AC   $name\n";
	print TAB "XX\n";
	print TAB "XX\n";
	print TAB "XX\n";
	print TAB "FH   Key             Location/Qualifiers\n";
	print TAB $tab_str;
	close TAB
}
