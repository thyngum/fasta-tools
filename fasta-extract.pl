#!/usr/bin/perl

# Reads a text file containing a list of sequence and coordinates, then extracts the
# subsequences specified by those coordinates from a multi-FASTA file and writes them
# in FASTA format.

#  fasta-extract.pl -l <list> -in <filename> [-out <output>]

#    -l <list>	Tab-delimited file with a subsequence specification per line, 
#             	with the format:
#
#            	ID	SOURCE	START	END
#
#            	Where SOURCE is a sequence in input file, and ID, START and END are
#            	the name and coordinates of the subsequence to be extracted from SOURCE.
#            	If start > end the subsequence is extracted from the reverse strand 
#            	of SOURCE (but using the positions of the forward strand). If no
#            	START and END is given, the whole source sequence is extracted.
#    -in    	Input multi-FASTA file
#    -out   	Ouput filename (outputs to SDTOUT by default).

# WARNING: This script loads the whole input file in memory. Be careful with
# large sequence files (e.g. larger than the available RAM)!!!

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use File::Spec;
use Getopt::Long;

# Import local package utils.pm
use FindBin;
use lib $FindBin::Bin;
use utils;

GetOptions ( 'in=s' => \$input, 'l=s' => \$list, 'out=s' => \$output );

die "Usage: fasta-extract.pl -l <list> -in <filename> [-out <output>].\n" if ( (! $input) || (! $list) );

$filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

$listfile = File::Spec->rel2abs($list);
die "File \'$listfile\' not found.\n" if ( ! -e $listfile );

if ( $output ) {
	$output_filename = File::Spec->rel2abs($output);
	$seqio_out = Bio::SeqIO->new(-file => ">$output_filename", 
	                             -format => 'fasta') or die $!;	
}

# Load input file in memory and open output file
%seqs = ();

$seqio_object = Bio::SeqIO->new(-file => $filename, 
                                -format => $format ) or die $!;								    

while ( $seq = $seqio_object->next_seq ) {
	$seqs{$seq->display_id()} = $seq->seq;
}							    

open LIST, $listfile;

while ( $line = <LIST> ) {
	chomp $line;
	if ( $line ne "" ) {		

		( $id, $source, $start, $end ) = split("\t", $line);

		unless ( $seqs{$source} ) {
			print STDERR "Warning: No sequece with ID \'$source\' in input file!\n";
		}
		else {
			$seq_obj = Bio::Seq->new(-seq => $seqs{$source} );
			if ( ! $start && ! $end ) {
				$seq_obj-> display_id($id);
				print ">$id\n"; ### TODO: Check if an ID is already in ouput, to avoid duplicate IDs!
				print sblock($seq_obj->seq());
				$seqio_output->write_seq($seq_obj) if ( $output );						
			}
			else {
				if ( $start <= $end ) {
					if ( $start > 0 and $end <= length($seq_obj->seq()) ) {
						$subseq = $seq_obj->subseq($start, $end);
						$subseq_obj = Bio::Seq->new(-seq => $subseq,
						                            -display_id => "$id",
						                            -desc => "    $source:$start..$end" );	
						print ">" . $subseq_obj->display_id() . $subseq_obj->desc() . "\n";
						print sblock($subseq_obj->seq());									
						$seqio_output->write_seq($subseq_obj) if ( $output );
					}
					else {
						print STDERR "Warning: sequence coordinates out of range for $id!\n";
					}					
				}
				else {
		
					# This is the case of a subsequence in the reverse strand.
					# Start and end are inverted, the subsequence is extracted, 
					# and then the reverse complement is written to the file.
					# Start and end coordinates are converted to the reverse strand by:
					#  rev = length of the original sequence - start/end + 1

					# TODO: Check validity of coordinates for subsequences!

					$subseq = $seq_obj->subseq($end, $start); 
					$tmp_obj = Bio::Seq->new(-seq => $subseq );
					$subseq_obj = $tmp_obj->revcom;
					$subseq_obj->display_id($id);
					$c_start = $seq_obj->length - $start + 1;
					$c_end = $seq_obj->length - $end + 1;
					$subseq_obj->desc("    $source:complement($c_start..$c_end)");
					print sblock($subseq_obj->seq());
					$seqio_output->write_seq($subseq_obj) if ( $output );					
				}
			}
		}
	}
}	
