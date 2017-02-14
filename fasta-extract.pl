#!/usr/bin/perl

# Reads a text file containing a list of sequence and coordinates, then extracts the
# subsequences specified by those coordinates from a multi-FASTA file and writes them
# in FASTA format.

#  fasta-extract.pl -l <list> -in <filename> [-out <output>]

#    -l <list>	Tab-delimited file with a subsequence specification per line, 
#             	with the format:
#
#            	SOURCE  [ID  START  END]
#
#            	Where SOURCE is a sequence in input file and ID, START and END are
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

die "Usage: fasta-extract.pl -l <list> -in <filename> [-out <output>].\n" unless ( $input and $list );

my $filename = File::Spec->rel2abs($input);
die "File \'$filename\' not found.\n" if ( ! -e $filename );

my $listfile = File::Spec->rel2abs($list);
die "File \'$listfile\' not found.\n" if ( ! -e $listfile );

if ( $output ) {
	$seqio_out = Bio::SeqIO->new(-file => ">$output", 
	                             -format => 'fasta') or die $!;	
}

# Hash input sequences
my %seqs = ();
my $seqio_in = Bio::SeqIO->new(-file => $filename, 
                            -format => $format ) or die $!;								    
while ( my $seq = $seqio_in->next_seq ) {
	$seqs{$seq->display_id} = $seq->seq;
}							    

# Process list file
open LIST, $listfile;

while ( my $line = <LIST> ) {
	chomp $line;
	if ( $line ne "" ) {

		( my $source, my $id, my $start, my $end ) = split /\t/, $line;

		if ( $id or $start or $end ) {
			unless ( $id and $start and $end ) {
				die "Format of list file is: SOURCE  [ID  START  END]!\n";
			}	
		}

		if ( $seqs{$source} ) {			
			my $tmp_seq = Bio::Seq->new(-seq => $seqs{$source});			                            
			unless ( $start and $end ) {
				if ( $output ) {
					$tmp_seq->display_id($source);
					$seqio_out->write_seq($tmp_seq);
				}
				else {
					print ">$source\n"; ### TO-DO: Check if an ID is already in ouput, to avoid duplicate IDs!
					print sblock( $seqs{$source} );
				}
			}
			else {
				if ( $start <= $end ) {
					if ( $start > 0 and $end <= length($tmp_seq->seq()) ) {
						$subseq = Bio::Seq->new(-seq => $tmp_seq->subseq($start, $end),
						                        -display_id => $id,
						                        -desc => "\t$source:$start..$end" );	
						if ( $output ) {
							$seqio_out->write_seq($subseq);
						}
						else {
							print ">$id\t$source:$start..$end\n";
							print sblock($subseq->seq());								
						}
					}
					else {
						print STDERR "Warning: sequence coordinates for $id are out of range!\n";
					}					
				}
				else {
		
					# This is the case of a subsequence in the reverse strand.
					# Start and end are inverted, the subsequence is extracted, 
					# and then the reverse complement is written to the file.
					# Start and end coordinates are converted to the reverse strand by:
					#  rev = length of the original sequence - start/end + 1

					# TODO: Check validity of coordinates for subsequences!
					$subseq = Bio::Seq->new(-seq => $tmp_seq->subseq($end, $start)->revcom->seq,
					                        -display_id => $id,
					                        -desc => "\t$source:complement($c_start..$c_end)" );					
					my $c_start = $tmp_seq->length - $start + 1;
					my $c_end = $tmp_seq->length - $end + 1;

					if ( $output ) {
						$seqio_out->write_seq($subseq_obj) if ( $output );
					}
					else {
						print ">$id\t$source:complement($c_start..$c_end)\n";
						print sblock($subseq_obj->seq());
					}					
				}
			}
		}
		else {
			print STDERR "Warning: Could not find a sequence \'$source\' in input file!\n";
		}
	}
}	
