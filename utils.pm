#!/usr/bin/perl

sub complement {

	# Returns the complement of a DNA sequence (chars other than A, T, G or C are 
	# replaced by N). Output sequence is in uppercase.
	# 	Usage: complement(sequence)

	my $seq = uc $_[0];

	$rev = "";
	for ( my $i = 0; $i < length($seq); $i++ ) {
		$char = substr $seq, $i, 1;
		if ( $char eq 'A' ) {
			$rev .= 'T';
		}
		elsif ( $char eq 'C' ) {
			$rev .= 'G';
		}
		elsif ( $char eq 'G' ) {
			$rev .= 'C';
		}
		elsif ( $char eq 'T' ) {
			$rev .= 'A';
		}						
		else {
			$rev .= 'N';
		} 
	}

	return $rev;
}


sub sblock {

	# Returns a sequence block with n chars per line (n = 60 by default)
	# 	Usage: sblock(sequence, n)

	my ( $seq, $n ) = @_;
	$n = 60 unless ( $n );

	$block = "";
	while ( my $chunk = substr($seq, 0, $n, "") ) {
		$block .= "$chunk\n";
	}

	return $block;
}

1;
