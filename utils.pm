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

sub basecounts {

	# Counts bases from a sequence and returns the counts as a list (A, C, G, T, other) or
	# (A, C, G, T, a, c, g, t, other)
	# 	Usage: basecounts(sequence, case-sensitive)

	my $seq = $_[0];
	my $case_sensitive = $_[1];

	my $len = length($seq);

	if ( $case_sensitive ) {
		my $A = $seq =~ tr/A//;
		my $C = $seq =~ tr/C//;
		my $G = $seq =~ tr/G//;
		my $T = $seq =~ tr/T//;
		my $a = $seq =~ tr/a//;
		my $c = $seq =~ tr/c//;
		my $g = $seq =~ tr/g//;
		my $t = $seq =~ tr/t//;
		my $sum = $A + $C + $G + $T + $a + $c + $g + $t;
		my $other = $len - $sum;

		return($A, $C, $G, $T, $a, $c, $g, $t, $other);
	}
	else {
		$seq = uc $seq;
		my $A = $seq =~ tr/A//;
		my $C = $seq =~ tr/C//;
		my $G = $seq =~ tr/G//;
		my $T = $seq =~ tr/T//;
		my $sum = $A + $C + $G + $T;
		my $other = $len - $sum;

		return($A, $C, $G, $T, $other);
	}
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
