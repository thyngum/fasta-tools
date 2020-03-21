else {

	# Processing a FASTQ file
	open FILE, $filename or die "Could not open input file.\n";

	$line_count = 0;
	while (<FILE>)  {  
		$line = $_;
		$line_count++;
		$start = substr $line, 0, 1;
		$dummy = substr $line, 0, 2;		
		if ( ( $start eq '@' ) && ( $dummy ne '@@' ) ) {
			$desc = $line;
			$seq = <FILE>;
			$plus = <FILE>;
			$qual = <FILE>;
			chomp $desc;
			chomp $seq;
			chomp $plus;
			# chomp $qual;
			$line_count = $line_count + 3;
			if ( $plus eq '+' ) {
				$count++;
				$len = length($seq);
				$size = $size + $len;
				$min = $len if ( !$min || $len < $min );
				$max = $len if ( !$max || $len > $max );
				
				$a_count += $seq =~ tr/a//;
				$t_count += $seq =~ tr/t//;
				$c_count += $seq =~ tr/c//;
				$g_count += $seq =~ tr/g//;
				$n_count += $seq =~ tr/n//;

				$UCseq = uc($seq);
				$A_count += $UCseq =~ tr/A//;
				$T_count += $UCseq =~ tr/T//;
				$C_count += $UCseq =~ tr/C//;
				$G_count += $UCseq =~ tr/G//;
				$N_count += $UCseq =~ tr/N//;
				$gap_count += () = $UCseq =~ /NN+/g;
				
				print substr($desc, 1) . "\t" . $len . "\n" if ( $print );
				
			}
			else {
				die "Incorrect FASTQ format around line " . ($line_count - 1) . ".\n" ;	
			}	
		}	
	}	
	close FILE;
}

