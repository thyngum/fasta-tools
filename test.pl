#!/usr/bin/perl

use utils;

$str = "ATGGggNRTaGG";

($A, $C, $G, $T, $other) = basecounts($str, 1);

print "A: $A\nC: $C\nG: $G\nT: $T\nOther: $other\n";
