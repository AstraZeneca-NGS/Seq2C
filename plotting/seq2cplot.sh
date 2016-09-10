#!/bin/bash

SAMPLE=`basename $1 .wg.cna.lr.txt`
perl -ne '{@a = split; $a[1] =~ s/:|-/\t/g; print join("\t", @a), "\n";}' $1 | joinLine -n -s -c 2 chr_loc_hg19.txt | perl -ne 'BEGIN { print join("\t", "Gene", "Chr", "Pos", "Log2Ratio"), "\n";} {$n++; @a = split; $a[2] += $a[9]; print join("\t", @a[0..1], $n, $a[8]), "\n" unless( $a[1] eq "chrY" );}' | plot.r3jpg 3 4 CNA.$SAMPLE 2
