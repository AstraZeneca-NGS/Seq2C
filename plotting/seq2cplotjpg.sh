#!/bin/bash

SAMPLE=$2
grep $SAMPLE $1 | perl -ne '{@a = split(/\t/, $_, -1); $a[2] = "chr$a[2]" unless( $a[2] =~ /chr/ ); print join("\t", @a);}' | joinLine -n -s -c 3 chr_loc_hg19.txt | perl -ne '{@a = split; $a[3] += $a[11]; print join("\t", @a[1..3], $a[10]), "\n" unless( $a[2] eq "chrY" );}' | sort -k3n | perl -ne 'BEGIN { print join("\t", "Gene", "Chr", "Pos", "Log2Ratio"), "\n";} { @a = split(/\t/); $a[2] = ++$n; print join("\t", @a);}' | plot.r3jpg 3 4 CNA.$SAMPLE 2
