#!/usr/bin/env perl

use warnings;
use strict;

while( <> ) {
    my @a = split(/\t/);
    next unless( /Del/ || /Amp/ );
    my $desc = /Whole/ ? "$a[11] of $a[11]" : "$a[10] of $a[10]";
    print join("\t", $a[0], "", "copy-number-alteration", $a[1], "NA", "-", "-", "-", "-", "-", $a[12], $desc, $a[12], /Del/ ? "loss" : "amplification", "-", "-", "-", "-", "-", "-", "-"), "\n";
}
