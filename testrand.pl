#!/usr/bin/perl -w
#

use strict;

my $N = shift;
my $E = shift;

for(my $e = 0; $e < $E; $e++) {
    my @cnt;
    for(my $i = 0; $i < $N; $i++) {
	my $rn = rand();
	if ($rn <= 0.2) {
	    $cnt[0]++;
	} elsif ($rn <= 0.4) {
	    $cnt[1]++;
	} elsif ($rn <= 0.7) {
	    $cnt[2]++;
	} elsif ($rn <= 0.85) {
	    $cnt[3]++;
	} else {
	    $cnt[4]++;
	}
    }
    print join("\t", @cnt), "\n";
}
