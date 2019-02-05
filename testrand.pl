#!/usr/bin/env perl
# Create matrix consists needed number of lines containg 5 random elements in each line. Sum of elements in each line
# will be the requested.
use Getopt::Std;
use strict;

our ($opt_h);
USAGE() unless( getopts( 'h:' ));
USAGE() if ( $opt_h );
my $N = shift;
my $E = shift;

for(my $e = 0; $e < $E; $e++) {
    my @cnt;
    push @cnt, (0) x (5); #Initialize array to avoid errors in 5 columns

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
sub USAGE {
    print STDERR <<USAGE;
    $0 sum_of_elements lines

    The program will create matrix consists needed number of lines containg 5 random elements in each line.
    Sum of elements in each line will be the requested.
    Arguments are:
    - sum_of_elements: sum of random elements in each line.
    - lines: number of lines

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
   exit(0);
}
