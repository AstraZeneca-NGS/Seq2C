#!/usr/bin/perl -w
#

use strict;
while( <> ) {
    chomp;
    my ($sample, $bam) = split(/\t/);
    open(SAMH, "samtools idxstats $bam |" );
    my $cnt = 0;
    while( <SAMH> ) {
        my @a = split(/\t/);
	$cnt += $a[2];
    }
    close( SAMH );
    print "$sample\t$cnt\n";
}

=INFO

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

=cut
