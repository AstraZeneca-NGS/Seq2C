#!/usr/bin/env perl
# Parse a list of refseq and check CDS coverage
use warnings;
use FindBin;
use lib "$FindBin::RealBin/libraries/";
use Getopt::Std;
use strict;

our ($opt_h, $opt_b, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_e, $opt_g, $opt_x, $opt_z, $opt_N, $opt_a);
USAGE() unless( getopts( 'hzb:s:e:S:E:n:c:g:x:N:a:' ) );
USAGE() if ( $opt_h );

my $BAM = $opt_b; # the bam file
my $sample = $1 if ( $BAM =~ /([^\/.]+)[\/]*.bam/ );
if ( $opt_n ) {
    $sample = $1 if ( $BAM =~ /$opt_n/ ); # set sample name by regex from BAM file name if -n option set
}
$sample = $opt_N if ( $opt_N ); #sample name to set
my $c_col = $opt_c ? $opt_c - 1 : 2; #chromosome
my $S_col = $opt_S ? $opt_S - 1 : 6; #region start (gene)
my $E_col = $opt_E ? $opt_E - 1 : 7; #region end (gene)
my $s_col = $opt_s ? $opt_s - 1 : 9; #segment start (exon)
my $e_col = $opt_e ? $opt_e - 1 : 10; #segment end (exon)
my $g_col = $opt_g ? $opt_g - 1 : 12; #gene

my $SPLICE = $opt_x ? $opt_x : 0; #number of bases to extend on each segment (exon)
my %regions; #hash of gene on (chr, CDS)

while( <> ) {
    next if ( /^track/i );
    next if ( /^browser/i );
    next if ( /^#/ );
    next if ( /^Sample/ );
    chomp;
    my @A = split(/\t/, $_, -1); #split bed file by columns
    # Amplicon/exon case: if option -a set, or if there are 8 columns in BED and start of segment more than start of
    # region and end of segment less then end of region.
    if ( $opt_a || (@A == 8 && $A[6] =~ /^\d+$/ && $A[7] =~ /^\d+$/ && $A[6] > $A[1] && $A[7] < $A[2]) ) {
	$opt_a = "10:0.95" unless($opt_a); #set number of base pairs to 10 and overlap fraction to 0.95 if were not set
	my ($chr, $start, $end, $gene, $score, $strand, $istart, $iend) = @A; #fill values from BED
	$start++ && $istart++; # if ( $opt_z ); #increment starts of region and segment
	push(@{ $regions{ $gene }->{ CDS } }, [$start, $end, $istart, $iend]); #push data to gene about coding region to the CDS
	$regions{ $gene }->{ chr } = $chr; #push data to gene about chromosome
    } else {
    # Non-amplicon case: of there are 4 columns in BED file and start of region less then end of region.
    # Segment start and end will be the same as start and end of region.
	($c_col, $S_col, $E_col, $g_col, $s_col, $e_col) = (0, 1, 2, 3, 1, 2) if ( (! $opt_c) && @A == 4 && $A[1] =~ /^\d+$/ && $A[2] =~ /^\d+$/ && $A[1] <= $A[2] );
	my ($chr, $cdss, $cdse, $gene) = @A[$c_col,$S_col,$E_col,$g_col]; #fill values from BED
	my @starts = split(/,/, $A[$s_col]); #if there are few starts or ends in one line, split them and add data about
	# coding region to array of CDS
	my @ends = split(/,/, $A[$e_col]);
	my @CDS = ();
	for(my $i = 0; $i < @starts; $i++) {
	    my ($s, $e) = ($starts[$i], $ends[$i]);
	    next if ( $cdss > $e ); # not a coding exon
	    # TODO: check this, how we will have comma if region and segment are the same.
	    last if ( $cdse < $s ); # No more coding exon
	    $s = $cdss if ( $s < $cdss ); # adjust start and end of coding region
	    $e = $cdse if ( $e > $cdse );
	    $s -= $SPLICE unless ( $s == $cdss ); #extend segment with splice length if splitted start/end and segment start/end differ
	    $e += $SPLICE unless ( $e == $cdse );
	    $s += 1 if ( $opt_z ); #increment start of segment for zero based numering
	    push(@CDS, [$s, $e]); #push data about coding region to the CDS array
	}
	push( @{ $regions{ $gene }->{ CDS } }, @CDS); #push data to gene about array of coding regions to the CDS
	$regions{ $gene }->{ chr } = $chr; #push data to gene about chromosome
    }
}

print join("\t", "Sample", "Gene", "Chr", "Start", "End", "Tag", "Length", "MeanDepth"), "\n";
my $bamhdr = `samtools view -H $BAM`; #bam header
my $genome = $bamhdr =~ /SN:chr\d+\t/ ? "hg" : "grch"; #determine reference by chromosome name
while( my ($gene, $r) = each %regions ) {
    my $exoncov; #coverage of exon
    my $CDS = $r->{ CDS };
    my $chr = $r->{ chr };
    my $tchr = $chr;
    if ( $genome eq "hg" ) { #fix the cromosome name
	$tchr = "chr$tchr" unless( $tchr =~ /^chr/ );
    } else {
	$tchr =~ s/^chr// if( $tchr =~ /^chr/ );
    }
    my $total = 0;
    my $gene_length = 0;
    my ($gene_start, $gene_end) = (500000000, 0);
    for(my $i = 0; $i < @{ $CDS }; $i++) {
	my ($START, $END, $ISTART, $IEND) = @{$CDS->[$i]}; #get info about current coding region.
	# For non amplicon/exon case ISTART and IEND will be unitialized
	$gene_length += $END-$START+1;
	$gene_start = $START if ( $START < $gene_start ); # adjust start and end of gene
	$gene_end = $END if ( $END > $gene_end );
	open(SAM, "samtools view $BAM $tchr:$START-$END |");
	$exoncov = 0;
    #parse SAM file
	while( <SAM> ) {
	    my @a = split(/\t/);
	    next if ( $a[1] & 0x800 && /\tSA:Z:/ ); # ignore supplementary alignments
	    my $dir = $a[1] & 0x10 ? "-" : "+"; #determine direction (reverse -, forward +)
	    my $start = $a[3]; #set start to alignment start
	    my @segs = $a[5] =~ /(\d+)[MD]/g; #get matched and deletion lengths from cigar
	    my $end = $start-1; $end += $_ foreach(@segs); #calculate end position
	    # Amplicon case
	    if ( $opt_a ) {
		my ($dis, $ovlp) = split(/:/, $opt_a); #split for distance and overlap fraction
		($dis, $ovlp) = (10, 0.95) unless($dis && $ovlp);
		my ($segstart, $segend) = ($start, $end); #set segment start and end to start and end of crrent read (matched+deleted part)
		if ( $a[5] =~ /^(\d+)S/ && $dir eq "+" ) { #forward softclips
		    my $ts1 = $segstart-$1 > $START ? $segstart-$1 : $START;
		    my $te1 = $segend < $END ? $segend : $END;
		    next unless( abs($ts1-$te1)/($segend-$segstart+$1) > $ovlp); #skip if current overlap is more then overlap fraction
		} elsif ($a[5] =~ /(\d+)S$/ && $dir eq "-" ) { #reverse softclips
		    my $ts1 = $segstart > $START ? $segstart : $START;
		    my $te1 = $segend+$1 < $END ? $segend+$1 : $END;
		    next unless( abs($te1-$ts1)/($segend+$1-$segstart) > $ovlp); #skip if current overlap is more then overlap fraction
		} else {
		    if ($a[6] eq "=" && $a[8]) { #if read and mate have the same chromosome name and inferred insert size exists
			($segstart, $segend) = $a[8] > 0 ? ($segstart, $segstart+$a[8]-1) : ($a[7], $a[7]-$a[8]-1);
		    }
		    my $ts1 = $segstart > $START ? $segstart : $START;
		    my $te1 = $segend < $END ? $segend : $END;
		    next unless( (abs($segstart - $START) <= $dis && abs($segend - $END) <= $dis ) && abs(($ts1-$te1)/($segend-$segstart)) > $ovlp);
		}
		$exoncov++;
		$total++;
	    } else { #Non-amplicon case: alignment length calculates as a smallest length from (region, read length)
		my $alignlen = ($END > $end ? $end : $END) - ($START > $start ? $START : $start)+1;
		$exoncov += $alignlen;
		$total += $alignlen;
	    }
	}
	close( SAM );
	# Print result for each coding region of gene
	if ( $opt_a ) {
	    print join("\t", $sample, $gene, $chr, $START, $END, "Amplicon", $END-$START+1, $exoncov), "\n";
	} else {
	    print join("\t", $sample, $gene, $chr, $START, $END, "Amplicon", $END-$START+1, sprintf("%.2f", $exoncov/($END-$START+1))), "\n";
	}
    }
    # Print result for each gene
    if ( $opt_a ) {
	print join("\t", $sample, $gene, $chr, $gene_start, $gene_end, "Whole-Gene", $gene_length, $total), "\n";
    } else {
	print join("\t", $sample, $gene, $chr, $gene_start, $gene_end, "Whole-Gene", $gene_length, sprintf("%.2f", $total/$gene_length)), "\n";
    }
}

sub USAGE {
    print STDERR <<USAGE;
    $0 [-hz] [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-o ori] [-d depth] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default
    input is IGV's one or more entries in refGene.txt, but can be regions

    -h Print this help
    -a int:float
       Indicate that it's PCR amplicon based calling.  Each line in region_info represents a PCR amplicon (including primers).
       Two numbers in option are parameter to decide whether a particular read or pairs belongs to the amplicon. First is the 
       number of base pairs.  The second is the fraction of overlapped portion to the length of read or pairs.  If the edges of 
       reads (paired for Illumina) are within defined base pairs to the edges of amplicons and overlapped portion greater then
       the fraction, it's considered belonging to the amplicon.  Suggested values are: 10:0.95.  When given a 6 column amplicon
       format BED files, it'll be set to 10:0.95 automatically, but can still be overwritten by -a option.
    -n name_reg
       The regular expression to extract sample name from bam filename
    -N name
       Mutual exclusive to -n.  Set the sample name to name
    -b bam
       The indexed BAM file
    -c chr
       The column for chr
    -S start
       The column for region start, e.g. gene start
    -E end
       The column for region end, e.g. gene end
    -s seg_starts
       The column for segment starts in the region, e.g. exon starts
    -e seg_ends
       The column for segment ends in the region, e.g. exon ends
    -g gene
       The column for gene name
    -x num
       The number of nucleotide to extend for each segment, default: 0
    -z 
       Indicate whether it's zero based numbering, default is 1-based

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
   exit(0);
}
