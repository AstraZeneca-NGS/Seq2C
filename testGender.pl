#!/usr/bin/env perl
# Generates gender data file. Can be used in cov2lr.pl with -G option to provide gender data for samples
use FindBin;
use lib "$FindBin::RealBin/libraries/";
use Stat::Basic;
use Statistics::TTest;
use Getopt::Std;
use strict;

our ($opt_B, $opt_d, $opt_H, $opt_x, $opt_L, $opt_b, $opt_N, $opt_n, $opt_y, $opt_h);
our $VERSION = '1.4.0';

getopts('hHyB:d:x:L:b:N:n:') || USAGE();
$opt_H && USAGE();

my $DIR = $opt_d ? $opt_d : "";
my $COVERAGE = $opt_x ? $opt_x : 10;
my $READLEN = $opt_L ? $opt_L : 100;
my $stat = new Stat::Basic;
my $ttest = new Statistics::TTest;
my @bedY= (); #arrays for chrX, chrY and all other regions
my @bedX= ();
my @bedA= ();
my $chrflag = 0;
my $BEDFile = $opt_B ? $opt_B : ($DIR ? "$DIR/hg19_gender.bed" : "hg19_gender.bed"); #Bed file contains gender data
open(BED, $BEDFile);
my $BEDSIZE = 0;
#Read bed data
while( <BED> ) {
    my @a = split;
    $BEDSIZE += $a[2] - $a[1];
    if ( $a[0] =~ /Y/ ) {
	push(@bedY, ["$a[0]:$a[1]-$a[2]", $a[2]-$a[1]]);
    } elsif ( $a[0] =~ /X/ ) {
	push(@bedX, ["$a[0]:$a[1]-$a[2]", $a[2]-$a[1]]);
    } else {
	push(@bedA, ["$a[0]:$a[1]-$a[2]", $a[2]-$a[1]]);
    }
}
close(BED);

my @sam2bam; #array of sample names and bam paths
if ( $opt_b ) {
    my $sample = `basename $opt_b .bam`; chomp $sample;
    $opt_n = "($opt_n)" if ( $opt_n && $opt_n !~ /\(/ );
    $sample = $1 if ( $opt_n && $opt_b =~ /$opt_n/ );
    $sample = $opt_N if ( $opt_N );
    push(@sam2bam,[$sample, $opt_b]);
} else {
    while( <> ) {
        my ($sample, $bam) = split;
	push(@sam2bam, [$sample, $bam]);
    }
}

@sam2bam = sort {$a->[0] cmp $b->[0]} @sam2bam; #sort by samples
#Output results
print join("\t", qw(Sample Gender chrY_dep chrX_log2 chrA_log2 p_value chrA/X_ratio)), "\n" if ( $opt_h );
foreach my $r (@sam2bam) {
    my ($sample, $bam) = @$r;
    my $total = totalRead($bam);
    my $depY = getcov($bam, \@bedY);
    print STDERR "Y: ", join(", ", @$depY), "\n" if ( $opt_y );
    my $depX = getcov($bam, \@bedX, 1);
    print STDERR "X: ", join(", ", @$depX), "\n" if ( $opt_y );
    my $depA = getcov($bam, \@bedA, 1);
    print STDERR "A: ", join(", ", @$depA), "\n" if ( $opt_y );
    my $medY = @$depY > 0 ? sprintf("%.3f", $stat->median($depY)) : 0; #median depth of Y chr
    my $gender = $medY > $COVERAGE ? "Male" : "Female";
    my $medA = sprintf("%.3f", $stat->median( $depA )); #median depth of autosomes
    my $medX = sprintf("%.3f", $stat->median( $depX )); #median depth of X chr
    my $p = "NA";
    if ( @$depA >= 3 && @$depX >= 3 ) {
	$ttest->load_data($depA, $depX);
	$p = sprintf("%.4f", $ttest->{ t_prob });
    }
    if ( $gender eq "Female" ) {
	if ( $p ne "NA" ) {
	    $gender = "X,-Y" if ( 2**($medA - $medX) > 1.5 && $p < 0.001 );
	} else {
	    $gender = "Unknown";
	}
    }
    print join("\t", $sample, $gender, $medY, $medX, $medA, $p, sprintf("%.2f", 2**($medA - $medX))), "\n";
}

#Calculate coverage for chromosomes (scaling log2 if flag is set)
sub getcov {
    my ($bam, $bedr, $flag) = @_;
    my @cov = ();
    foreach my $beda (@$bedr) {
	my ($bed, $size) = @$beda;
        if ( $chrflag ) {
	    $bed = "chr$bed" unless( $bed =~ /^chr/ );
	} else {
	    $bed =~ s/^chr// if( $bed =~ /^chr/ );
	}
	my $cnt = `samtools view $bam $bed | wc -l`;
	chomp $cnt;
	my $dep = sprintf("%.1f", $cnt * $READLEN / $size);
	if ( $flag ) {
	    $dep = log($dep+1)/log(2);
	    push(@cov, $dep) if ( $cnt > 5 ); # Ignore non-chrY regions without coverage (likely deleted)
	} else { 
	    push(@cov, $dep);
	}
    }
    return \@cov;
}

#Calculate total count of mapped reads
sub totalRead {
    my $bam = shift;
    open(CNT, "samtools idxstats $bam |");
    my $total = 0;
    while( <CNT> ) {
	my @a = split;
        $total += $a[2]; #add mapped reads count
	$chrflag = 1 if ( /^chr/ );
    }
    close( CNT);
    return $total;
}

sub USAGE {
    print <<USAGE;
USAGE:
    $0 [-hHy] [-B gender_BED] [-d dir] [-x cov] [-L len] [-b bam] [-n regex] [-N name] sample2bam.txt
    Version: $VERSION
    The program will calculate the chrY coverage and make prediction of genders.  The BAM file need to be
    targeted with chrY genes, exome or WGS.  Otherwise, it might make wrong prediction.

    Prerequisite:
    1. Program 'samtools' need to be in your path.
    2. Need to install modules Stat::Basic and Statistics::TTest, which are included in seq2c package, and
       make them available in your PERL5LIB environment.
    3. The BED files.  hg19 and hg38 BED files are included.

    Output:
    The output contains 7 columns:
    1. The name of the sample
    2. Predicted gender.  Possible values are:
       2.1 Male	 It is likely a male
       2.2 Female	 It is likely a female
       2.3 X,-Y	 It is likely only one X chr, but without Y chr.
       2.4 Unknown	 The information is not enough to determine the gender
    3. ChrY median depth
    4. ChrX median depth
    5. ChrA (autosome) median depth
    6. p_value of t-test between chrX and chrA
    7. The ratio between chrA/chrX (0.75-1.25 for Female, 1.65-2.35 for Male, roughly)

    Limitations:
    1. The bed file used for gender test should be targeted for targeted panels. No issue for WXS or WGS.
    2. The chrY in bed file should be unique to chrY
    3. For pure samples (e.g. cancer cell lines), males lost chrY or females lost one chrX can't not be
       differentiated for type "X,-Y".  It should be less an issue for clinical sequencing, as there're 
       almost always have normal cells in the sample.
    4. For pure samples, males who lose chrY and but with amplified chrX can be mistaken as female.

    Options:
    -H Print this help usage.
    -h Print the header.
    -y Print intermediate results (debuging purpose only).
    -B chrY_BED
       The BED file with CDS of chrY specific genes. Default to hg19_gender.bed in seq2c base directory.

    -d dir
       The installed seq2c base directory, if it is not in the path

    -x coverage
       The expected depth of coverage.  Specify it unless the depth is > 10x.  Default: 10.

    -L read_length
       The read length.  Default: 100

    -b bam_file (optional)
       The bam file.  For single sample prediction only, and no need to have sample2bam.txt file.

    -n regex
       Used with -b option.  The regular expression to extract sample name

    -N string
       The sample name

    sample2bam.txt is a file of two columns (tab delimited), first one is the sample name, the 2nd one is the BAM file path.

EXAMPLES
    1. Given a BAM file aligned to hg19, sample.bam, determine the gender:

       testGender.pl -b sample.bam -B hg19_gender.bed -N sample_name

    2. If you have a file (samples.txt) with list of BAM files (tab delimited), such as:
       sample1	/path/sample1.bam
       sample2	/path/sample2.bam
       sample3	/path/sample3.bam

       Then run the following command, assuming aligned to hg19:

       cat samples.txt | testGender.pl -B hg19_gender.bed

AUTHOR
    Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
    exit(0);
}
