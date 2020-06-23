Seq2C
=====
Seq2C is Copy-Number Variations (CNV) variant caller. Seq2C normalizes the coverage from targeted sequencing to CNV log2ratio, 
and detects CNVs in regions that have abnormally high or low read depths compared to the rest.
## Getting started
Requirements:
1. Program 'samtools' need to be in your PATH.
2. Perl version 5.8 and higher.
3. You need to add PATH variables to Seq2C scripts, description in[this section](#add-path-variables-for-scripts).
### Getting source code
The Seq2C source code is located at [https://github.com/AstraZeneca-NGS/Seq2C](https://github.com/AstraZeneca-NGS/Seq2C).

To load the project, execute the following command:
```
git clone https://github.com/AstraZeneca-NGS/Seq2C.git
```

### Add PATH variables for scripts
Add short paths to the scripts seq2cov.pl, cov2lr.pl, lr2gene.pl, bam2reads.pl and others to your OS PATH. 
In Linux it can be done by adding this line to your `.bashrc` file: 

`export PATH=path_to_Seq2C_scripts_folder:$PATH`

For using plotting script, also add the plotting folder:  
`export PATH=path_to_Seq2C_scripts_folder/plotting:$PATH`
## Run Seq2C

Usage: 
```
seq2c.sh sample2bam.txt bed [control_sample1:[control_sample2]]
```
Parameters:
* `sample2bam.txt` - required, file must contains name of sample and paths to bam files with TAB delimiter (`\t`). Must contain
one empty line in end of file.
The example of the file:
```
NA12878.unmapped    /samples/NA12878.unmapped.bam
NA12889.mapped  /samples/NA12889.mapped.bam

```
* `bed` - required, path to BED file with TAB delimiter. Can be simple or amplicon. 
For simple mode it must contain 4 columns:
```
2   45904950    45904999    gene
```
For amplicon mode it must contain 8 columns. 7 column must contain amplicon start (must be more than start of region), 
8 column must contain amplicon end (must be less than end of region). 
```
2   45904950    45904999    gene    .   .   45904958    45904994
```
* `control_samples` - optional, for multiple controls, separate them using colon (`:`).

seq2c.sh script starts scripts from this repository in this order:
1. **seq2cov.pl** - Calculate candidate variance for a given region(s) in an indexed BAM file.
2. **bam2reads.pl** - Reads each read from sample and BAM file with samtools, creates statistics in `read_stats.txt` file.
3. **cov2lr.pl** - Normalizes the coverage from targeted sequencing to CNV log2 ratio.
4. **lr2gene.pl** - Converts a coverage file with CNV log2 ratio to copy number profile.

Output result (`seq2c_results.txt` file) contains columns:
1. **Sample** - sample name from original sample2bam file
2. **Gene** - gene name from BED file/region info 
3. **Chr** - chromosome name from BED file/region info 
4. **Start** - start of the segment/gene from BED file/region info 
5. **End** - end of the segment/gene from BED file/region info 
6. **Length** - length of the segment/gene
7. **Log2ratio** - CNV log2 ratio of normalized median depth by sample
8. **Sig** - significance (0 if AMP or DEL)
9. **BP_Whole** - type of CNV: Whole if AMP or DEL, else empty
10. **Amp_Del** - AMP (amplified) if log2 ratio more then -A option, DEL (deleted) if log2 ratio less then -D option 
11. **Ab_Seg** -  affected segments on gene
12. **Total_Seg** - total count of segments on gene
13. **Ab_log2ratio** - log2 ratio of normalized median depth by sample
14. **Log2r_Diff** - difference between CNV log2 ratio
15. **Ab_Seg_Loc** - segment location: ALL if AMP or DEL, else empty
16. **Ab_Samples** - abberation samples count for segments. Output only for BP calls that are “good”, i.e. wasn’t filtered by thresholds.
17. **Ab_Samples_Pcnt** - abberation samples fraction from total number of samples.

## Command line options for scripts
You can run scripts separately if you want transfer additional command (for example, gender files or tumor purities).
Parameters for scripts are listed in this chapter.
### seq2cov.pl
Usage:  
`seq2cov.pl [-hz] [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-o ori] [-d depth] region_info`

The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default input
 is IGV's one or more entries in refGene.txt, but can be regions (BED file).

Arguments are:
-    `region_info`  
    Required.  IGV's one or more entries in refGene.txt, but can be regions (BED file).
    
Options are:
-   `-h`  
    Print this help
-   `-a int:float`  
       Indicate that it's PCR amplicon based calling.  Each line in region_info represents a PCR amplicon (including primers).
       Two numbers in option are parameter to decide whether a particular read or pairs belongs to the amplicon.  
       First is the number of base pairs.  
       The second is the fraction of overlapped portion to the length of read or pairs.  
       If the edges of reads (paired for Illumina) are within defined base pairs to the edges of 
       amplicons and overlapped portion greater then the fraction, it's considered belonging to the amplicon.   
       Suggested values are: 10:0.95.  When given a 8 column amplicon format BED files, it'll be set to 10:0.95 automatically, but can still be overwritten by -a option.
-   `-n name_reg`  
       The regular expression to extract sample name from bam filename
-  `-N name`  
       Mutual exclusive to -n.  Set the sample name to name
-   `-b bam`  
       The indexed BAM file
-   `-c chr`  
       The column for chr
-   `-S start`  
       The column for region start, e.g. gene start
-   `-E end`  
       The column for region end, e.g. gene end
-   `-s seg_starts`  
       The column for segment starts in the region, e.g. exon starts
-   `-e seg_ends`  
       The column for segment ends in the region, e.g. exon ends
-   `-g gene`  
       The column for gene name
-   `-x num`  
       The number of nucleotide to extend for each segment, default: 0
-   `-z `  
       Indicate whether it's zero based numbering, default is 1-based

### bam2reads.pl 
The program produces a file containing number of mapped or sequenced reads for samples.   
Usage:  
`bam2reads.pl sample2bam.txt`

Arguments are:
* `sample2bam.txt` file must contains name of sample and paths to bam files with TAB delimiter (`\t`). 

### cov2lr.pl
Usage: 
`cov2lr.pl [-aH] [-c control] mapping_reads.txt coverage.txt`

The program will convert a coverage file to copy number profile.

Arguments are:
-    `mapping_reads.txt`  
    Required.  A file containing # of mapped or sequenced reads for samples.  At least two columns. 
    Can be created by `bam2reads.pl` script.
    First is the sample name, 2nd is the number of mapped or sequenced reads.
-   `coverage.txt`  
    The coverage output file from seq2cov.pl script.  Can also take from standard in or more than one file.

Options are:
-   `-a`  
    Indicate this is amplicon or exon based calling.  By default, it'll aggregate at gene level.
-   `-M`  
    Indicate to adjust the MAD when transforming the distribution.  Default: no, or just simple linear function.
    If not sure, don't use this option.  It might have better performance when cohort size is over 30.
-   `-c sample(s)`  
    Specify the control sample(s), if aplicable.  Multiple controls are allowed, which are separated by ":"
-   `-F double`  
   The failed factor for individual amplicons.  If (the 80th percentile of an amplicon depth)/(the global median depth)
   is less than the argument, the amplicon is considered failed and won't be used in calculation.  Default: 0.2.
-   `-G Gender`  
   Take a file of gender information.  Two columns, first is sample name, second is either M or F.  If not provided,
   the program will try to guess.
 -  `-Y chrYratio`  
   For gender testing, if chrY is designed.  Default: 0.15.  If chrY is carefully designed, such as Foundation's assay,
   default is good.  For exome, the number should be higher, such as 0.3.
-   `-Z`  
    Indicate to output the frozen_file and all parameters into file Seq2C.frozen.txt

### lr2gene.pl
The program will convert a coverage file to copy number profile.  The default parameters are designed for 
detecting such aberrations for high tumor purity samples, such as cancer cell lines.  
For clinical samples, many parameters need to be adjusted.  

Usage:   
`lr2gene.pl [-aPH] [-cy] [-F float] [-s min_amplicon_#] [-A float] [-D float] cov2lr.txt`

Arguments are:
-    `cov2lr.txt` 
    The coverage output file from cov2lr.pl script.  Can also take from standard in or more than one file.

Options are:
-    `-c`  
        Indicate that control sample is used for normalization
-    `-y`  
        Debugging mode
-    `-s int`  
       The minimum consecutive amplicons to look for deletions and amplifications.  Default: 1.  Use with caution
       when it's less than 3.  There might be more false positives.  Though it has been successfully applied with
       option "-s 1" and identified one-exon deletion of PTEN and TP53 that were confirmed by RNA-seq. 
         
For whole gene:
-    `-A float`
       Minimum log2 ratio for a whole gene to be considered amplified.  Default: 1.50
-    `-D float`
       Minimum log2 ratio for a whole gene to be considered deleted.  Default: -2.00   
       
For < 3 exons:
-    `-E float`
       Minimum mean log2 ratio difference for <3 exon deletion/amplification to be called.  Default: 1.25 
-    `-M float`
       When considering partial deletions less than 3 exons/amplicons, the minimum MAD value, in addition to -E,
       before considering it to be amplified or deleted.  Default: 10
-    `-d float`
       When considering >=3 exons deletion/amplification within a gene, the minimum differences between the log2 of two segments.
       Default: 0.5
-    `-p float (0-1)`
       The p-value for t-test when consecutive exons/amplicons are >= 3.  Default: 0.00001
       
For break point in the middle of the gene:
-    `-e float`
       When considering breakpoint in the middle of a gene, the minimum number of exons.  Default: 5       
-    `-t float`
       When considering breakpoint in the middle of a gene, the minimum differences between the log2 of two segments.
       Default: 0.4
-    `-P float (0-1)`
       The p-value for t-test when the breakpoint is in the middle with min exons/amplicons >= [-e].  Default: 0.000001

For cohort level aberrations:
-    `-R float (0-1)`
       If a breakpoint has been detected more than "float" fraction of samples, it's considered false positive and removed.
       Default: 0.03, or 3%.  Use in combination with -N
-    `-N int`
       If a breakpoint has been detected more than "int" samples, it's considered false positives and removed.
       Default: 5.  Use in combination with -R.
       
## Additional scripts
Here is the information about scripts that can be used to generate gender files or convert the results of Seq2C.
### seq2c2fm.pl
The program will parse seq2c output and make calls for each gene and output in the format compatible with OncoPrint.
It has the option to provide the purity so that log2ratio thresholds will be adjusted accordingly.  
By default, it calls genes that are homozygously deleted or amplified with >= 6 copies.

Usage:  
`seq2c2fm.pl [-g] [-e exons] [-n reg] [-N num] [-A amp] [-a amp] [-D del] [-d del] [-p purity_file] [-P purity] lr2gene_output`

Arguments are:
-  `lr2gene_output`  
    The output file from seq2c.sh script (or lr2gene.pl).  Can also take from standard in or more than one file.

Options:
-   `-k`  
    Print header
-   `-g`  
    Whether to output copy gains [4-5] copies.  Default: no.
-   `-p file`  
   A file contains the tumor purity for all samples.  Two columns, first is sample name, second is the purity in % or fraction [0-1].
-   `-P double`  
   The purity.  Default: 1 or 100%, as is for cell lines.  If set, all samples will assume to have the same purity.
-   `-n regex`  
   The regular expression to extract sample names.  Default: none.
-   `-N num`  
   If an breakpoint has been called in >= num of samples, it's deemed false positive.  Default: 5
-   `-e exons`  
   Minimum number of exons/amplicon.  Default: 1  
   
For whole gene:  
-   `-D log2ratio`  
   The log2ratio to determine that a gene is homozygously deleted.  Default: -2.0  
-   `-A log2ratio`  
   The log2ratio to determine that a gene is amplified.  Default: 1.45.  
   
For < 3 exons:  
-   `-d log2ratio`   
   The minimum log2ratio to determine that 1-2 exons are deleted. Should be lower than [-d] to reduce false positives. Default: -2.5
-   `-a log2ratio`  
   The minimum log2ratio to determine that 1-2 exons are amplified. Should be larger than [-a] to reduce false positives. Default: 1.75
   
For gains:
-   `-G Genes`  
   List of genes, seperated by ":", for which gain will be captured.  Default: MYC

Output:  
Some columns are filled with blanks, "NA" or "-" for compatibility with OncoPrint format.
1. Sample 
2. Always ""
3. Variant_type: always "copy_number_alteration"
4. Gene 
5. Always "NA"
6. Always "-"
7. Always "-"
8. Segment "chr:start"
9. Always "-"
10. Always "-"
11. Transf_LogRatio: transformed log_ratio as (2^LogRatio)*2 
12. Segments: for BP - affected segments of total segments, for Whole - both numbers are total segments.  
13. LogRatio 
14. Alteration in format "amplification/gain/loss"  
15. Always "-" 
16. Always "-" 
17. Always "-" 
18. Always "-" 
19. Always "-" 
20. Always "-" 
21. Always "-" 
22. Alteration in format "Amplification/Gain/Deletion"

### testGender.pl
The program will calculate the chrY coverage and make prediction of genders.  The BAM file need to be
targeted with chrY genes, exome or WGS.  Otherwise, it might make wrong prediction.  
Output can be used in cov2lr.pl script with -G option.  

Usage:  
`testGender.pl [-hHy] [-B gender_BED] [-d dir] [-x cov] [-L len] [-b bam] [-n regex] [-N name] sample2bam.txt`

Output:
The output contains 7 columns:
1. The name of the sample
2. Predicted gender.  Possible values are:
   2.1 Male	 It's likely a male
   2.2 Female	 It's likely a female
   2.3 X,-Y	 It's likely only one X chr, but without Y chr.
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
-   `-H` 
    Print this help usage.
-   `-h` 
    Print the header.
-   `-y`  
    Print intermediate results (debuging purpose only).
-   `-B chrY_BED`  
    The BED file with CDS of chrY specific genes. Default to hg19_gender.bed in seq2c base directory.
-   `-d dir`   
    The installed seq2c base directory, if it's not in the path
-   `-x coverage`  
    The expected depth of coverage.  Specify it unless the depth is > 10x.  Default: 10.
-   `-L read_length`  
    The read length.  Default: 100
-   `-b bam_file (optional)`  
    The bam file.  For single sample prediction only, and no need to have sample2bam.txt file.
-   `-n regex`  
    Used with -b option.  The regular expression to extract sample name
-   `-N string`  
    The sample name

EXAMPLES
1. Given a BAM file aligned to hg19, sample.bam, determine the gender:
  `testGender.pl -b sample.bam -B hg19_gender.bed -N sample_name`
2. If you have a file (samples.txt) with list of BAM files (tab delimited), then run the following command, 
assuming aligned to hg19:  
   `cat samples.txt | testGender.pl -B hg19_gender.bed`
   
### testrand.pl
The program will create matrix consists needed number of lines containg 5 random elements in each line.
Sum of elements in each line will be the requested.

Usage:  
`testrand.pl sum_of_elements lines`

Arguments are:  
-   `sum_of_elements`  
    Sum of random elements in each line.
-   `lines`  
    Number of lines in constructed matrix.


Written by Zhongwu Lai, AstraZeneca, Boston, USA
