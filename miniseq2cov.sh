#!/bin/bash

# Usage: minicheckCov.sh bam sample bed N
SEQ2COV=/group/cancer_informatics/tools_resources/NGS/seq2c/seq2cov.pl
while getopts ":N:b:B:p:P:" opt; do
    case $opt in
      N)
        SAMPLE=$OPTARG
	;;
      b)
        BAM=$OPTARG
	;;
      B)
        BED=$OPTARG
	;;
      p)
        SEQ2COV=$OPTARG
	;;
      P)
        N=$OPTARG
	;;
    esac
done

$SEQ2COV -N $SAMPLE -b $BAM ${BED}.$N > ${SAMPLE}_cov.txt.$N
touch seq2cov.done.$N

