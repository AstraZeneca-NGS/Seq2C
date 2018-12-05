#!/bin/bash

# Given a bam and bed file, calculate the coverage
# The script is used to submit jobs to grid

SEQ2COV="/users/kdld047/work/NGS/NGS/seq2c/seq2cov.pl"
MINISEQ2COV="/users/kdld047/work/NGS/NGS/seq2c/miniseq2cov.sh"
THREADS=6
while getopts ":N:b:B:p:P:m:" opt; do
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
        MINISEQ2COV=$OPTARG
	;;
      P)
        THREADS=$OPTARG
	;;
      m)
        SEQ2COV=$OPTARG
	;;
      \?)
        echo "Invalid option: -$OPTARG"
	exit 1
	;;
    esac
done

if [ ! $SAMPLE -o ! $BAM -o ! $BED ]
  then
      echo "Usage: $0 -N sample_name -B BED_file -b bam_file -P processes\n"
fi

splitBed.pl $BED $THREADS
BEDBASE=`basename $BED`
for n in $(eval echo "{1..$THREADS}"); do
    $MINISEQ2COV -p $SEQ2COV -N $SAMPLE -b $BAM -B $BEDBASE -P $n &
done
waitVardict.pl seq2cov $THREADS
head -1 ${SAMPLE}_cov.txt.1 > ${SAMPLE}_cov.txt
cat ${SAMPLE}_cov.txt.[0-9]* | grep -v Length > ${SAMPLE}_cov.txt
rm $BEDBASE.* ${SAMPLE}_cov.txt.[0-9]* seq2cov.done.[0-9]*

