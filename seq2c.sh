#!/bin/bash

# Usage: seq2c.sh sample2bam.txt bed control_samples seq2copt sge_opt

SAM2BAM=$1
BED=$2

CONTROL=$3  # Optional control sample names. For multiple controls, separate them using :
SEQ2COPT=$4

SGE_OPT=$5 # e.g. "-q batch.q"

SAMTOOLS_OPT='samtools'
if [ $6 ]
    then
    SAMTOOLS_OPT=$6
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

OPT=""
if [ $CONTROL ]
    then
    CONS="-c $CONTROL"
    OPT="-c"
fi

if [ -n "$SGE_OPT" ]; then
  COUNTER=0
  while read i; do 
    let COUNTER=COUNTER+1
    a=(${i//\\t/}) 
    qsub $SGE_OPT -pe smp 1 -cwd -V -N seq2c_sam${COUNTER} -S /bin/bash ${DIR}/seq2cov_wrap.sh ${a[1]} ${a[0]} $BED $COUNTER ${DIR}/seq2cov.pl $SAMTOOLS_OPT
  done < $SAM2BAM
  perl ${DIR}/waitVardict.pl seq2c $COUNTER
  cat cov.txt.* > cov.txt
else
  cat $SAM2BAM | while read i; do a=(${i//\\t/}); perl ${DIR}/seq2cov.pl -z -b ${a[1]} -N ${a[0]} $BED; done > cov.txt
fi
perl ${DIR}/bam2reads.pl $SAM2BAM > read_stats.txt

#echo cov2lr.pl -a $CONS read_stats.txt cov.txt lr2gene.pl $OPT 
perl ${DIR}/cov2lr.pl -a $CONS read_stats.txt cov.txt | perl ${DIR}/lr2gene.pl $OPT $SEQ2COPT > seq2c_results.txt



