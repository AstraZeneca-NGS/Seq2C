#!/bin/bash
WORKSPACE="$HOME/IdeaProjects/Seq2C/tests/"

DIR_BAM="$WORKSPACE/bam"
BAM_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00306/exome_alignment/HG00306.mapped.ILLUMINA.bwa.FIN.exome.20120522.bam"
BAM=$(echo $BAM_URL | sed 's#.*/##')
BAM_PATH="$DIR_BAM/$BAM"

BED="bed/hs37d5_gender_part.bed"
SAMPLE2BAM="samples/HG00306_gender_sample.txt"
OUTPUT="seq2c_results.txt"
COV="cov.txt"

#---
# Download data files
#---
cd $DIR_BAM
# BAM and BAI download
echo Downloading BAM HG00306
if [ ! -f "$BAM_PATH" ]; then
    wget -nc $BAM_URL -O $BAM_PATH || true
    wget -nc $BAM_URL.bai -O $BAM_PATH.bai || true
fi

cd ..
echo Run HG00306 test without control sample.
seq2c.sh $SAMPLE2BAM $BED

#Compare results for gender prediction
echo Compare results of HG00306 gender prediction
ACTUAL_HG00306_GENDER="gender_predicted.txt"
EXPECTED_HG00306_GENDER="expected/expected_HG00306_gender.txt"
DIFF_FILE_HG00306_GENDER="diff_HG00306_gender.txt"

if diff <(sort $EXPECTED_HG00306_GENDER) <(sort $ACTUAL_HG00306_GENDER) > $DIFF_FILE_HG00306_GENDER
then
	echo "OK: no differences in gender prediction";
else
	echo "ERROR: files have differences in gender prediction!"
	exit 1;
fi

#Compare results for seq2c
echo Compare results of HG00306
ACTUAL_HG00306="seq2c_results_HG00306.txt"
mv $OUTPUT $ACTUAL_HG00306
ACTUAL_COV_HG00306="cov_HG00306.txt"
cp $COV $ACTUAL_COV_HG00306
EXPECTED_HG00306="expected/expected_seq2c_HG00306.txt"
EXPECTED_COV_HG00306="expected/expected_cov_HG00306.txt"
DIFF_FILE_HG00306="diff_HG00306.txt"
DIFF_FILE_COV_HG00306="diff_cov_HG00306.txt"

if diff <(sort $EXPECTED_COV_HG00306) <(sort $ACTUAL_COV_HG00306) > $DIFF_FILE_COV_HG00306
then
	echo "OK: no differences in coverage results";
else
	echo "ERROR: files have differences in coverage results!"
	exit 1;
fi

if diff <(sort $EXPECTED_HG00306) <(sort $ACTUAL_HG00306) > $DIFF_FILE_HG00306
then
	echo "OK: no differences in seq2c results";
else
	echo "ERROR: files have differences in seq2c results!"
	exit 1;
fi
