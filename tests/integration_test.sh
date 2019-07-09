#!/bin/bash -eu
set -o pipefail

#-------------------------------------------------------------------------------
# IMPORTANT: Before run you have to set paths to the BAM files of 
#            GTL_16_5 and controls in samples/GTL_16_5_sample2bam.txt
#-------------------------------------------------------------------------------

OUTPUT="seq2c_results.txt"
OUTPUT_DIR="/tmp/seq2c_test"

SCRIPT_DIR=`cd "$(dirname $0)"; pwd -P`
SEQ2C_PATH="$SCRIPT_DIR/.."
TEST_DIR="$SCRIPT_DIR"
echo "SCRIPT_DIR=$SCRIPT_DIR"
echo "SEQ2C_PATH=$SEQ2C_PATH"

# Add seq2c to path
PATH="$PATH:$SEQ2C_PATH"

#---
# Compare results
# Arguments: test_name expected_resutls actual_results
#---
function diff_results {
	TEST_NAME="$1"
	EXPECTED="$2"
	RESULTS="$3"

	DIFF_FILE="$OUTPUT_DIR/test_$TEST_NAME.diff.txt"
	echo "Test $TEST_NAME: Compare results, creating diff file '$DIFF_FILE'"

	if diff <(sort $EXPECTED) <(sort $RESULTS) > $DIFF_FILE
	then
		echo "Test $TEST_NAME: OK"
	else
		echo "Test $TEST_NAME: ERROR"
		exit 1;
	fi

	test_zeros $TEST_NAME $RESULTS
}

#---
# Check for absence of bad genes and presence of top good genes with ratios in results
# Arguments: test_name expected_bad_genes expected_good_genes actual_results
#---
function diff_genes {
	TEST_NAME="$1"
	EXPECTED_BAD_GENES="$2"
	EXPECTED_TOP_GENES="$3"
	RESULTS="$4"

	# Create empty diff file
	DIFF_FILE="$OUTPUT_DIR/test_$TEST_NAME.diff.txt"
	> $DIFF_FILE
	echo "Test $TEST_NAME: Creating diff file '$DIFF_FILE'"

	# Check if there are no bad genes in results
	echo "Test $TEST_NAME: Compare bad genes, append differences to '$DIFF_FILE'"
	RESULTS_GENES="$OUTPUT_DIR/test_$TEST_NAME.genes.txt"
	# AWK to get list of unique genes from results and decrease time
	awk 'unique[$2]++ {print $2;}' $RESULTS > $RESULTS_GENES
	while read BAD_GENE; do
	    if grep -P "$BAD_GENE\n" $RESULTS_GENES > /dev/null; then
			echo "Bad gene presence in results: $BAD_GENE" >> $DIFF_FILE
		fi
	done < $EXPECTED_BAD_GENES

	# Check if there are all good top genes with ratios and lengths in results
	echo "Test $TEST_NAME: Compare top good genes, append differences to '$DIFF_FILE'"
	while read GOOD_GENE; do
		if grep -P "$GOOD_GENE" $RESULTS > /dev/null; then
			continue
		else
			echo "Good gene absence or another results: $GOOD_GENE" >> $DIFF_FILE
		fi
	done < $EXPECTED_TOP_GENES

	if [ ! -s $DIFF_FILE ]; then
		echo "Test $TEST_NAME: OK"
	else
		echo "Test $TEST_NAME: ERROR"
		exit 1;
	fi

	test_zeros $TEST_NAME $RESULTS
}

#---
# Check for zeros (there must be only "0" and no "-0" in the output)
# Arguments: test_name actual_results
#---
function test_zeros {
    TEST_NAME="$1"
    RESULTS="$2"
	if grep -P "\t-0\t" $RESULTS > /dev/null; then
	    echo "Test $TEST_NAME: ERROR: result has -0 values"
		exit 1;
    fi
}

#---
# Test 1
#---
function test_1 {
	echo "Test 1: Run sample1 test without control sample"
	SAMPLE2BAM="$TEST_DIR/samples/Colo_sample1_hg19_chr20_sample2bam.txt"
	BED="$TEST_DIR/bed/hg19_chr20_Colo_sample1.bed"
	ACTUAL="$OUTPUT_DIR/test_1.seq2c_results_sample1.txt"
	EXPECTED="$TEST_DIR/expected/expected_Colo_sample1.txt"

	$SEQ2C_PATH/seq2c.sh $SAMPLE2BAM $BED
	mv $OUTPUT $ACTUAL
	diff_results 1 "$EXPECTED" "$ACTUAL"
}

#---
# Test 2
#---
function test_2 {
	echo "Test 2: Run sample1 test with control sample"
	CONTROL_SAMPLE="Colo829.chr20.18.Sample1"
	SAMPLE2BAM="$TEST_DIR/samples/Colo_sample1_hg19_chr20_sample2bam.txt"
	BED="$TEST_DIR/bed/hg19_chr20_Colo_sample1.bed"
	ACTUAL="$OUTPUT_DIR/test_2.seq2c_results_sample1_control.txt"
	EXPECTED="$TEST_DIR/expected/expected_Colo_sample1_control.txt"

	$SEQ2C_PATH/seq2c.sh $SAMPLE2BAM $BED $CONTROL_SAMPLE
	mv $OUTPUT $ACTUAL
	diff_results 2 "$EXPECTED" "$ACTUAL"
}

#---
# Test 3
#---
function test_3 {
	echo "Test 3: Run GTL_16_5 test with control sample"
	CONTROL_SAMPLE="control1:control2:control3:control4:control5:control6:"
	SAMPLE2BAM="$TEST_DIR/samples/GTL_16_5_sample2bam.txt"
	BED="$TEST_DIR/bed/panel_az600_chr7.bed"
	ACTUAL="$OUTPUT_DIR/test_3.seq2c_results_GTL_16_5.txt"
	EXPECTED="$TEST_DIR/expected/expected_GTL_16_5.txt"

	$SEQ2C_PATH/seq2c.sh $SAMPLE2BAM $BED $CONTROL_SAMPLE
	mv $OUTPUT $ACTUAL
	diff_results 3 "$EXPECTED" "$ACTUAL"
}

#---
# Test 4 (Test for checking bad, good genes on 2 samples from dfj example seq2c-coverage).
# COVERAGE variable must be set to a path to /seq2c-coverage_2samples.tsv
#---
function test_4 {
	echo "Test 4: Run dfj test for 2 samples 2-444961 and 2-444962"

	READ_MAPPING="$TEST_DIR/coverage/seq2c-read_mapping.txt"
	COVERAGE="$TEST_DIR/coverage/seq2c-coverage_2samples.tsv"

	ACTUAL="$OUTPUT_DIR/test_4.seq2c_results_dfj_2samples.txt"
	EXPECTED_BAD="$TEST_DIR/expected/bad_genes_2samples.txt"
	EXPECTED_TOP="$TEST_DIR/expected/top_genes_2samples.txt"

	echo "Test 4: Running cov2lr.pl script"
	$SEQ2C_PATH/cov2lr.pl $READ_MAPPING $COVERAGE > cov2lr.results.txt
	echo "Test 4: Running lr2gene.pl script"
	$SEQ2C_PATH/lr2gene.pl cov2lr.results.txt > seq2c_results.txt

	mv $OUTPUT $ACTUAL
	diff_genes 4 "$EXPECTED_BAD" "$EXPECTED_TOP" "$ACTUAL"
}

#---
# Test 5 (Test for checking bad, good genes on all samples from dfj example seq2c-coverage).
# COVERAGE variable must be set to a path to /seq2c-coverage.tsv
#---
function test_5 {
	echo "Test 5: Run dfj test for all samples"

	READ_MAPPING="$TEST_DIR/coverage/seq2c-read_mapping.txt"
	COVERAGE="$TEST_DIR/coverage/seq2c-coverage.tsv"

	ACTUAL="$OUTPUT_DIR/test_5.seq2c_results_dfj_all_samples.txt"
	EXPECTED_BAD="$TEST_DIR/expected/bad_genes_all.txt"
	EXPECTED_TOP="$TEST_DIR/expected/top_genes_all.txt"

	echo "Test 5: Running cov2lr.pl script"
	$SEQ2C_PATH/cov2lr.pl $READ_MAPPING $COVERAGE > cov2lr.results.txt
	echo "Test 5: Running lr2gene.pl script"
	$SEQ2C_PATH/lr2gene.pl cov2lr.results.txt > seq2c_results.txt

	mv $OUTPUT $ACTUAL
	diff_genes 5 "$EXPECTED_BAD" "$EXPECTED_TOP" "$ACTUAL"
}

#---
# Test 6 (Test for exons with the same start coordinates).
# The end of sigseg must contain the end coordinate from the last exon in a list of segments:
# in this case end is 154115315 for chrX F8A1 gene.
#---
function test_6 {
	echo "Test 6: Run HG00306 sample for the case with the same exon starts in chrX F8A1"

	SAMPLE2BAM="$TEST_DIR/samples/HG00306_exon_starts.txt"
	BED="$TEST_DIR/bed/hs37d5_gender_exon_order.bed"
	ACTUAL="$OUTPUT_DIR/test_6.seq2c_results_HG00306_exons.txt"
	EXPECTED="$TEST_DIR/expected/expected_HG00306_exons.txt"

	$SEQ2C_PATH/seq2c.sh $SAMPLE2BAM $BED
	mv $OUTPUT $ACTUAL
	diff_results 6 "$EXPECTED" "$ACTUAL"
}

#---
# Main: Run all tests
#---
mkdir "$OUTPUT_DIR" || true

test_1
test_2
test_3
test_4
test_5
test_6

