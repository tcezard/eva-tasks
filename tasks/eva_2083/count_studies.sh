#!/bin/bash

set -e

INPUT=$1
COUNT_VAR=${INPUT%.csv}_count_variants.txt
COUNT_STUDY=${INPUT%.csv}_count_study.txt
COUNT_ASSEMBLY=${INPUT%.csv}_count_assembly.txt

awk -F ',' '{a[$1"\t"$2"\t"$3]++} END{for (v in a){print v"\t"a[v]}}' $INPUT > $COUNT_VAR

awk -F '\t' '{a[$1"\t"$2]+=$4; b[$1"\t"$2]++} END{for (v in a){print v"\t"b[v]"\t"a[v]}}' $COUNT_VAR > $COUNT_STUDY

awk -F '\t' '{a[$2]+=$4; b[$2]+=$3; c[$2]++} END{for (v in a){print v"\t"c[v]"\t"b[v]"\t"a[v]}}' $COUNT_STUDY > $COUNT_ASSEMBLY
