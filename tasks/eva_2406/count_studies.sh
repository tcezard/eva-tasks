#!/bin/bash

set -e

INPUT=$1
COUNT_VAR=${INPUT%.csv}_count_variants.txt
COUNT_STUDY=${INPUT%.csv}_count_study.txt
COUNT_ASSEMBLY=${INPUT%.csv}_count_assembly.txt

# Remove the end of the timestamp to only keep the date
# Then count the number of varant per study.
# The if statement is here to account for the fact that some dbSNP studies have commas in there name.
awk -F 'T' 'BEGIN{OFS="T"} {NF--; print }' $INPUT | \
awk -F ',' '{if (NF>5){study=$3; for (i=4; i<NF-1; i++){study=study","$i};$3=study; $5=$NF}; a[$1"\t"$2"\t"$3"\t"$5]++} END{for (v in a){print v"\t"a[v]}}' > $COUNT_VAR