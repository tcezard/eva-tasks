#!/bin/bash

# obtain the info of all the studies, just once
#curl https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all > all_studies.txt

cat all_studies.txt | /nfs/production3/eva/software/jq/1.5/jq -r ".response[].result[] | select(.id==\"$1\") | .assemblyAccession"

