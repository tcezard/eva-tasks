#!/bin/bash

study=$1
if ls /nfs/production3/eva/data/${1}/60_eva_public/*accessioned* 1> /dev/null 2>&1; then
    : # no operation, file exists
else
    echo "no accessioned files available in eva/data/ folder for study $1, (/nfs/production3/eva/data/${1}/60_eva_public/)" >&2
    exit 1
fi

submittedVariantAccession=$(zcat /nfs/production3/eva/data/${1}/60_eva_public/*accessioned* | grep -m 1 "^[^#]" | cut -f 3 | sed s/ss//)

is_a_number_regex='^[0-9]+$'
if ! [[ $submittedVariantAccession =~ $is_a_number_regex ]] ; then
    echo "error: couldn't find SS accession in eva/data/${1}/60_eva_public/ : $submittedVariantAccession" >&2
    exit 1
fi

curl -s https://www.ebi.ac.uk/eva/webservices/identifiers/v1/submitted-variants/${submittedVariantAccession}
