#!/bin/bash

study=$1
if ls /nfs/ftp/pub/databases/eva/${1}/*accessioned* 1> /dev/null 2>&1; then
    : # no operation, file exists
else
    echo "no accessioned files available in FTP folder for study $1, (/nfs/ftp/pub/databases/eva/${1}/)" >&2
    exit 1
fi

submittedVariantAccession=$(zcat /nfs/ftp/pub/databases/eva/${1}/*accessioned* | grep -m 1 "^[^#]" | cut -f 3 | sed s/ss//)

is_a_number_regex='^[0-9]+$'
if ! [[ $submittedVariantAccession =~ $is_a_number_regex ]] ; then
    echo "error: couldn't find SS accession in /nfs/ftp/pub/databases/eva/${1}/ : $submittedVariantAccession" >&2
    exit 1
fi

curl -s https://www.ebi.ac.uk/eva/webservices/identifiers/v1/submitted-variants/${submittedVariantAccession}
