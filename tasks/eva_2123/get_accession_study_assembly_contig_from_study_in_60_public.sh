
#!/bin/bash

study=$1

set -o pipefail

variant=$(./get_variant_from_study_in_60_public.sh $study)

if [ -z "$variant" ]
then
    echo "error_couldnt_fetch_variant"
    exit 1
fi

echo $variant | /nfs/production3/eva/software/jq/1.5/jq -r   '(.[0].accession|tostring) + " " + .[].data.projectAccession + " " + .[].data.referenceSequenceAccession + " " + .[].data.contig' 

