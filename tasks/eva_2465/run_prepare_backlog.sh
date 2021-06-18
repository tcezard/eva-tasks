set -e
ELOAD=$1
cd /nfs/production3/eva/submissions/ELOAD_${ELOAD}
/nfs/production3/eva/software/eva-submission/production/bin/prepare_backlog_study.py --eload ${ELOAD}
/nfs/production3/eva/software/eva-submission/production/bin/validate_submission.py --eload ${ELOAD} --report
