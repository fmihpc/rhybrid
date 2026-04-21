#!/usr/bin/env bash
output_file="resource_usage.txt"
date >$output_file 2>&1
pwd >>$output_file 2>&1
echo "" >>$output_file 2>&1
jobid=""
# try to find a slurm-*.out file
file=$(ls slurm-*.out 2>/dev/null | head -n 1)
if [[ -n "$file" ]]; then
    # extract N from slurm-N.out
    jobid=$(echo "$file" | sed -n 's/^slurm-\([0-9]\+\)\.out$/\1/p')
elif [[ -f jobid.txt ]]; then
    # read first line
    jobid=$(head -n 1 jobid.txt)
else
    echo "slurm.*.out nor jobid.txt not found, cannot determine job ID"
    exit 0
fi
# check jobid is numberic
[[ $jobid =~ ^[0-9]+$ ]] || {
    echo "error: invalid job ID: $jobid" >&2
    exit 1
}
echo "Job ID: $jobid" >>$output_file 2>&1
echo "" >>$output_file 2>&1
echo "seff -d $jobid" >>$output_file 2>&1
echo "" >>$output_file 2>&1
seff -d $jobid >>$output_file 2>&1
echo "" >>$output_file 2>&1
echo "sacct --long -j $jobid" >>$output_file 2>&1
echo "" >>$output_file 2>&1
sacct --long -j $jobid >>$output_file 2>&1
