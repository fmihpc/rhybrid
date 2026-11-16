#!/usr/bin/env bash
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
# check jobid is numeric
[[ $jobid =~ ^[0-9]+$ ]] || {
    echo "error: invalid job ID: $jobid" >&2
    exit 1
}

output_file="resource_usage_$jobid.txt"
echo "BEGIN : get_post_job_statistic.sh" >$output_file 2>&1
pwd >>$output_file 2>&1
date >$output_file 2>&1
echo "Job ID: $jobid" >>$output_file 2>&1

wait
echo "" >>$output_file 2>&1
echo "" >>$output_file 2>&1
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : sacct : basics        ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
sacct --format=JobID,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem -j $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END : sacct : basics          ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>1&
echo "" >>$output_file 2>1&
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : sacct : details       ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
sacct --format=JobID,State,ExitCode,NTasks,MaxRSS,AveRSS,MaxRSSNode,MaxRSSTask -j $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END : sacct : details         ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>1&
echo "" >>$output_file 2>1&
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : seff                  ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
seff -d $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END : seff                    ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>&1
echo "" >>$output_file 2>&1
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : sacct --long          ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
sacct --long -j $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END   : sacct --long          ============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>&1
echo "" >>$output_file 2>&1
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : sacct --format=ALL   =============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
sacct --format=ALL -j $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END   : sacct --format=ALL   =============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>&1
echo "" >>$output_file 2>&1
wait

echo "=======================================================================================" >>$output_file 2>&1
echo "============================ BEGIN : scontrol show        =============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
scontrol show job $jobid >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1
echo "============================ END   : scontrol show        =============================" >>$output_file 2>&1
echo "=======================================================================================" >>$output_file 2>&1

wait
echo "" >>$output_file 2>&1
echo "" >>$output_file 2>&1
wait

date >>$output_file 2>&1
echo "END   : get_post_job_statistic.sh" >>$output_file 2>&1
