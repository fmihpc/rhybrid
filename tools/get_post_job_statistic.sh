#!/usr/bin/env bash
set -euo pipefail
trap 'echo "SCRIPT ERROR: line $LINENO: $BASH_COMMAND" >&2' ERR
SECONDS=0
script_name="${BASH_SOURCE[0]##*/}"
script_path="$(realpath -- "${BASH_SOURCE[0]}")"

# ============================================================
# Determine Slurm job ID
# ============================================================

# Use job ID from command-line argument if provided.
jobid="${1:-}"

if [[ -z "$jobid" ]]; then
    matches=()

    # Search for slurm-N.out files.
    for file in slurm-*.out; do
        if [[ -f "$file" && "$file" =~ ^slurm-([0-9]+)\.out$ ]]; then
            matches+=("${BASH_REMATCH[1]}")
        fi
    done

    if (( ${#matches[@]} == 1 )); then
        jobid="${matches[0]}"

    elif (( ${#matches[@]} > 1 )); then
        echo "ERROR: Multiple slurm-N.out files found." >&2
        echo "Please specify a job ID explicitly." >&2
        exit 1

    elif [[ -f jobid.txt ]]; then
        IFS= read -r jobid < jobid.txt || true

    else
        echo "ERROR: Neither slurm-N.out nor jobid.txt found." >&2
        exit 1
    fi
fi

# Validate job ID.
if [[ ! "$jobid" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Invalid job ID: $jobid" >&2
    exit 1
fi

# ============================================================
# Initialize report
# ============================================================

output_file="resource_usage_${jobid}.txt"
separator="======================================================================================="

{
    echo "Time stamp: $(date -R)"
    echo
    echo "Script file: $script_path"
    echo
    echo "START: $script_name"
    echo
    echo "Current folder: $(pwd)"
    echo "Job ID: $jobid"
} > "$output_file"

# ============================================================
# Function: Execute command and log its output
# ============================================================

run_section() {
    local title="$1"
    local status
    shift

    {
        echo
        echo
        echo "$separator"
        echo "BEGIN: $title"
        echo "$separator"

        if "$@"; then
            :
        else
            status=$?
            echo "WARNING: $title failed with exit status $status"
        fi

        echo "$separator"
        echo "END: $title"
        echo "$separator"

    } >> "$output_file" 2>&1
}

# ============================================================
# Collect Slurm statistics
# ============================================================

run_section "sacct: basics" \
    sacct \
    --format=JobID,State,Elapsed,Submit,Start,End,AllocCPUS,ReqMem \
    -j "$jobid"

run_section "sacct: details" \
    sacct \
    --format=JobID,State,ExitCode,NTasks,MaxRSS,AveRSS,MaxRSSNode,MaxRSSTask \
    -j "$jobid"

run_section "seff" \
    seff -d "$jobid"

run_section "sacct --long" \
    sacct --long -j "$jobid"

run_section "sacct --format=ALL" \
    sacct --format=ALL -j "$jobid"

run_section "scontrol show job" \
    scontrol show job "$jobid"

# ============================================================
# Finish report
# ============================================================

{
    echo
    echo "FINISHED: $script_name"
    echo
    echo "Time stamp: $(date -R)"
    elapsed=$SECONDS
    printf 'Elapsed time: %02d:%02d:%02d\n' "$((elapsed / 3600))" "$(((elapsed % 3600) / 60))" "$((elapsed % 60))"
    echo
} >> "$output_file"

echo "Report saved to: $output_file"


