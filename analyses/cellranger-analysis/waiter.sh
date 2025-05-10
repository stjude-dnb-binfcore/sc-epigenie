#!/bin/bash

#BSUB -P run_cellranger
#BSUB -J waiter
#BSUB -q standard
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -o waiter.out
#BSUB -e waiter.err
#BSUB -cwd "."

########################################################################
# Read container path and file
rootdir=$(realpath "./../..")
echo "$rootdir" 

########################################################################
# Set up variables
jobdir="${rootdir}/analyses/cellranger-analysis"

########################################################################
# Define array with sample_prefix
mapfile -t sample_prefix < <(grep '^ *-' ${rootdir}/project_parameters.Config.yaml | sed 's/^ *- *//')

# Print the sample_prefix for debugging
echo "Extracted sample_prefix:"
for p in "${sample_prefix[@]}"; do
    echo "$p"
done

########################################################################
# Function to check if there are any running jobs with the title pattern `WT123`, `KNOCKOUT456`, etc.
check_jobs() {
    # Loop through each prefix and check for jobs
    for job_prefix in "${sample_prefix[@]}"; do
        pattern="${job_prefix}[0-9]+"
        # Check for jobs using the pattern
        if bjobs | grep -E "$pattern" > /dev/null; then
            echo "Found matching job: $pattern"
            return 1  # Jobs still running
        fi
    done
    return 0  # No matching jobs
}


# Loop until no matching jobs are found
while true; do
    if check_jobs; then
        echo "$(date +"%D"): $(date +"%T"): No matching jobs found. Submitting the job..." >> waiter.log
        # Submit your job here
        bsub -P summarize_CellRanger -q standard -n 1 -R "rusage[mem=2GB]" -R "span[hosts=1]" -J j2 -o ${jobdir}/j2.out -e ${jobdir}/j2.err ${jobdir}/j2.bsub
        break
    else
        echo "$(date +"%D"): $(date +"%T"): Matching jobs found. Sleeping for 10 seconds..." >> waiter.log
        sleep 10
    fi
done 