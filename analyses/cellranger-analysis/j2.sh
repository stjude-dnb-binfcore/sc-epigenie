#!/bin/bash

set -e
set -o pipefail

########################################################################
# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Read root path
rootdir=$(realpath "./../..")
echo "$rootdir"

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
metadata_dir=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'metadata_dir:' | awk '{print $2}')
metadata_dir=${metadata_dir//\"/}  # Removes all double quotes
echo "$metadata_dir"  # Output: This is a string with quotes.

cellranger_parameters=$(cat ${rootdir}/project_parameters.Config.yaml | grep 'cellranger_parameters:' | awk '{print $2}')
cellranger_parameters=${cellranger_parameters//\"/}  # Removes all double quotes
echo "$cellranger_parameters"  # Output: This is a string with quotes.


########################################################################
# Create directories to save output files to
module_dir=$rootdir/analyses/cellranger-analysis
echo "$module_dir"

results_dir=results-test4
echo "$results_dir"

logs_dir=$module_dir/$results_dir/01_logs
echo "$logs_dir"

########################################################################
# Path to the TSV file
SAMPLES_FILE="${metadata_dir}"/project_metadata.tsv

# Check if the TSV file exists
if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Error: TSV file not found at $SAMPLES_FILE"
  exit 1
fi


########################################################################
# STEP 2 ###############################################################
########################################################################
# Create directories to save output files to
mkdir -p ./$results_dir/03_cellranger_count_summary/${cellranger_parameters}

# Define the output CSV file where all data will be combined
output_file="${module_dir}/$results_dir/03_cellranger_count_summary/${cellranger_parameters}/QC_Summary_CellRanger_Report.csv"

# Read and process each sample from the SAMPLES_FILE
tail -n +2 "$SAMPLES_FILE" | sort -t$'\t' -k2,2 | cut -f 1-3 | while IFS=$'\t' read -r ID SAMPLE FASTQ; do
  # Skip the header row by checking if SAMPLE is "SAMPLE"
  if [ "$SAMPLE" == "SAMPLE" ]; then
    continue
  fi

  # Process each sample here
  echo "Processing Sample ID: $ID, Sample: $SAMPLE, FASTQ: $FASTQ"

  # Define the summary file for the current sample
  summary_file="${module_dir}/$results_dir/02_cellranger_count/${cellranger_parameters}/${ID}/outs/summary.csv"
  echo "$summary_file"

  # Check if the summary file exists
  if [ -f "$summary_file" ]; then
    # Extract the header from the summary.csv and save it to the output file (if not already added)
    if [ ! -f "$output_file" ]; then
      # Extract header from the first line of summary.csv
      header=$(head -n 1 "$summary_file")
      echo "$header" > "$output_file"  # Create the output file with the new header
      echo "Created new output file with header."
    fi

    # Extract the entire content starting from the second line (NR>1)
    awk 'NR>1 {print $0}' "$summary_file" >> "$output_file"

    # Check if the command was successful
    if [ $? -eq 0 ]; then
      echo "Completed successfully: Reading summary.csv table for sample $ID."
    else
      echo "Error: Failed to read summary.csv table for sample $ID."
    fi
  else
    echo "Warning: $summary_file does not exist. Skipping sample $ID."
  fi

done

########################################################################
# Save TSV FILE
# Define the output TSV file where all data will be combined
cat "${output_file}" | sed 's/,/\t/g' > "${results_dir}/03_cellranger_count_summary/${cellranger_parameters}/QC_Summary_CellRanger_Report.tsv"

# Delete the output file if it exists
rm -f "$output_file"
########################################################################
