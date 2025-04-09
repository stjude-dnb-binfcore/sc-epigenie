#!/bin/bash

set -e
set -o pipefail

########################################################################
# Load modules
#module load python/3.9.9
module load cellranger-atac/2.1.0
########################################################################

# Set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Read root path
rootdir=$(realpath "./../..")
echo "$rootdir"

########################################################################
# Read multiple values and assign them to variables by parsing yaml file
metadata_dir=$(cat ../../project_parameters.Config.yaml | grep 'metadata_dir:' | awk '{print $2}')
metadata_dir=${metadata_dir//\"/}  # Removes all double quotes
echo "$metadata_dir"  # Output: This is a string with quotes.

genome_reference_path=$(cat ../../project_parameters.Config.yaml | grep 'genome_reference_path:' | awk '{print $2}')
genome_reference_path=${genome_reference_path//\"/}  # Removes all double quotes
echo "$genome_reference_path"  # Output: This is a string with quotes.

cellranger_parameters=$(cat ../../project_parameters.Config.yaml | grep 'cellranger_parameters:' | awk '{print $2}')
cellranger_parameters=${cellranger_parameters//\"/}  # Removes all double quotes
echo "$cellranger_parameters"  # Output: This is a string with quotes.


########################################################################
# Create directories to save output files to
mkdir -p ./results/01_logs
mkdir -p ./results/02_cellranger_count
mkdir -p ./results/02_cellranger_count/${cellranger_parameters}

#mkdir -p ./input

module_dir=$rootdir/analyses/cellranger-analysis
echo "$module_dir"

logs_dir=$module_dir/results/01_logs
echo "$logs_dir"

########################################################################
# If your `project_metadata` is not in `*.txt` file format
# use the following code line to convert it
#cat "${metadata_dir}"/project_metadata.tsv | sed 's/,/\t/g' > ./input/project_metadata.txt

# Path to the TSV file
SAMPLES_FILE="${metadata_dir}"/project_metadata.tsv

# Check if the TSV file exists
if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Error: TSV file not found at $SAMPLES_FILE"
  exit 1
fi

########################################################################
# Single-Library Analysis with cellranger-atac count
# https://software.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count
# To generate single cell accessibility counts for a single library
#  #--output-dir=${module_dir}/results/02_cellranger_count/${cellranger_parameters}/ \

cd ./results/02_cellranger_count/${cellranger_parameters}

# Loop through each line of the TSV file (skipping the header)
# Using `cut` to extract the first three columns (ID, SAMPLE, FASTQ)
# This assumes your TSV is tab-separated. If it's space-separated, adjust the delimiter accordingly.
# Sort the samples by the second column (SAMPLE) alphabetically, and process
tail -n +2 "$SAMPLES_FILE" | sort -t$'\t' -k2,2 | cut -f 1-3 | while IFS=$'\t' read -r ID SAMPLE FASTQ; do
  # Skip the header row by checking if SAMPLE is "SAMPLE"
  if [ "$SAMPLE" == "SAMPLE" ]; then
    continue
  fi

  # Process each sample here
  echo "Processing Sample ID: $ID, Sample: $SAMPLE, FASTQ: $FASTQ"

  # Run cellranger-atac count for each sample
  echo "Running cellranger-atac count for sample $SAMPLE..."
   # Submit job to LSF with appropriate options
  bsub -n 8 -M 64000 -R "rusage[mem=64000]" -o "${logs_dir}/${ID}.out" -e "${logs_dir}/${ID}.err" \
    cellranger-atac count \
      --id="${ID}" \
      --reference="${genome_reference_path}" \
      --fastqs="${FASTQ}" \
      --sample="${SAMPLE}" \
      --localcores=8 \
      --localmem=64 \
      --jobmode=lsf


  # Check if the command was successful
  if [ $? -eq 0 ]; then
    echo "cellranger-atac count completed successfully for sample $SAMPLE."
  else
    echo "Error: cellranger-atac count failed for sample $SAMPLE."
  fi

done