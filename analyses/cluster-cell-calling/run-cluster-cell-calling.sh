#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 


# Set up modules
module load R/4.4.0
#module load pandoc/1.19.2.1
module load pandoc/3.2
module load texlive/20240326 #to load LaTeX on a compute node

module load python/3.9.9


# Activate virtual environment
python -m venv py39


################################################################################################################
# Run module
Rscript --vanilla run-cluster-cell-calling-step1.R
Rscript --vanilla run-cluster-cell-calling-step2.R
################################################################################################################
