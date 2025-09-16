#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Set up modules
module load R/4.4.0
#module load pandoc/1.19.2.1
module load pandoc/3.2
module load texlive/20240326 # to load LaTeX on a compute node

module load python/2.7.5 # required for macs2
module load  macs2/2.1.1 

################################################################################################################
# Run module
Rscript --vanilla run-integration-with-scrna-seq-data.R
################################################################################################################
