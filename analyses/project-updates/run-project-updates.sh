#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

################################################################################################################
# Run all scripts of the module
Rscript --vanilla run-project-updates.R
################################################################################################################
