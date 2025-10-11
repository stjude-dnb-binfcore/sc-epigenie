#!/bin/bash

set -e
set -o pipefail

########################################################################
pwd
#echo $PATH
#echo $LD_LIBRARY_PATH
#which python

########################################################################
# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 
pwd

# Create dir
mkdir -p reports

################################################################################################################
# Run other dependencies
python3 --version
fastqc --version
multiqc --version
pandoc --version
tex --version
#macs3 --version
macs2 --version

# Run R script with all R packages
Rscript -e "rmarkdown::render('01-run-test-packages.Rmd', clean = TRUE,
                              output_dir = file.path('./reports/'),
                              output_file = c(paste('Report-', 'run-test-packages', '-', Sys.Date(), sep = '')),
                              output_format = 'all')"

# Run R script with all lsf problematic R packages
Rscript -e "rmarkdown::render('02-run-lsf-problematic-packages.Rmd', clean = TRUE,
                              output_dir = file.path('./reports/'),
                              output_file = c(paste('Report-', 'run-lsf-problematic-packages', '-', Sys.Date(), sep = '')),
                              output_format = 'all')"

Rscript --vanilla 03-run-lsf-problematic-packages.R
