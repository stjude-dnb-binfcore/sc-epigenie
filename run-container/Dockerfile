# -------------------------------------------------------------------
# Base image
# -------------------------------------------------------------------
FROM rocker/rstudio:4.4.0
# Set default locale
ENV LC_ALL=C

# -------------------------------------------------------------------
# Environment variables for R and virtualenv
# -------------------------------------------------------------------
# Ensure virtual env is on the PATH so we're using the right venv binaries
ENV VIRTUAL_ENV=/opt/container-venv
ENV PATH=/opt/container-venv/bin:$PATH
# Ensure virtual env is on the PATH for LSF TO load the Package with a Specific Library Path
# This fixes the issue with the cc libraries
ENV R_LIBS_USER=/home/user/R/x86_64-pc-linux-gnu-library/4.4


# -------------------------------------------------------------------
# Install Cloudflare certificates so Docker build isn't interrupted by
# certificate errors. Combination of certificates should be valid until 2029.
#

#COPY "./cloudflare_certs/original-certificate.pem" /usr/local/share/ca-certificates/2024-certificate.crt
#COPY "./cloudflare_certs/SJ-2025-2029-CF-WARP-certificate.pem" /usr/local/share/ca-certificates/2025-2029-certificate.crt
#RUN update-ca-certificates
#RUN apt-get update && apt-get install -y --no-install-recommends  ca-certificates gnupg
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Python virtual environment
# -------------------------------------------------------------------
RUN apt-get -y update && apt-get -y install gdebi-core python3-pip python3-venv
RUN python3 -m venv /opt/container-venv
RUN . /opt/container-venv/bin/activate
# Avoids BiocManager errors in Rhtslib; "fatal error: lzma.h: No such file or directory" etc
RUN apt-get -y install liblzma-dev libbz2-dev

# Avoids BiocManager errors in bluster; "libglpk.so.40: cannot open shared object file: No such file or directory"
RUN apt-get -y install libglpk-dev

# Avoids BiocManager errors in infercnv: "libjags.so.4: cannot open shared object file: No such file or directory"
RUN apt-get update && \
    apt-get install -y r-cran-rjags

# Avoids errors in "Rscript -e 'devtools::install_github("welch-lab/RcppPlanc")'"
RUN apt-get -y install cmake libhdf5-dev

# Needed for pipeline stages
RUN apt-get update && \
    apt-get install -y openjdk-8-jre-headless fastqc curl

RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install Python deps
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
      pandas==2.2.3 \
      numpy==1.23.5 \
      scipy==1.14.1 \
      leidenalg \
      multiqc==1.25 \
      MACS2==2.2.9.1
    
 


# -------------------------------------------------------------------
# R setup
# -------------------------------------------------------------------
# RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

# -------------------------------------------------------------------
# install unix utils required by R libraries and not already installed above
# Needed for R libraries
RUN apt-get update && apt-get install -y libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxml2-dev \
    libgsl-dev \
    libgit2-dev \
    libudunits2-dev \
    libgdal-dev \
    gdal-bin \
    libhdf5-dev \
    g++ \
    git \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

## Core R package managers
RUN Rscript -e 'install.packages(c("BiocManager","devtools","usethis"))'

# Pin Bioconductor version to 3.20
RUN Rscript -e 'BiocManager::install(version="3.20")'

# -------------------------------------------------------------------
# R packages: core dependencies
# -------------------------------------------------------------------

# install dependencies from Guangchuang Yu's lab dependencies and libraries
RUN Rscript -e 'remotes::install_version("ggplot2", version = "3.5.0", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'install.packages("yulab.utils", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'install.packages("ggfun", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'install.packages("ggtangle", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'BiocManager::install("enrichplot", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("clusterProfiler", ask=FALSE, update=FALSE, verbose=TRUE)'

# Double-check after installation
# https://github.com/bioconda/bioconda-recipes/blob/master/recipes/bioconductor-enrichplot/meta.yaml
RUN Rscript -e 'if (!requireNamespace("enrichplot", quietly = TRUE)) stop("enrichplot not installed!")'
RUN Rscript -e 'if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("clusterProfiler not installed!")'


# -------------------------------------------------------------------

# install bioconductor packages
RUN Rscript -e 'install.packages("matrixStats", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'BiocManager::install("celldex", ask=FALSE, update=FALSE, verbose=TRUE)'


# -------------------------------------------------------------------
RUN Rscript -e 'BiocManager::install(c( \
  "infercnv","scDblFinder","numbat","miQC","scater", \
  "SingleCellExperiment","SingleR","ggbio","TxDb.Hsapiens.UCSC.hg19.knownGene", \
  "org.Hs.eg.db","org.Mm.eg.db", \
  "BSgenome.Hsapiens.UCSC.hg19","EnsDb.Hsapiens.v75", \
  "BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v86", \
  "BSgenome.Mmusculus.UCSC.mm10","EnsDb.Mmusculus.v79", \
  "JASPAR2020"), ask=FALSE, update=FALSE, verbose=TRUE)'

RUN Rscript -e 'BiocManager::install("TFBSTools", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("motifmatchr", ask=FALSE, update=FALSE, verbose=TRUE)'

# Double-check after installation
RUN Rscript -e 'if (!requireNamespace("TFBSTools", quietly = TRUE)) stop("TFBSTools not installed!")'
RUN Rscript -e 'if (!requireNamespace("motifmatchr", quietly = TRUE)) stop("motifmatchr not installed!")'


# -------------------------------------------------------------------
# GitHub packages (pinned where necessary)
# -------------------------------------------------------------------
RUN Rscript -e 'devtools::install_github("stuart-lab/signac", ref="develop", upgrade="never")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-wrappers@community-vignette", upgrade="never")' && \
    Rscript -e 'remotes::install_version("Seurat", version="4.4.0", upgrade="never")'
RUN Rscript -e 'remotes::install_github("cole-trapnell-lab/monocle3", upgrade="never")' && \
    Rscript -e 'remotes::install_github("cole-trapnell-lab/cicero-release", ref="monocle3", upgrade="never")'
#
#
# -------------------------------------------------------------------
# zlskidmore:update-docker
# RcppPlanc requires c++ headers for hdf5 which are not available in jammy by default, 
# I have modified things such that they are downloaded via an r-package, but in order to do so a specific cmake version is required, 
# which is installed also --no-install-recommends

# install newer cmake, needed for RcppPlanc
RUN apt-get update && apt-get install -y wget g++ make \
   && wget https://cmake.org/files/v3.28/cmake-3.28.3.tar.gz \
   && tar -zxvf cmake-3.28.3.tar.gz \
   && cd cmake-3.28.3 && ./bootstrap && make -j$(nproc) && make install \
   && cd .. && rm -rf cmake-3.28.3*
 
RUN Rscript -e 'devtools::install_github("theAeon/HighFive", upgrade="never")' # needed for RcppPlanc
RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc@5744d5f351f3dd1d45e4a59b9d2e4a0e8fdab6e7", upgrade="never")'
# https://github.com/welch-lab/RcppPlanc/issues/3
#RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc", force = TRUE, verbose = TRUE)'
# -------------------------------------------------------------------
#
#


RUN Rscript -e 'devtools::install_github("welch-lab/liger", dependencies=TRUE, force=TRUE, upgrade="never")'
#RUN Rscript -e 'devtools::install_github("YuLab-SMU/ggtree", dependencies=TRUE)'
#RUN Rscript -e 'devtools::install_github("kharchenkolab/numbat", dependencies=TRUE)'
RUN Rscript -e 'devtools::install_github("SGDDNB/ShinyCell",  upgrade="never")'
RUN Rscript -e 'devtools::install_github("igordot/scooter", upgrade="never")'


# -------------------------------------------------------------------
# Additional R CRAN packages
# -------------------------------------------------------------------
RUN Rscript -e 'install.packages(c( \
  "clustree", "cowplot", "data.table", "flexmix", "forcats", "fs", "future", "GGally", \
  "ggh4x", "ggpmisc", "ggrepel", "ggthemes", "harmony", "igraph", "irlba", \
  "knitr", "optparse", "patchwork", "purrr", "RColorBrewer", "remotes", "reshape2", "rlist", "R.utils", \
  "SeuratObject", "shiny", "SoupX", "stringr", "tidytext", "tidyverse", "tinytex", "yaml", \
  "shinyhelper", "DT", "ggdendro", "ragg", "pkgdown"), \
  repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'

RUN Rscript -e 'install.packages("flextable", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'if (!requireNamespace("flextable", quietly = TRUE)) stop("flextable not installed!")'


#
#
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Install texlive full version
# This step can break easily. The code is correct, user needs to keep bulding until it passes this step (very time consuming).
# I do not know why it breaks, but there must be some issue from Cloudflare that interferes with permissions and proper installation.
# -------------------------------------------------------------------
# zlskidmore:update-docker
# The texlive installation required some modification, it works on my end, however critically the cpu architecture matters here, I've compiled my image on an arm cpu, 
# if this image is to be used on intel chips you'll need to build this image on an intel machine and modify line 286 to point to an x86_64 path, 
# tex is automatically checking the cpu architecture during the image build and is only compiling for one based on that. 
# (note you can probably emulate x86 via docker if you need to).
# -------------------------------------------------------------------
ENV TEXLIVE_INSTALL_NO_CONTEXT_CACHE=1 \
    NOPERLDOC=1

# Install prerequisites
RUN apt-get update && \
    apt-get install -y wget unzip tar make fontconfig perl openjdk-8-jre \
    libgetopt-long-descriptive-perl libdigest-perl-md5-perl libncurses5 \
    python3-pygments && \
    rm -rf /var/lib/apt/lists/*

# Install TeX Live 2024 with full scheme
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz -O /tmp/install-tl-unx.tar.gz && \
    tar -xzf /tmp/install-tl-unx.tar.gz -C /tmp && \
    cd /tmp/install-tl-* && \
    echo "selected_scheme scheme-full" > texlive.profile && \
    echo "TEXDIR /usr/local/texlive/2024" >> texlive.profile && \
    ./install-tl --profile=texlive.profile && \
    rm -rf /tmp/install-tl-* /tmp/install-tl-unx.tar.gz

# Set PATH for TeX Live 2024
ENV PATH="/usr/local/texlive/2024/bin/x86_64-linux:$PATH"

# Optional: update tlmgr (not strictly needed for scheme-full, but still good practice)
RUN tlmgr option repository http://mirror.ctan.org/systems/texlive/tlnet && \
    tlmgr update --self

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#

#
# -------------------------------------------------------------------
# Run virtual end on container startup- adjust as needed.
# ENTRYPOINT . /opt/container-venv/bin/activate
# -------------------------------------------------------------------
#
#
