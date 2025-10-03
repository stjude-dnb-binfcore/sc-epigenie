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
      numpy==2.2.1 \
      scipy==1.14.1 \
      leidenalg \
      multiqc==1.25 \
      macs3==3.0.0a7

# -------------------------------------------------------------------
# R setup
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
#RUN Rscript -e 'BiocManager::install("DOSE", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggtree", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("clusterProfiler", ask=FALSE, update=FALSE, verbose=TRUE)'

# install bioconductor packages
RUN Rscript -e 'install.packages("matrixStats", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest")'
RUN Rscript -e 'BiocManager::install("celldex", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("infercnv", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("scDblFinder", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("numbat", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("miQC", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("scater", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("SingleCellExperiment", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("SingleR", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("ggbio", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("org.Hs.eg.db", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("org.Mm.eg.db", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("EnsDb.Hsapiens.v75", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("EnsDb.Hsapiens.v86", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("EnsDb.Mmusculus.v79", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("JASPAR2020", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("GOSemSim", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggfun", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggnewscale", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggplot2", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggrepel", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("ggtangle", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("igraph", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("magrittr", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("plyr", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("purrr", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("RColorBrewer", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("reshape2", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("rlang", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("scatterpie", ask=FALSE, update=FALSE, verbose=TRUE)'

# Install TFBSTools + motifmatchr with explicit dependencies
#RUN Rscript -e 'BiocManager::install("Biostrings", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("GenomicRanges", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("IRanges", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("S4Vectors", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("GenomeInfoDb", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("BiocGenerics", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("seqLogo", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("TFMPvalue", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("DirichletMultinomial", ask=FALSE, update=FALSE, verbose=TRUE)'
#RUN Rscript -e 'BiocManager::install("pwalign", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("TFBSTools", ask=FALSE, update=FALSE, verbose=TRUE)'
RUN Rscript -e 'BiocManager::install("motifmatchr", ask=FALSE, update=FALSE, verbose=TRUE)'

# -------------------------------------------------------------------
# Test installation
# -------------------------------------------------------------------
RUN Rscript -e "library(enrichplot); library(clusterProfiler)"
# -------------------------

# -------------------------------------------------------------------
# GitHub packages (pinned where necessary)
# -------------------------------------------------------------------
RUN Rscript -e 'devtools::install_github("stuart-lab/signac", ref="develop", upgrade="never")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-wrappers@community-vignette", upgrade="never")' && \
    Rscript -e 'remotes::install_version("Seurat", version="4.4.0", upgrade="never")'
RUN Rscript -e 'remotes::install_github("cole-trapnell-lab/monocle3", upgrade="never")' && \
    Rscript -e 'remotes::install_github("cole-trapnell-lab/cicero-release", ref="monocle3", upgrade="never")'

# install newer cmake, needed for RcppPlanc
RUN apt-get update && apt-get install -y wget g++ make \
   && wget https://cmake.org/files/v3.28/cmake-3.28.3.tar.gz \
   && tar -zxvf cmake-3.28.3.tar.gz \
   && cd cmake-3.28.3 && ./bootstrap && make -j$(nproc) && make install \
   && cd .. && rm -rf cmake-3.28.3*
 
RUN Rscript -e 'devtools::install_github("theAeon/HighFive", upgrade="never")' # needed for RcppPlanc
RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc@5744d5f351f3dd1d45e4a59b9d2e4a0e8fdab6e7", upgrade="never")'
RUN Rscript -e 'devtools::install_github("welch-lab/liger", dependencies=TRUE, force=TRUE, upgrade="never")'
#RUN Rscript -e 'devtools::install_github("YuLab-SMU/ggtree", dependencies=TRUE)'
#RUN Rscript -e 'devtools::install_github("kharchenkolab/numbat", dependencies=TRUE)'
RUN Rscript -e 'devtools::install_github("SGDDNB/ShinyCell",  upgrade="never")'
#RUN Rscript -e 'devtools::install_github("igordot/scooter", upgrade="never")'
#RUN Rscript -e 'devtools::install_github("YuLab-SMU/enrichplot")'

# Double-check after installation
RUN Rscript -e 'if (!requireNamespace("enrichplot", quietly = TRUE)) stop("enrichplot not installed!")'

# -------------------------------------------------------------------
# Additional R CRAN packages
# -------------------------------------------------------------------
RUN Rscript -e 'install.packages("clustree", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("cowplot", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("data.table", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("flexmix", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("flextable", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("forcats", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("fs", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("future", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("GGally", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("ggh4x", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
# Rscript -e 'install.packages("ggplot2", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("ggpmisc", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("ggrepel", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("ggthemes", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("harmony", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("igraph", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("irlba", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("knitr", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("optparse", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("patchwork", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("purrr", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("RColorBrewer", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("reshape2", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("rlist", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("R.utils", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("SeuratObject", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("shiny", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("SoupX", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("stringr", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("tidytext", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("tidyverse", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("tinytex", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("yaml", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("shinyhelper", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("DT", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
RUN Rscript -e 'install.packages("ggdendro", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("ragg", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'
#RUN Rscript -e 'install.packages("pkgdown", repos="https://packagemanager.posit.co/cran/__linux__/jammy/latest", upgrade="never")'

# make sure flextable is installed
RUN Rscript -e 'if (!requireNamespace("flextable", quietly = TRUE)) stop("flextable not installed!")'


# -------------------------------------------------------------------
# Install texlive full version
# This step can break easily. The code is correct, user needs to keep bulding until it passes this step (very time consuming).
# I do not know why it breaks, but there must be some issue from Cloudflare that interferes with permissions and proper installation.
# -------------------------------------------------------------------
ENV TEXLIVE_INSTALL_NO_CONTEXT_CACHE=1 \
    NOPERLDOC=1

# Install TexLive prerequisites
RUN apt-get update && \
    apt-get install -y wget unzip tar \
    make fontconfig perl openjdk-8-jre \
    libgetopt-long-descriptive-perl libdigest-perl-md5-perl \
    libncurses5 python3-pygments \
    && rm -rf /var/lib/apt/lists/*

# Install TeX Live (scheme-basic)
WORKDIR /
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz -O /tmp/install-tl-unx.tar.gz && tar -xzf /tmp/install-tl-unx.tar.gz -C /tmp
RUN rm -f /tmp/install-tl-unx.tar.gz
RUN TEX_TMP=$(find /tmp/install-tl* -mindepth 0 -maxdepth 0 -type d) && \
    cd $TEX_TMP && \
    echo "selected_scheme scheme-basic" > texlive.profile && \
    echo "TEXDIR /usr/local/texlive/2025" >> texlive.profile && \
    cd $TEX_TMP && ./install-tl --profile=texlive.profile
RUN cd / && rm -rf /tmp/install-tl-*

# Add TeX Live binaries to PATH
ENV PATH="/usr/local/texlive/2025/bin/aarch64-linux/:$PATH"

# Install missing LuaLaTeX math packages
RUN tlmgr update --self && \
    tlmgr install unicode-math fontspec xunicode lualatex-math luaotfload

# set ENTRYPOINT for container
#ENTRYPOINT . /opt/container-venv/bin/activate

