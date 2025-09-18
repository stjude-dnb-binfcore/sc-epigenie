FROM rocker/rstudio:4.4.0


# -------------------------------------------------------------------
# Environment setup
# -------------------------------------------------------------------
# Set default locale
ENV LC_ALL=C

# Ensure virtual env is on the PATH so we're using the right venv binaries
ENV VIRTUAL_ENV=/opt/container-venv
ENV PATH=/opt/container-venv/bin:$PATH

# Ensure virtual env is on the PATH for LSF TO load the Package with a Specific Library Path
# This fixes the issue with the cc libraries
ENV R_LIBS_USER=/home/user/R/x86_64-pc-linux-gnu-library/4.4

# -------------------------------------------------------------------
# Core system dependencies
# -------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    gdebi-core python3-pip python3-venv \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
    liblzma-dev libbz2-dev libglpk-dev libhdf5-dev \
    cmake openjdk-8-jre-headless fastqc curl git \
    build-essential autoconf automake gettext libtool \
    libgit2-dev samtools wget unzip tar \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libncurses5-dev libncursesw5-dev zlib1g-dev \
    r-cran-rjags jags time sudo perl \
    libgetopt-long-descriptive-perl libdigest-perl-md5-perl \
    python3-pygments \
    && rm -rf /var/lib/apt/lists/*
    
    
 
#
# Install Cloudflare certificates so Docker build isn't interrupted by
# certificate errors. Combination of certificates should be valid until 2029.
#

#COPY "./cloudflare_certs/original-certificate.pem" /usr/local/share/ca-certificates/2024-certificate.crt
#COPY "./cloudflare_certs/SJ-2025-2029-CF-WARP-certificate.pem" /usr/local/share/ca-certificates/2025-2029-certificate.crt
#RUN update-ca-certificates


# -------------------------------------------------------------------
# Python virtual environment
# -------------------------------------------------------------------
RUN python3 -m venv /opt/container-venv && \
    . /opt/container-venv/bin/activate
# Install system deps (Python + build essentials)
RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-pip \
    python3-venv \
    build-essential \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

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
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" \
    >> /usr/local/lib/R/etc/Rprofile.site

# Core R package managers
RUN Rscript -e 'install.packages(c("BiocManager","devtools","remotes","usethis"))'

# Pin Bioconductor version to 3.20
RUN Rscript -e 'BiocManager::install(version="3.20")'

# -------------------------------------------------------------------
# R packages: core dependencies
# -------------------------------------------------------------------
RUN Rscript -e 'BiocManager::install(c( \
  "celldex","infercnv","scDblFinder","numbat","miQC","scater", \
  "SingleCellExperiment","SingleR","ggbio","TxDb.Hsapiens.UCSC.hg19.knownGene", \
  "org.Hs.eg.db","org.Mm.eg.db", \
  "BSgenome.Hsapiens.UCSC.hg19","EnsDb.Hsapiens.v75", \
  "BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v86", \
  "BSgenome.Mmusculus.UCSC.mm10","EnsDb.Mmusculus.v79", \
  "JASPAR2020"), ask=FALSE, update=FALSE, verbose=TRUE)'

# This is for the sc-epigenie/library(enrichplot)
RUN Rscript -e 'BiocManager::install("enrichplot")'
RUN Rscript -e 'if (!requireNamespace("enrichplot", quietly = TRUE)) stop("enrichplot not installed!")'

# This is for the sc-epigenie/library(clusterProfiler)
RUN Rscript -e 'BiocManager::install("clusterProfiler")'
RUN Rscript -e 'if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("clusterProfiler not installed!")'

# Install TFBSTools + motifmatchr with explicit dependencies
RUN Rscript -e 'BiocManager::install(c( \
  "Biostrings","GenomicRanges","IRanges","S4Vectors","GenomeInfoDb","BiocGenerics", \
  "seqLogo","TFMPvalue","DirichletMultinomial","pwalign"), ask=FALSE, update=FALSE, verbose=TRUE)' && \
    Rscript -e 'BiocManager::install(c("TFBSTools","motifmatchr"), ask=FALSE, update=FALSE, verbose=TRUE)'

# -------------------------------------------------------------------
# GitHub packages (pinned where necessary)
# -------------------------------------------------------------------
RUN Rscript -e 'devtools::install_github("stuart-lab/signac", ref="develop")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-wrappers@community-vignette")' && \
    Rscript -e 'remotes::install_version("Seurat", version="4.4.0")'
RUN Rscript -e 'remotes::install_github("cole-trapnell-lab/monocle3")' && \
    Rscript -e 'remotes::install_github("cole-trapnell-lab/cicero-release", ref="monocle3")'
RUN Rscript -e 'devtools::install_github("welch-lab/RcppPlanc@5744d5f351f3dd1d45e4a59b9d2e4a0e8fdab6e7")'
RUN Rscript -e 'devtools::install_github("welch-lab/liger", dependencies=TRUE, force=TRUE)'
RUN Rscript -e 'devtools::install_github("YuLab-SMU/ggtree", dependencies=TRUE)'
RUN Rscript -e 'devtools::install_github("kharchenkolab/numbat", dependencies=TRUE)'
RUN Rscript -e 'devtools::install_github("SGDDNB/ShinyCell")'
RUN Rscript -e 'devtools::install_github("igordot/scooter")'

# -------------------------------------------------------------------
# Additional R CRAN packages
# -------------------------------------------------------------------
RUN Rscript -e 'install.packages(c( \
  "clustree","cowplot","data.table","flexmix","flextable","forcats","fs","future","GGally", \
  "ggh4x","ggplot2","ggpmisc","ggrepel","ggthemes","harmony","igraph","irlba", \
  "knitr","optparse","patchwork","purrr","RColorBrewer","remotes","reshape2","rlist","R.utils", \
  "SeuratObject","shiny","SoupX","stringr","tidytext","tidyverse","tinytex","yaml", \
  "shinyhelper","DT","ggdendro","ragg","pkgdown"))'
  

#
#
# -------------------------------------------------------------------
# Install texlive full version
# This step can break easily. The code is correct, user needs to keep bulding until it passes this step (very time consuming).
# I do not know why it breaks, but there must be some issue from Cloudflare that interferes with permissions and proper installation.
#
ENV TEXLIVE_INSTALL_NO_CONTEXT_CACHE=1 \
    NOPERLDOC=1

# Install TexLive prerequisites
RUN apt-get update && \
    apt-get install -y wget unzip tar \
    make fontconfig perl openjdk-8-jre libgetopt-long-descriptive-perl \
    libdigest-perl-md5-perl libncurses5 \
    python3-pygments && \
    rm -rf /var/lib/apt/lists/*

    
# Install TeX Live
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz -O /tmp/install-tl-unx.tar.gz && \
    tar -xzf /tmp/install-tl-unx.tar.gz -C /tmp && \
    cd /tmp/install-tl-* && \
    echo "selected_scheme scheme-basic" > texlive.profile && \
    echo "TEXDIR /usr/local/texlive/2024" >> texlive.profile && \
    ./install-tl --profile=texlive.profile && \
    rm -rf /tmp/install-tl-* /tmp/install-tl-unx.tar.gz

# Dynamically find tlmgr bin dir and add to PATH
RUN TEXLIVE_BIN=$(find /usr/local/texlive/2024/bin -type d -name "x86_64-*" | head -n1) && \
    echo "export PATH=$TEXLIVE_BIN:\$PATH" >> /etc/profile && \
    echo "PATH=$TEXLIVE_BIN:\$PATH" >> /etc/environment && \
    ln -s $TEXLIVE_BIN/tlmgr /usr/local/bin/tlmgr

ENV PATH="/usr/local/texlive/2024/bin/x86_64-linux:$PATH"
# -------------------------------------------------------------------
#
#



# -------------------------------------------------------------------
# Entrypoint (activates venv)
# -------------------------------------------------------------------
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

