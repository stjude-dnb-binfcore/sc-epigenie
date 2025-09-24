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
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

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




# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Troubleshoot
# Converted the `meta.yaml` to `dockerfile`
# https://github.com/bioconda/bioconda-recipes/blob/master/recipes/bioconductor-enrichplot/meta.yaml

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(c('DOSE>=4.0.0,<4.1.0','ggtree>=3.14.0,<3.15.0','GOSemSim>=2.32.0,<2.33.0','aplot>=0.2.1','ggfun>=0.1.7','ggnewscale','ggplot2','ggrepel>=0.9.0','ggtangle>=0.0.4','igraph','magrittr','plyr','purrr','RColorBrewer','reshape2','rlang','scatterpie','yulab.utils>=0.1.6','enrichplot','clusterProfiler'), version='3.20', update=FALSE)"

# -------------------------------------------------------------------
# Test installation
# -------------------------------------------------------------------
RUN R -e "library(enrichplot); library(clusterProfiler)"
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


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
RUN Rscript -e 'devtools::install_github("YuLab-SMU/enrichplot")'

# Double-check after installation
RUN Rscript -e 'if (!requireNamespace("enrichplot", quietly = TRUE)) stop("enrichplot not installed!")'


# -------------------------------------------------------------------
# Additional R CRAN packages
# -------------------------------------------------------------------
RUN Rscript -e 'install.packages(c( \
  "clustree","cowplot","data.table","flexmix","flextable","forcats","fs","future","GGally", \
  "ggh4x","ggplot2","ggpmisc","ggrepel","ggthemes","harmony","igraph","irlba", \
  "knitr","optparse","patchwork","purrr","RColorBrewer","remotes","reshape2","rlist","R.utils", \
  "SeuratObject","shiny","SoupX","stringr","tidytext","tidyverse","tinytex","yaml", \
  "shinyhelper","DT","ggdendro","ragg","pkgdown"))'
  

RUN Rscript -e 'if (!requireNamespace("flextable", quietly = TRUE)) stop("flextable not installed!")'


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
    make fontconfig perl openjdk-8-jre \
    libgetopt-long-descriptive-perl libdigest-perl-md5-perl \
    libncurses5 python3-pygments \
    && rm -rf /var/lib/apt/lists/*

# Install TeX Live (scheme-basic)
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

# Install missing LuaLaTeX math packages
RUN tlmgr update --self && \
    tlmgr install unicode-math fontspec xunicode lualatex-math luaotfload

# Optional: update all packages to latest (can take time)
# RUN tlmgr update --all
# -------------------------------------------------------------------
#
#


# -------------------------------------------------------------------
# Run virtual end on container startup- adjust as needed.
# -------------------------------------------------------------------
ENTRYPOINT . /opt/container-venv/bin/activate


