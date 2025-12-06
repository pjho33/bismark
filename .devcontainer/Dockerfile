# ==============================================================================
# NGS Methylation Pipeline - Dev Container
# Python 3.11 + R 4.3 + Bioconductor + Bismark tools
# ==============================================================================

FROM rocker/r-ver:4.3.2

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------------------------
# 1. System Dependencies
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools
    build-essential \
    gcc \
    g++ \
    gfortran \
    make \
    cmake \
    # Python
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip \
    # R dependencies
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libhdf5-dev \
    # Bioinformatics tools
    fastqc \
    cutadapt \
    samtools \
    bowtie2 \
    # Utilities
    wget \
    curl \
    git \
    vim \
    htop \
    pigz \
    parallel \
    # Cleanup
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# 2. Python Setup (Miniconda)
# ------------------------------------------------------------------------------
ENV CONDA_DIR=/opt/conda
ENV PATH="${CONDA_DIR}/bin:${PATH}"

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniconda.sh \
    && conda update -n base -c defaults conda -y \
    && conda clean -afy

# Python packages
RUN pip install --no-cache-dir \
    pyyaml \
    pandas \
    numpy \
    biopython \
    multiqc

# ------------------------------------------------------------------------------
# 3. Bismark Installation
# ------------------------------------------------------------------------------
ENV BISMARK_VERSION=0.24.2
RUN cd /opt \
    && wget -q https://github.com/FelixKrueger/Bismark/archive/refs/tags/v${BISMARK_VERSION}.tar.gz \
    && tar -xzf v${BISMARK_VERSION}.tar.gz \
    && rm v${BISMARK_VERSION}.tar.gz \
    && ln -s /opt/Bismark-${BISMARK_VERSION}/bismark /usr/local/bin/ \
    && ln -s /opt/Bismark-${BISMARK_VERSION}/bismark_genome_preparation /usr/local/bin/ \
    && ln -s /opt/Bismark-${BISMARK_VERSION}/bismark_methylation_extractor /usr/local/bin/ \
    && ln -s /opt/Bismark-${BISMARK_VERSION}/deduplicate_bismark /usr/local/bin/ \
    && ln -s /opt/Bismark-${BISMARK_VERSION}/coverage2cytosine /usr/local/bin/

# Trim Galore
ENV TRIM_GALORE_VERSION=0.6.10
RUN cd /opt \
    && wget -q https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/${TRIM_GALORE_VERSION}.tar.gz \
    && tar -xzf ${TRIM_GALORE_VERSION}.tar.gz \
    && rm ${TRIM_GALORE_VERSION}.tar.gz \
    && ln -s /opt/TrimGalore-${TRIM_GALORE_VERSION}/trim_galore /usr/local/bin/

# ------------------------------------------------------------------------------
# 4. R Packages (Bioconductor + CRAN)
# ------------------------------------------------------------------------------
RUN R -e "install.packages(c('yaml', 'data.table', 'ggplot2', 'pheatmap', 'RColorBrewer', 'gtools', 'tidyverse', 'rlang', 'remotes'), repos='https://cloud.r-project.org/')"

# Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('bsseq', 'DSS', 'GenomicRanges', 'rhdf5', 'HDF5Array'), ask=FALSE, update=FALSE)"

# ------------------------------------------------------------------------------
# 5. User Setup
# ------------------------------------------------------------------------------
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# ------------------------------------------------------------------------------
# 6. Workspace Setup
# ------------------------------------------------------------------------------
WORKDIR /workspace

# Set default shell
SHELL ["/bin/bash", "-c"]

USER $USERNAME

CMD ["bash"]
