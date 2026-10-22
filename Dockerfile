# Build from the workspace root:
# docker build -f scripts/CNNregR/Dockerfile -t cnnregr:0.1.0 .

# -----------------------------
# Builder stage
# -----------------------------
FROM rocker/r-ver:4.4.1 AS builder

ARG MAKEFLAGS="-j2"
ARG CNNREG_VERSION=0.1.2
ARG CNNREGR_REF=v0.1.0

ENV DEBIAN_FRONTEND=noninteractive \
    MAKEFLAGS=${MAKEFLAGS} \
    CONDA_DIR=/opt/conda \
    RENV_PROJECT=/project/CNNregR \
    RENV_PATHS_ROOT=/opt/renv \
    PATH=/opt/conda/bin:${PATH}

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
    libxml2-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e 'install.packages("renv", repos = "https://cloud.r-project.org")'

# Install Conda and create the pinned CNNreg environment.
ENV CONDA_VERSION=25.9.1-3
RUN wget --quiet \
      "https://repo.anaconda.com/miniconda/Miniconda3-py39_${CONDA_VERSION}-Linux-x86_64.sh" \
      -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p "${CONDA_DIR}" \
    && rm /tmp/miniconda.sh \
    && conda clean --all --yes

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

COPY scripts/CNNreg/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml \
    && conda run --no-capture-output -n cnnreg_env \
       python -m pip install --no-cache-dir --no-deps "CNNreg==${CNNREG_VERSION}" \
    && conda clean --all --yes

# Copy only the CNNregR renv project files before restoring dependencies.
WORKDIR /project/CNNregR
COPY scripts/CNNregR/.Rprofile ./.Rprofile
COPY scripts/CNNregR/renv ./renv
COPY scripts/CNNregR/renv.lock ./renv.lock

RUN R -q -e 'renv::restore(prompt = FALSE)' \
    && R -q -e "renv::install('mwang159/CNNregR@${CNNREGR_REF}', prompt = FALSE)"

RUN R -q -e 'library(CNNregR); stopifnot(packageVersion("CNNregR") == "0.1.0")' \
    && conda run --no-capture-output -n cnnreg_env \
    python -c 'from importlib.metadata import version; import torch; print("CNNreg", version("CNNreg")); print("PyTorch", torch.__version__)'

# -----------------------------
# Final stage
# -----------------------------
FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/conda \
    RENV_PROJECT=/project/CNNregR \
    RENV_PATHS_ROOT=/opt/renv \
    R_PROFILE_USER=/project/CNNregR/.Rprofile \
    PATH=/opt/conda/envs/cnnreg_env/bin:/opt/conda/bin:${PATH} \
    RETICULATE_PYTHON=/opt/conda/envs/cnnreg_env/bin/python

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libfontconfig1 \
    libfreetype6 \
    libfribidi0 \
    libgit2-1.5 \
    libharfbuzz0b \
    libicu72 \
    libjpeg62-turbo \
    libpng16-16 \
    libssl3 \
    libtiff6 \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e 'install.packages("renv", repos = "https://cloud.r-project.org")'

COPY --from=builder /project /project
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/renv /opt/renv

RUN rm -rf /opt/conda/pkgs/*

WORKDIR /project/CNNregR
CMD ["R"]
