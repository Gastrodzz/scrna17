# =============================================================================
# scRNA-seq Analysis Platform — production image
# =============================================================================
# Base: Bioconductor's official image. It ships a matching R + Bioconductor
# release, all the system libraries the Bioconductor stack needs, AND OpenBLAS
# (critical for Seurat matrix performance to match a Mac). Pin the tag to the
# Bioconductor release that matches your R version:
#   RELEASE_3_21  ->  R 4.5   (this project was developed on R 4.5.0)
# Change it if you upgrade R/Bioconductor.
FROM bioconductor/bioconductor_docker:RELEASE_3_21

# --- Extra system libraries a few CRAN/Bioc packages need on top of the base ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      libglpk-dev \
      libgdal-dev \
      libudunits2-dev \
      libcairo2-dev \
      libxt-dev \
      libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Fast, reproducible package installs via Posit Public Package Manager -----
# Prebuilt Linux binaries (no compiling) + a frozen snapshot for reproducibility.
RUN echo 'options(\
    repos = c(P3M = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"),\
    BioC_mirror = "https://packagemanager.posit.co/bioconductor",\
    Ncpus = max(1L, parallel::detectCores()),\
    timeout = 600)' >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /app

# Install R dependencies FIRST as a cached layer (only re-runs when the manifest
# changes). Uses the project's existing, curated dependency manifest.
COPY install_dependencies.R /app/install_dependencies.R
RUN Rscript /app/install_dependencies.R

# --- Copy the application (data/ and sessions/ are excluded via .dockerignore
#     and provided at runtime through a mounted volume so they persist) --------
COPY . /app

# The app writes its SQLite auth DB + per-project session RDS under /app/data.
# Mount a host volume there (see deploy/README.md) so data survives restarts.
RUN mkdir -p /app/data /app/sessions

EXPOSE 3838

# ShinyProxy (and plain `docker run -p 3838:3838`) connect to this port.
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
