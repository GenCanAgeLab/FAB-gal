FROM rocker/r-ver:4.5.2

# System dependencies
RUN apt-get update -y && apt-get install -y \
    make \
    pkg-config \
    libpcre2-dev \
    libbz2-dev \
    liblzma-dev \
    libicu-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    pandoc \
    zlib1g-dev \
    libjpeg-dev \
    libpng-dev \
    default-jdk \
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# rJava configuration
RUN R CMD javareconf

# R and renv setup
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/

# Change to workdir
WORKDIR /srv/shiny-server/

# Copy renv
COPY Shiny/renv /srv/shiny-server/renv
COPY Shiny/.Rprofile /srv/shiny-server/.Rprofile
COPY Shiny/renv.lock /srv/shiny-server/renv.lock

# Env configuration
ENV RENV_PATHS_CACHE=/root/.cache/R/renv
ENV RENV_CONFIG_CACHE_SYMLINKS=false
ENV HOME=/home/

# Run renv::restore()
RUN  R -e 'source("renv/activate.R")'
RUN --mount=type=cache,id=renv-cache,target=/root/.cache/R/renv R -e 'renv::restore()'

# Pre-download the Bio-Formats Java archive so it is baked into the image
RUN R -e 'library(RBioFormats)'

# Copy Shiny app
COPY Shiny/ij /srv/shiny-server/ij
COPY Shiny/www /srv/shiny-server/www
COPY Shiny/app.R /srv/shiny-server/app.R
COPY Shiny/Functions.R /srv/shiny-server/Functions.R
COPY Example_images/ /srv/Example_images/

# Expose and run
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3838)"]