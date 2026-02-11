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
RUN echo "options(renv.config.pak.enabled = FALSE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site

# Install renv
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version = "1.1.4")'

# Change to workdir
WORKDIR /srv/shiny-server/

# Copy renv config files
COPY Shiny/renv.lock renv.lock
COPY Shiny/.Rprofile .Rprofile
COPY Shiny/renv/activate.R renv/activate.R
COPY Shiny/renv/settings.json renv/settings.json

# Force not to use symlinks
ENV RENV_CONFIG_CACHE_SYMLINKS=false

# 3. Run renv::restore()
RUN --mount=type=cache,id=renv-cache,target=/root/.cache/R/renv R -e 'renv::restore()'

# 4. Copy app
COPY Shiny/ /srv/shiny-server/

# 5. Exponer y arrancar
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3838)"]