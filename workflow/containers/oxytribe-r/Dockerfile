FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('data.table', 'purrr', 'tibble'), \
    repos='https://cloud.r-project.org/', \
    Ncpus=parallel::detectCores())"

WORKDIR /workspace
