# hash:sha256:260debda408b10eca068c2a6b2000d414e2d42ad12decc53dd66d242eec191f9
FROM registry.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libpoppler-cpp-dev=0.62.0-2ubuntu2.14 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("dplyr", "1.1.4")' \
    && Rscript -e 'remotes::install_version("haven", "2.5.4")' \
    && Rscript -e 'remotes::install_version("lmtest", "0.9-40")' \
    && Rscript -e 'remotes::install_version("lubridate", "1.9.3")' \
    && Rscript -e 'remotes::install_version("moments", "0.14.1")' \
    && Rscript -e 'remotes::install_version("pdftools", "3.4.0")' \
    && Rscript -e 'remotes::install_version("psych", "2.3.6")' \
    && Rscript -e 'remotes::install_version("purrr", "1.0.2")' \
    && Rscript -e 'remotes::install_version("readr", "2.1.5")' \
    && Rscript -e 'remotes::install_version("sandwich", "3.0-2")' \
    && Rscript -e 'remotes::install_version("sentimentr", "2.9.0")' \
    && Rscript -e 'remotes::install_version("stargazer", "5.2.3")' \
    && Rscript -e 'remotes::install_version("stringr", "1.5.1")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.1")' \
    && Rscript -e 'remotes::install_version("webshot", "0.5.5")'

RUN Rscript -e 'remotes::install_github( \
        "willschulz/rechat", \
        ref = "b68ee000e295cee1442c8a0d1fe3e73a92cae728")'

COPY postInstall /
RUN /postInstall
