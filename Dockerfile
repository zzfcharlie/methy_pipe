FROM rocker/tidyverse:latest

RUN R -e "install.packages('BiocManager', repos = 'https://cran.rstudio.com/')"

RUN R -e "install.packages(c('data.table', 'dplyr', 'tibble', 'readxl', 'fmsb', 'scales', 'ggplot2', 'plotly', 'reticulate'), repos = 'https://cran.rstudio.com/')"

RUN R -e "BiocManager::install(c('EpiDISH', 'sesame', 'BiocParallel'), ask = FALSE)"

RUN R -e "BiocManager::install('preprocessCore', configure.args='--disable-threading', force = TRUE); BiocManager::install(c('EpiDISH', 'sesame', 'BiocParallel'), ask = FALSE)"

RUN R -e "install.packages('optparse', repos = 'https://cran.rstudio.com/')"

COPY ./dnaMethyAge.tar.gz /project/dnaMethyAge.tar.gz
COPY ./meffonym.zip /project/meffonym.zip

RUN R -e "devtools::install_local('/project/dnaMethyAge.tar.gz', dependencies = TRUE)"
RUN R -e "remotes::install_local('/project/meffonym.zip')"

COPY ./ExperimentHub /root/.cache/R/ExperimentHub

WORKDIR /project


ENTRYPOINT ["Rscript", "main.R"]

CMD ["--config=config.R", "--idat-dir=./data", "--sample-mapping-path=./data/sample_mapping_test.xlsx", "--sample-info-path=./data/sample_info_test.csv", "--output-dir=./output"]
