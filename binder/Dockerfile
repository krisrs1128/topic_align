FROM r-base:4.1.2

RUN apt-get update
RUN apt-get install -y libgsl-dev vim libxml2-dev
RUN apt-get install -y libssl-dev libcurl4-gnutls-dev pandoc
RUN Rscript -e "install.packages('Barycenter')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('magrittr')"
RUN Rscript -e "install.packages('MCMCpack')"
RUN Rscript -e "install.packages('philentropy')"
RUN Rscript -e "install.packages('purrr')"
RUN Rscript -e "install.packages('readr')"
RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('scales')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('superheat')"
RUN Rscript -e "install.packages('tibble')"
RUN Rscript -e "install.packages('tidyr')"
RUN Rscript -e "install.packages('topicmodels')"

RUN Rscript -e "install.packages('devtools')"
RUN apt-get install -y git
RUN git clone https://github.com/lasy/alto.git
RUN Rscript -e "devtools::install('alto')"
ENV R_HOME=/usr/bin/R/
