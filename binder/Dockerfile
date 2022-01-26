FROM rocker/tidyverse:4.1.2

RUN apt-get update
RUN apt-get install -y libgsl-dev vim libxml2-dev libxt6
RUN apt-get install -y libssl-dev libcurl4-gnutls-dev
RUN Rscript -e "install.packages('Barycenter')"
RUN Rscript -e "install.packages('MCMCpack')"
RUN Rscript -e "install.packages('philentropy')"
RUN Rscript -e "install.packages('scales')"
RUN Rscript -e "install.packages('superheat')"
RUN Rscript -e "install.packages('topicmodels')"

RUN Rscript -e "install.packages('devtools')"
RUN apt-get install -y git
RUN git clone https://github.com/lasy/alto.git
RUN Rscript -e "devtools::install('alto')"
