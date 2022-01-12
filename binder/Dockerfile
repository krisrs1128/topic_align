FROM r-base:4.1.2

RUN apt-get update
RUN apt-get install -y libgsl-dev vim
RUN Rscript -e "install.packages('Barycenter')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('devtools')"
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

RUN git clone https://github.com/lasy/alto.git
RUN cd alto/vignettes/
RUN Rscript -e "devtools::install()"

