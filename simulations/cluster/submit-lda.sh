#!/usr/bin/env bash

R CMD install alto

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  Rscript -e "rmarkdown::render('simulations/sim_lda.Rmd', params = list(id = $RUN))"

mv *.csv $_CONDOR_SCRATCH_DIR/
