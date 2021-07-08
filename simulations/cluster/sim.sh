#!/usr/bin/env bash

tar -zxvf topic_align.tar.gz
cd topic_align/simulations/

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  Rscript -e "rmarkdown::render('$script', params = list(id = $RUN))"
done;

mv *.csv $_CONDOR_SCRATCH_DIR/
