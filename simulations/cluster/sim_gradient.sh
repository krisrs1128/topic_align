#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 1; ++k )); do
  export RUN=$( expr 1 '*' "$id" + "$k")
  for alpha in $(seq 0.0 0.1 1); do
    Rscript -e "rmarkdown::render('sim_gradient.Rmd', params=list(id = $RUN, alpha=$alpha, method='product'))"
    Rscript -e "rmarkdown::render('sim_gradient.Rmd', params=list(id = $RUN, alpha=$alpha, method='transport'))"
  done;
done;

tar -zcvf gradient-$id.tar.gz gradient/
