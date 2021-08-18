#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 1; ++k )); do
  export RUN=$( expr 1 '*' "$id" + "$k")
  for alpha in $(seq 0.0 0.1 1); do
    Rscript -e "rmarkdown::render('background-noise.Rmd', params=list(id = $RUN, alpha=$alpha, method='product', N=250, V=1000, n_models=10))"
    Rscript -e "rmarkdown::render('background-noise.Rmd', params=list(id = $RUN, alpha=$alpha, method='transport', N=250, V=1000, n_models=10))"
  done;
done;

tar -zcvf gradient-$id.tar.gz gradient/