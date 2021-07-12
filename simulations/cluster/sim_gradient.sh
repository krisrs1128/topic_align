#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 2; ++k )); do
  export RUN=$( expr 2 '*' "$id" + "$k")
  for alpha in $(seq 0.0 0.05 1); do
    Rscript -e "rmarkdown::render('sim_gradient.Rmd', params=list(id = $RUN, alpha=$alpha))"
  done;
done;

tar -zcvf gradient-$id.tar.gz gradient/
