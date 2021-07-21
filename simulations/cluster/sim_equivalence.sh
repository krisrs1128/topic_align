#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 1; ++k )); do
  export RUN=$( expr 2 '*' "$id" + "$k")
  for S in $(seq 10 20 200); do
    Rscript -e "rmarkdown::render('sim_equivalence.Rmd', params=list(rep = $RUN, S=$S, method='product'))"
    Rscript -e "rmarkdown::render('sim_equivalence.Rmd', params=list(rep = $RUN, S=$S, method='transport'))"
  done;
done;

tar -zcvf equivalence-$id.tar.gz equivalence/
