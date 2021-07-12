#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 2; ++k )); do
  export RUN=$( expr 2 '*' "$id" + "$k")
  for S in $(seq 10 10 200); do
    Rscript -e "rmarkdown::render('sim_equivalence.Rmd', params=list(rep = $RUN, S=$S))"
  done;
done;

tar -zcvf equivalence-$id.tar.gz equivalence/
