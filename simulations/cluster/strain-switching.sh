#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 1; ++k )); do
  export RUN=$( expr 2 '*' "$id" + "$k")
  for S in $(seq 10 20 240); do
    Rscript -e "rmarkdown::render('strain-switching.Rmd', params=list(rep = $RUN, S=$S, method='product', N=250, V=1000, n_models=10, K=5, save = TRUE, perplexity = TRUE))"
    Rscript -e "rmarkdown::render('strain-switching.Rmd', params=list(rep = $RUN, S=$S, method='transport', N=250, V=1000, n_models=10, K=5, save = TRUE, perplexity = TRUE))"
  done;
done;

tar -zcvf equivalence-$id.tar.gz equivalence/
