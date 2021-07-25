#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  Rscript -e "rmarkdown::render('sim_lda.Rmd', params=list(id = $RUN, method='gradient', N=250, V=1000, n_models=10))"
  Rscript -e "rmarkdown::render('sim_lda.Rmd', params=list(id = $RUN, method='transport', N=250, V=1000, n_models=10))"
done;

tar -zcvf lda-$id.tar.gz lda/
