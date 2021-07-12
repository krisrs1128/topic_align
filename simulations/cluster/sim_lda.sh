#!/usr/bin/env bash

cp /srv/rlibs/alto/doc/*Rmd .

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  Rscript -e "rmarkdown::render('sim_lda.Rmd', params=list(id = $RUN))"
done;

tar -zcvf lda-$id.tar.gz lda/
