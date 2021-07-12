#!/usr/bin/env bash

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  cp /srv/rlibs/alto/doc/$script .
  Rscript -e "rmarkdown::render('sim_lda.Rmd', params=list(id = $RUN))"
done;

tar -zcvf lda-$id.tar.gz lda/
