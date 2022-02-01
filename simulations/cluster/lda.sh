#!/usr/bin/env bash

cp /alto/vignettes/*.Rmd .

for (( k = 0; k < 10; ++k )); do
  export RUN=$( expr 10 '*' "$id" + "$k")
  Rscript -e "rmarkdown::render('lda.Rmd', params=list(id = $RUN, method='product', K=5, N=250, V=1000, n_models=10, save = TRUE, perplexity = TRUE))"
  Rscript -e "rmarkdown::render('lda.Rmd', params=list(id = $RUN, method='transport', K=5, N=250, V=1000, n_models=10, save = TRUE, perplexity = TRUE))"
done;

tar -zcvf lda-$id.tar.gz lda/
