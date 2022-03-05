#!/bin/bash
#
#SBATCH --job-name=alto-lda-asymptotic
#SBATCH --output=alto-lda-asymptotic.out
#SBATCH --error=alto-lda-asymptotic.err
#SBATCH --partition=hns,stat,normal
#SBATCH --qos=normal
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL


singularity exec alto.sif Rscript -e "rmarkdown::render('run_simulations_for_asymptotic_behavior.Rmd', params = list(min_N = 5, n_N = 10, n_iters = 50))"

tar -zcvf lda-asymptotic.tar.gz lda-asymptotic/*topics*
