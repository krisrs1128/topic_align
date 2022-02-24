This repository contains code for reproducing the analysis of sections 4 and 6 in [Multiscale Analysis of Count Data through Topic Alignment](https://arxiv.org/abs/2109.05541). All experiments call the [alto](<https://lasy.github.io/alto>) package, which exports the functions necessary for topic alignment.

#### Repository Structure

This repository is organized as follows,
* `simulations`
  - `cluster`: Scripts to launch parallel runs of the `lda.Rmd`, `background-noise.Rmd`, and `strain-switching.Rmd` vignettes, each corresponding to one simulation section in the manuscript. They are designed to be launched in an [HTCondor](<)> high-throughput computing environment. The results of these simulations are directories of CSV files, each giving diagnostic measures for one experimental configuration. Saved versions of these computer cluster outputs are available at [this link](https://uwmadison.box.com/shared/static/accan0ji43d7l97bia2fa9ai5zu826xm.gz).
  - `many_ldas.Rmd`: Simulates multiple replicates of alignments from both true LDA and null multinomial models, respecitvely. Supplementary Figures 4 - 9 are created using this analysis.
  - `simulation_functions.R, visualization_functions.R`: Define the generative mechanisms and visualization helpers used within the simulation experiments, but which are not generic enough to include within `alto`.
  - `visualize-lda.Rmd, visualize-background.Rmd, visualize-switching.Rmd`: Output Figures 3 - 5 in the main manuscript, visualizing the CSV files output by the scripts in `cluster`. A number of additional visualizations are given, summarizing the dependence of alignment properties on experimental and model parameters. These Rmarkdown documents have a YAML parameter specifying the path where the simulation data (the tarball downloaded from the link above) have to be stored.
* `notes`: More detailed summaries, calculations, and experiments written over the course of the project, but which didn't end up in the final manuscript.
* `Dockerfile`: Specifies a Docker image that install all the packages needed to reproduce the analyses in this repository. The pre-built image can be viewed [here](https://hub.docker.com/r/krisrs1128/alignment).

#### Help
The authors welcome questions about this code either through an [issue](https://github.com/krisrs1128/topic_align/issues) within this repository or an email.
