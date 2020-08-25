# Bayes Multiple Measures

Code for Edwards et. al 2020.

Contact owner of this repository for questions regarding the analysis.  In order to run our analyses, you should install:

* `rstan` 2.21.2 (assuming you already have a Stan installation)

* `tidyverse` 1.3.0

* `tidybayes` 2.1.1

* and `cowplot` 1.0.0 and `patchwork` 1.0.1 to create some of our figures.

The folder structure is as follows:

* `analysis` contains R scripts to fit models, simulate coverage, and perform sensitivity analyses, and make figures.

* `stan` contains stan files for models as well as any existing compiled models.

* `fits` contains mode summary outputs from `analysis/fit_combined.R` and `analysis/fit_single.R`.

* `sensitivity` contains results from our sensitivity analysis.

* `coverage` contains results from our coverage simulation.

