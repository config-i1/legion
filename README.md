# legion
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/legion)](https://cran.r-project.org/package=legion)
[![Downloads](http://cranlogs.r-pkg.org/badges/legion)](https://cran.r-project.org/package=legion)
[![R-CMD-check](https://github.com/config-i1/legion/actions/workflows/test.yml/badge.svg)](https://github.com/config-i1/legion/actions/workflows/test.yml)

The package _legion_ implements several multivariate models for purposes of forecasting.

![hex-sticker of the legion package for R](https://github.com/config-i1/legion/blob/master/man/figures/legion-web.png?raw=true)

Here is the list of the included functions:

1. ves - Vector Exponential Smoothing.
2. vets - Vector ETS with PIC taxonomy.
3. auto.vets - Automatic selection of restrictions for VETS.
4. sim.ves - simulates data from VES.
5. oves - occurrence state space vector exponential smoothing model.

Available methods:

1. AIC, BIC, AICc, BICc;
2. coefficients;
3. fitted;
4. forecast;
5. actuals;
6. logLik;
7. modelType - type of the estimated model;
8. nobs;
9. nparam - number of the estimated parameters in the model;
10. nvariate - number of series in the model;
11. residuals - the residuals of the model (et in case of additive and log(1+et) for the multiplicative ones);
12. rstandard, rstudent - standardised and studentised residuals;
13. outlierdummy - extracts outliers in the model and creates dummy variables for them;
14. plot - produces several plots for diagnostics purposes. See the documentation for plot.legion();
15. print;
16. sigma;
17. simulate;
18. summary;

## Installation

The stable version of the package is available on CRAN, so you can install it by running:
> install.packages("legion")

A recent, development version, is available via github and can be installed using "remotes" in R. First, make sure that you have remotes:
> if (!require("remotes")){install.packages("remotes")}

and after that run:
> remotes::install_github("config-i1/legion")