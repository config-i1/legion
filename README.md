# legion
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/legion)](https://cran.r-project.org/package=legion)
[![Downloads](http://cranlogs.r-pkg.org/badges/legion)](https://cran.r-project.org/package=legion)
[![ko-fi](https://ivan.svetunkov.ru/ko-fi.png)](https://ko-fi.com/G2G51C4C4)

The package _legion_ implements several multivariate models for purposes of forecasting.

Here is the list of the included functions:

1. ves - Vector Exponential Smoothing.
2. vets - Vector ETS with PIC taxonomy.
3. sim.ves - simulates data from VES.
4. oves - occurrence state space vector exponential smoothing model.

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
10. residuals - the residuals of the model (et in case of additive and log(1+et) for the multiplicative ones);
11. plot - produces several plots for diagnostics purposes. See the documentation for plot.legion();
12. print;
13. sigma;
14. simulate;
15. summary;

## Installation

<!-- The stable version of the package is available on CRAN, so you can install it by running:
> install.packages("legion") -->

A recent, development version, is available via github and can be installed using "devtools" in R. First, make sure that you have devtools:
> if (!require("devtools")){install.packages("devtools")}

and after that run:
> devtools::install_github("config-i1/legion")