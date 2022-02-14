---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "14 February 2022"
output: html_document
---
## Version
This is an initial submission of ``legion`` package, v0.1.1.


## Test environments
* local ubuntu 20.04.3, R 4.1.2
* github actions
* win-builder (devel and release)
* rhub with rhub::check_for_cran() command


## R CMD check results
>* checking installed package size ... NOTE
>  installed size is  5.4Mb
>  sub-directories of 1Mb or more:
>    libs   4.5Mb
>
> 0 errors ✓ | 0 warnings ✓ | 1 note x


## win-builder check results
>Found the following (possibly) invalid URLs:
>  URL: https://doi.org/10.2307/2533213
>    From: inst/doc/ves.html
>    Status: 403
>    Message: Forbidden

URL works and leads to the paper of Bedrick & Tsai (1994). Model Selection for Multivariate Regression in Small Samples.


## rhub checks
### Fedora Linux, R-devel, clang, gfortran
PREPERROR with "Error: Bioconductor version '3.13' requires R version '4.1'; R version is too new". Seems to be an issue with rhub rather than with the package.

### Windows Server 2022, R-devel, 32/64 bit, Fedora Linux, R-devel, clang, gfortran
>Found the following (possibly) invalid DOIs:
>  DOI: 10.1177/1471082X0901000401
>    From: DESCRIPTION
>    Status: Service Unavailable
>    Message: 503

DOI is correct, the paper is available online.

>* checking installed package size ... NOTE
>  installed size is  5.4Mb
>  sub-directories of 1Mb or more:
>    libs   4.5Mb

Does not look important.

### Debian Linux, R-devel, GCC ASAN/UBSAN
PREPERROR due to "ERROR: compilation failed for package ‘forecast‘". Not clear, what is wrong with forecast package.

### Windows Server 2022, R-devel, 32/64 bit
>* checking package dependencies ... ERROR
>Package suggested but not available: 'doMC'

doMC is not available for Windows.


### Debian Linux, R-devel, GCC ASAN/UBSAN
> ERROR: compilation failed for package ‘Rcpp’

Rcpp failed to compile on this server for whatever reason

## Downstream dependencies
This is a new package, so there are no reverse dependencies.
