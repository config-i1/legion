---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "15 May 2021"
output: html_document
---
## Version
This is an initial submission of ``legion`` package, v0.1.0.

I have fixed all the issues outlined by Julia Haider and by Gregor Seyer.

## Test environments
* local ubuntu 20.04, R 4.0.5
* win-builder (devel and release)
* rhub with rhub::check_for_cran() command

## R CMD check results
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## win-builder check results
>* checking package dependencies ... NOTE
>Package suggested but not available for checking: 'doMC'

This is expected, because doMC is not available for Windows.

>Possibly mis-spelled words in DESCRIPTION:
>  ETS (15:21, 16:57)
>  VES (16:89)
>  al (15:94)
>  de (15:82)
>  et (15:91)

Everything is correct.

>Found the following (possibly) invalid URLs:
>  URL: https://doi.org/10.2307/2533213
>    From: inst/doc/ves.html
>    Status: 403
>    Message: Forbidden

URL works and leads to the paper of Bedrick & Tsai (1994). Model Selection for Multivariate Regression in Small Samples.

## rhub checks
### Fedora Linux, R-devel, clang, gfortran
PREPERROR with "Error: Bioconductor version '3.13' requires R version '4.1'; R version is too new". Seems to be an issue with rhub rather than with the package.

### Ubuntu Linux 20.04.1 LTS, R-release, GCC; Windows Server 2008 R2 SP1, R-devel, 32/64 bit
>Possibly mis-spelled words in DESCRIPTION:
>  ETS (15:21, 16:57)
>  VES (16:89)
>  al (15:94)
>  de (15:82)
>  et (15:91)

Everything is correct.

>* checking installed package size ... NOTE
>  installed size is  5.2Mb
>  sub-directories of 1Mb or more:
>    libs   4.4Mb

Does not look important.

### Debian Linux, R-devel, GCC ASAN/UBSAN
PREPERROR due to "ERROR: compilation failed for package ‘forecast‘". Not clear, what is wrong with forecast package.

### Windows Server 2008 R2 SP1, R-devel, 32/64 bit
>* checking package dependencies ... ERROR
>Package suggested but not available: 'doMC'

doMC is not available for Windows.

## Downstream dependencies
This is a new package, so there are no reverse dependencies.
