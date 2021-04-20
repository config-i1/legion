---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "20 April 2021"
output: html_document
---
## Version
This is an initial submission of ``legion`` package, v0.1.0.

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


## rhub checks
### Fedora Linux, R-devel, clang, gfortran
PREPERROR with "Error: Bioconductor version '3.13' requires R version '4.1'; R version is too new". Seems to be an issue with rhub rather than with the package.

### Debian Linux, R-devel, GCC ASAN/UBSAN



## Downstream dependencies
This is a new package, so there are no reverse dependencies.
