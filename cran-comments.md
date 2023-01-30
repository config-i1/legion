---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "30 January 2023"
output: html_document
---
## Version
This is an initial submission of ``legion`` package, v0.1.2.


## Test environments
* local ubuntu 22.04.1, R 4.2.2
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
>* checking package dependencies ... ERROR
>Package suggested but not available: 'doMC'

doMC is not available for Windows.


## rhub checks
### Fedora Linux, R-devel, clang, gfortran
>Found the following (possibly) invalid DOIs:
>  DOI: 10.1177/1471082X0901000401
>    From: DESCRIPTION
>    Status: Service Unavailable
>    Message: 503

DOI is correct, the paper is available online.

### Windows Server 2022, R-devel, 32/64 bit
>* checking package dependencies ... ERROR
>Package suggested but not available: 'doMC'

doMC is not available for Windows.


## Downstream dependencies
This is a new package, so there are no reverse dependencies.
