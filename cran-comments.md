---
title: "Cran Comments"
author: "Ivan Svetunkov"
date: "03 February 2025"
output: html_document
---
## Version
This is an initial submission of ``legion`` package, v0.2.1.


## Test environments
* local ubuntu 24.10, R 4.4.1
* github actions
* win-builder (devel and release)


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


## Downstream dependencies
There are no reverse dependencies.
