## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, fig.path='Figs/', fig.show='hold',
                      warning=FALSE, message=FALSE)
library(smooth)

## ----artificialData------------------------------------------------------
y <- ts(c(rpois(20,0.25),rpois(20,0.5),rpois(20,1),rpois(20,2),rpois(20,3),rpois(20,5)))

## ----iETSFExample1-------------------------------------------------------
oETSFModel1 <- oes(y, occurrence="fixed", h=10, holdout=TRUE)
oETSFModel1
plot(oETSFModel1)

## ----iETSFExample2-------------------------------------------------------
es(y, "MMN", occurrence="fixed", h=10, holdout=TRUE, silent=FALSE)

## ----iETSOExample1-------------------------------------------------------
oETSOModel <- oes(y, model="MMN", occurrence="o", h=10, holdout=TRUE)
oETSOModel
plot(oETSOModel)

## ----iETSOExample2-------------------------------------------------------
es(y, "MMN", occurrence="o", oesmodel="MMN", h=10, holdout=TRUE, silent=FALSE)

## ----iETSOExample3, eval=FALSE-------------------------------------------
#  es(y, "MMN", occurrence=oETSOModel, h=10, holdout=TRUE, silent=FALSE)

## ----iETSIExample1-------------------------------------------------------
oETSIModel <- oes(y, model="MMN", occurrence="i", h=10, holdout=TRUE)
oETSIModel
plot(oETSIModel)

## ----iETSIExample2-------------------------------------------------------
es(y, "MMN", occurrence="i", oesmodel="MMN", h=10, holdout=TRUE, silent=FALSE)

## ----iETSIExample3, eval=FALSE-------------------------------------------
#  es(y, "MMN", occurrence=oETSIModel, h=10, holdout=TRUE, silent=FALSE)

## ----iETSDExample1-------------------------------------------------------
oETSDModel <- oes(y, model="MMN", occurrence="d", h=10, holdout=TRUE)
oETSDModel
plot(oETSDModel)

## ----iETSDExample2-------------------------------------------------------
es(y, "MMN", occurrence=oETSDModel, h=10, holdout=TRUE, silent=FALSE)

## ----iETSGExample1-------------------------------------------------------
oETSGModel1 <- oesg(y, modelA="MNN", modelB="AAN", h=10, holdout=TRUE)
oETSGModel1
plot(oETSGModel1)

## ----iETSGExample2-------------------------------------------------------
oETSGModel2 <- oes(y, model="MNN", occurrence="g", h=10, holdout=TRUE)
oETSGModel2
plot(oETSGModel2)

## ----iETSGExample3-------------------------------------------------------
es(y, "MMN", occurrence="g", oesmodel="MMN", h=10, holdout=TRUE, silent=FALSE)

## ----iETSAExample1-------------------------------------------------------
oETSAModel <- oes(y, model="MNN", occurrence="a", h=10, holdout=TRUE)
oETSAModel
plot(oETSAModel)

## ----iETSGRoundedExample-------------------------------------------------
es(rpois(100,0.3), "MNN", occurrence="g", oesmodel="MNN", h=10, holdout=TRUE, silent=FALSE, interval=TRUE, rounded=TRUE)

