# Flexible models for causal mediation analysis:<br/> an introduction to the R package medflex
### Pre-meeting workshop at UK Causal Inference Meeting 2017

Dear participant,

welcome to the github page of pre-meeting workshop B at the UK-CIM!
On this page, you can find all updated material, code and instructions for the workshop.

Slides and demo R code have just been uploaded.

Please bring your own laptop to this workshop if you wish to follow along with the R demo.
To have your laptop all set, please make sure that you have 

* the latest R and Rstudio release installed locally

* the latest development release of the medflex package installed (as it happens the latest CRAN release contains a bug that is fixed in the development release)
```R
devtools::install_github('jmpsteen/medflex')
```
This requires that you also have the `devtools` R package installed.
```R
install.packages('devtools')
```

* other recommended R packages installed: `mediation`, `car`, `multcomp`
 
If you already wish to obtain more background and detailed information on the use of the `medflex` package, please consult the companion paper in Journal of Statistical Software: https://www.jstatsoft.org/article/view/v076i11