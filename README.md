# R package to measure various phenomenological ecological portfolio effects.

A financial portfolio metaphor is often used to describe how population
diversity can make a group of populations more stable, much like how
investing in a diverse set of assets makes a financial portfolio more
stable. This package contains functions to estimate phenomenological
ecological portfolio effects. Currently, the functions focus on
single-species metapopulations.

This package accompanies the paper:  
Anderson, S.C., A.B. Cooper, N.K. Dulvy. 2013. Methods in Ecology and
Evolution. In Press.

You can install the package with:

```r
# Install these first if needed:
# install.packages(c("plyr", "reshape", "MuMIn", "robustbase", "devtools"))
devtools::install_github("ecofolio", username="seananderson")
```

If you're using Windows then you may need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) first.

You can load the package, read the vignette, and view the help pages
with:

```r
library(ecofolio)
vignette("ecofolio")
help(package = "ecofolio")
```
