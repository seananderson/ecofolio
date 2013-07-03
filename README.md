`ecofolio` is an R package to estimate ecological portfolio effects.
Presently it focuses mainly on metapopulations.

A financial portfolio metaphor is often used to describe how population
diversity can make a group of populations more stable, much like how
investing in a diverse set of assets makes a financial portfolio more
stable. For example, we can treat the abundance of salmon in one stream
as financial stock value and salmon within a river catchment as
financial portfolios. 

If a group of populations are stressed and react differently then the
risk of collapse or decline of the group may be lowered, similar to
a diversified financial portfolio. This risk reduction has been referred
to as the "portfolio effect".

To cite this package please cite the accompanying paper:  
Anderson, S.C., A.B. Cooper, N.K. Dulvy. 2013. Ecological prophets:
Quantifying metapopulation portfolio effects. Methods in Ecology and
Evolution. In Press.

You can install `ecofolio` with the following R code:

```r
# Install these first if needed:
# install.packages(c("plyr", "reshape", "MuMIn", "robustbase", "devtools"))
devtools::install_github("ecofolio", username="seananderson")
```

If you're using Windows then you may need to install
[Rtools](http://cran.r-project.org/bin/windows/Rtools/) first.

You can load the package, read the vignette, and view the help pages
with:

```r
library(ecofolio)
vignette("ecofolio")
help(package = "ecofolio")
```

If you have problems installing the package or using the functions,
please email sean "at" seananderson.
