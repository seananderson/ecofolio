ecofolio is an R package to estimate ecological portfolio effects. Currently it
focuses on metapopulations or metapopulation-like systems.

A financial portfolio metaphor is often used to describe how population
diversity can make a group of populations more stable, much like how investing
in a diverse set of assets makes a financial portfolio more stable. For
example, we can treat the abundance of salmon in one stream as financial stock
value and the abundance of salmon within a river catchment as a financial
portfolio. If a group of populations are stressed and react differently then
the risk of collapse or decline of the group may be lowered, similar to
a diversified financial portfolio. This risk reduction has been referred to as
the "portfolio effect".

To cite this package please cite the accompanying paper:  

Anderson, S.C., A.B. Cooper, N.K. Dulvy. 2013. Ecological prophets:
Quantifying metapopulation portfolio effects. Methods in Ecology and
Evolution. In Press. doi: 10.1111/2041-210X.12093.

The published version of the paper is available [here](http://dx.doi.org/10.1111/2041-210X.12093). And a pre-print is available [here](http://seananderson.ca/papers/Anderson_etal_2013_ecological_prophets.pdf) [PDF].

You can install ecofolio with the following:

```r
# Install these first if needed:
pak::pak("seananderson/ecofolio")
```

You can load the package, read the vignette, and view the help pages with:

```r
library(ecofolio)
vignette("ecofolio")
help(package = "ecofolio")
```

If you have problems installing the package or using the functions, please
email sean "at" seananderson.ca.
