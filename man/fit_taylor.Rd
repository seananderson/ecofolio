\name{fit_taylor}
\alias{fit_taylor}
\title{Fit Taylor's power law}
\usage{
  fit_taylor(x, ci = FALSE, na.rm = FALSE)
}
\arguments{
  \item{x}{A matrix or dataframe of abundance or biomass
  data. The columns should represent different
  subpopulations or species. The rows should represent the
  values through time.}

  \item{ci}{A Logical value indicating whether 95\%
  confidence intervals should be calculated for the z value
  (the exponent in Taylor's power law).}

  \item{na.rm}{A logical value indicating whether \code{NA}
  values should be row-wise deleted.}
}
\description{
  Fits Taylor's power law to the temporal mean and variance
  in log-log space and returns the coefficients.
}
\details{
  The function returns: A list containing the constant
  \code{c} value and the exponent \code{z} in Taylor's
  power law equation. If confidence intervals were
  requested then the list will also contain \code{ci} with
  the 95\% confidence intervals on the z value.
}
\examples{
dat = data.frame(x1 = rlnorm(20), x2 = rlnorm(20), x3 = rlnorm(20))
fit_taylor(dat)
}
\references{
  Taylor, L. 1961. Aggregation, Variance and the Mean.
  Nature 189:732–735. doi: 10.1038/189732a0.

  Taylor, L., I. Woiwod, and J. Perry. 1978. The
  Density-Dependence of Spatial Behaviour and the Rarity of
  Randomness. J. Anim. Ecol. 47:383–406.

  Taylor, L., and I. Woiwod. 1982. Comparative Synoptic
  Dynamics. I. Relationships Between Inter- and
  Intra-Specific Spatial and Temporal Variance/Mean
  Population Parameters. J. Anim. Ecol. 51:879–906.
}

