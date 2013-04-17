\name{pe_mv}
\alias{pe_mv}
\title{Estimate the mean-variance portfolio effect}
\usage{
  pe_mv(x,
    fit_type = c("linear", "linear_robust", "quadratic", "linear_quad_avg", "linear_detrended", "loess_detrended"),
    ci = FALSE, boot = FALSE, boot_reps = 1000)
}
\arguments{
  \item{x}{A matrix of abundance or biomass data. The
  columns should represent different subpopulations or
  species. The rows should represent the values through
  time.}

  \item{fit_type}{Type of model to fit to the
  log(variance)-log(mean) data.  Options are: \itemize{
  \item \code{linear}: linear regression, \item
  \code{linear_robust}: robust linear regression \item
  \code{quadratic}: quadratic regression \item
  \code{linear_quad_avg}: AICc-weighted model averaging of
  linear and quadratic regression \item
  \code{linear_detrended}: detrend the time series with a
  linear model before estimating z \item
  \code{loess_detrended}: detrend the time series with a
  loess smoother before estimating z }}

  \item{ci}{Logical value describing whether a 95\%
  confidence interval should be calculated and returned
  (defaults to \code{TRUE}).}

  \item{boot}{Logical value (defaults to \code{FALSE}).
  Determines whether the confidence interval should be
  calculated using the (bias-adjusted) bootstrap instead of
  using the parametric confidence interval from the linear
  model fit.}

  \item{boot_reps}{Number of bootstrap repetitions.}
}
\value{
  A numeric value representing the portfolio effect that
  takes into account the mean-variance relationship. If
  confidence intervals were requested then a list is
  returned with the portfolio effect (pe) and 95\%
  confidence interval (ci).
}
\description{
  Takes a matrix of abundance or biomass data and returns
  various estimates of the mean-variance portfolio effect.
  Options exist to fit various mean-variance models and to
  detrend the time series data.
}
\details{
  This version of the portfolio effect consists of dividing
  the CV of a theoretical single population (single asset
  system) that has the same overall mean but with the
  variance scaled according to the mean-variance
  relationship by the CV of the combined total population.
  The calculation of the portfolio CV is the same as in
  \code{\link{pe_avg_cv}} but the calculation of the single
  asset system CV is different.
}
\examples{
dat = data.frame(x1 = rnorm(20, 10), x2 = rnorm(20, 10), x3 = rnorm(20,10))
pe_mv(dat)
}
\references{
  Doak, D., D. Bigger, E. Harding, M. Marvier, R. O'Malley,
  and D.  Thomson. 1998. The Statistical Inevitability of
  Stability-Diversity Relationships in Community Ecology.
  Amer. Nat. 151:264–276.

  Tilman, D., C. Lehman, and C. Bristow. 1998.
  Diversity-Stability Relationships: Statistical
  Inevitability or Ecological Consequence? Amer.  Nat.
  151:277–282.

  Tilman, D. 1999. The Ecological Consequences of Changes
  in Biodiversity: A Search for General Principles. Ecology
  80:1455–1474.

  Taylor, L. 1961. Aggregation, Variance and the Mean.
  Nature 189:732–735.  doi: 10.1038/189732a0.

  Taylor, L., I. Woiwod, and J. Perry. 1978. The
  Density-Dependence of Spatial Behaviour and the Rarity of
  Randomness. J. Anim. Ecol. 47:383–406.
}

