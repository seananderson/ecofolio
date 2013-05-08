\name{plot_mv}
\alias{plot_mv}
\title{Plot mean-variance relationship}
\usage{
  plot_mv(x, show = c("linear", "quadratic", "robust"),
    col = c("#D95F02", "#1B9E77", "#E7298A"),
    lty = c(1, 1, 1), pch_single_asset = c(1, 5, 6),
    ci = FALSE, pch_subpops = 21, pch_port = 4,
    add_z = TRUE, xlab = "log(mean)",
    ylab = "log(variance)", ...)
}
\arguments{
  \item{x}{A matrix or dataframe of abundance or biomass
  data. The columns should represent different
  subpopulations or species. The rows should represent the
  values through time.}

  \item{show}{A vector of character objects indicating
  which mean-variance models to show.}

  \item{col}{Colour for the mean-variance model fit. A
  vector of length 3 with the three values corresponding to
  \code{linear}, \code{quadratic}, \code{robust}.}

  \item{lty}{Line type for the mean-variance model fit. A
  vector of length 3 with the three values corresponding to
  \code{linear}, \code{quadratic}, \code{robust}.}

  \item{pch_single_asset}{Point type for the extrapolated
  "single-asset" portfolio. A vector of length 3 with the
  three values corresponding to \code{linear},
  \code{quadratic}, \code{robust}.}

  \item{ci}{Add a confidence interval around the model fit?
  Only appears for the linear fit option.}

  \item{pch_subpops}{Point type for the subpopulations.}

  \item{pch_port}{Point type for the portfolio.}

  \item{add_z}{Logical. Taylor's power law z value (based
  on a linear model?).}

  \item{xlab}{Label for x-axis.}

  \item{ylab}{Label for y-axis.}

  \item{...}{Other values to pass to \code{plot}.}
}
\description{
  Creates a scatter plot of the time series log(variance)
  vs. log(mean). Shows a model fit to the mean-variance
  data and an extrapolation to the size of the
  metapopulation. The linear version of this model is
  referred to as Taylor's power law.
}
\examples{
data(pinkbr)
par(mfrow = c(1,3))
plot_mv(pinkbr[,-1], show = "linear")
mtext("Linear")
plot_mv(pinkbr[,-1], show = "quadratic", add_z = FALSE)
mtext("Quadratic")
plot_mv(pinkbr[,-1], show = "robust", add_z = FALSE)
mtext("Robust linear")
}
\references{
  Doak, D., D. Bigger, E. Harding, M. Marvier, R. O'Malley,
  and D. Thomson. 1998. The Statistical Inevitability of
  Stability-Diversity Relationships in Community Ecology.
  Amer. Nat. 151:264-276.

  Tilman, D., C. Lehman, and C. Bristow. 1998.
  Diversity-Stability Relationships: Statistical
  Inevitability or Ecological Consequence? Amer. Nat.
  151:277-282.

  Tilman, D. 1999. The Ecological Consequences of Changes
  in Biodiversity: A Search for General Principles. Ecology
  80:1455-1474.

  Taylor, L. 1961. Aggregation, Variance and the Mean.
  Nature 189:732-735. doi: 10.1038/189732a0.

  Taylor, L., I. Woiwod, and J. Perry. 1978. The
  Density-Dependence of Spatial Behaviour and the Rarity of
  Randomness. J. Anim. Ecol. 47:383-406.
}

