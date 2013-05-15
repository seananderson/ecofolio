\name{pe_avg_cv}
\alias{pe_avg_cv}
\title{Estimate the average-CV portfolio effect}
\usage{
  pe_avg_cv(x,
    detrending = c("not_detrended", "linear_detrended", "loess_detrended"),
    ci = FALSE, boot_reps = 500, na.rm = FALSE)
}
\arguments{
  \item{x}{A matrix or dataframe of abundance or biomass
  data. The columns should represent different
  subpopulations or species. The rows should represent the
  values through time.}

  \item{detrending}{Character value describing if (and how)
  the time series should be detrended before estimating the
  portfolio effect. Defaults to not detrending.}

  \item{ci}{Logical value (defaults to \code{FALSE}).
  Should a 95\% confidence interval be calculated using a
  bootstrap procedure? Returns the bias-corrected (bca)
  version of the bootstrap confidence interval.}

  \item{boot_reps}{Number of bootstrap replicates.}

  \item{na.rm}{A logical value indicating whether \code{NA}
  values should be row-wise deleted.}
}
\value{
  A numeric value representing the average-CV portfolio
  effect. If confidence intervals were requested then a
  list is returned with the portfolio effect \code{pe} and
  95\% bootstrapped confidence interval \code{ci}.
}
\description{
  Takes a matrix of abundance or biomass data and returns
  various estimates of the average-CV portfolio effect.
  Options exist to detrend the time series data.
}
\details{
  This version of the portfolio effect consists of dividing
  the mean of the coefficient of variations (CV) of all
  individual subpopulations (assets) by the CV of the
  combined total population.
}
\examples{
data(pinkbr)
pe_avg_cv(pinkbr[,-1], ci = TRUE)
pe_avg_cv(pinkbr[,-1], detrending = "loess_detrended", ci = TRUE)
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

  Schindler, D., R. Hilborn, B. Chasco, C. Boatright, T.
  Quinn, L. Rogers, and M. Webster. 2010. Population
  diversity and the portfolio effect in an exploited
  species. Nature 465:609-612. doi: 10.1038/nature09060.
}

