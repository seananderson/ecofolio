\name{thibaut_connolly_dsr}
\alias{thibaut_connolly_dsr}
\title{Calculate the Thibaut and Connolly diversity-stability relationship}
\usage{
  thibaut_connolly_dsr(x, synchrony, z, overyielding)
}
\arguments{
  \item{x}{A matrix or dataframe of abundance or biomass
  data. The columns should represent different
  subpopulations or species. The rows should represent the
  values through time.}

  \item{synchrony}{The Loreau and de Mazencourt synchrony
  index. See \code{\link{synchrony}}.}

  \item{z}{Taylor's power law exponent from variance = c *
  mean^z. See \code{\link{fit_taylor}}.}

  \item{overyielding}{The overyielding coefficient. This
  reflects the increase in abundance or biomass with
  increasing diversity and is relevant to community
  portfolio effects.}
}
\value{
  A list containing the CV of the observed community
  (portfolio) \code{cv_p} expected CV of the community in
  monoculture \code{cv_1} and the ratio of the monoculture
  CV to the observed community or portfolio CV \code{pe}.
}
\description{
  Calculate the Thibaut and Connolly diversity-stability
  relationship
}
\examples{
dat = data.frame(x1 = rnorm(20, 10), x2 = rnorm(20, 10), x3 = rnorm(20,10))
thibaut_connolly_dsr(dat, synchrony = 0.7, z = 2, overyielding = 1)
}
\references{
  Thibaut, L.M. & Connolly, S.R. (2013). Understanding
  diversity-stability relationships: towards a unified
  model of portfolio effects. Ecology Letters, 16, 140-150.
}

