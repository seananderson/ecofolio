#' Calculate the Thibaut and Connolly diversity-stability relationship
#' 
#' @param x A matrix or dataframe of abundance or biomass data. The
#'   columns should represent different subpopulations or species. The
#'   rows should represent the values through time.
#' @param synchrony The Loreau and de Mazencourt synchrony index. See
#' \code{\link{synchrony}}.
#' @param z Taylor's power law exponent from variance = c * mean^z.
#' See \code{\link{fit_taylor}}.
#' @param overyielding The overyielding coefficient. This reflects the
#' increase in abundance or biomass with increasing diversity and is
#' relevant to community portfolio effects.
#' @return A list containing the CV of the observed community
#' (portfolio) \code{cv_p} expected CV of the community in monoculture
#' \code{cv_1} and the ratio of the monoculture CV to the observed
#' community or portfolio CV \code{pe}.
#' @examples
#' dat = data.frame(x1 = rnorm(20, 10), x2 = rnorm(20, 10), x3 = rnorm(20,10))
#' thibaut_connolly_dsr(dat, synchrony = 0.7, z = 2, overyielding = 1)
#' @export

thibaut_connolly_dsr <- function(x, synchrony, z, overyielding) {
  cv_p <- cv(rowSums(x))
  n <- ncol(x)
  pe_inv <- sqrt(synchrony) * sqrt(n^((2 - z) * overyielding))
  cv1 <- cv_p / pe_inv
  list(cv_p = cv_p, cv_1 = cv1, pe = 1/pe_inv)
}

