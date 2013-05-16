#' Calculate the Loreau and de Mazencourt synchrony index
#'
#' Takes a matrix of abundance or biomass data and returns the
#' synchrony index. The synchrony index varies between 0 (maximally
#' asynchronous) and 1 (maximally synchronous) and is independent of
#' the number of species or subpopulations.
#'
#' @param x A matrix or dataframe of abundance or biomass data. The
#' columns should represent different subpopulations or species. The
#' rows should represent the values through time.
#'
#' @references
#' Loreau, M. & de Mazancourt, C. (2008). Species synchrony and its
#' drivers: neutral and nonneutral community dynamics in fluctuating
#' environments. Amer. Nat., 172, E48-66.
#'
#' Loreau, M. (2010). From Populations to Ecosystems: Theoretical
#' Foundations for a New Ecological Synthesis. Princeton University
#' Press, Princeton, NJ.
#' 
#' Thibaut, L.M. & Connolly, S.R. (2013). Understanding
#' diversity-stability relationships: towards a unified model of
#' portfolio effects. Ecology Letters, 16, 140-150.
#' @export

synchrony <- function (x) {
  var(rowSums(x)) / (sum(apply(x, 2, sd)) ^ 2)
}

