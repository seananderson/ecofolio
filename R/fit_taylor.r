#' Fit Taylor's power law
#'
#' Fits Taylor's power law to the temporal mean and variance in log-log space
#' and returns the coefficients.
#'
#' The function returns: A list containing the constant \code{c} value and the
#' exponent \code{z} in Taylor's power law equation. If confidence intervals
#' were requested then the list will also contain \code{ci} with the 95\%
#' confidence intervals on the z value.
#'
#' @param x A matrix or dataframe of abundance or biomass data. The columns
#' should represent different subpopulations or species. The rows should
#' represent the values through time.
#' @param ci A Logical value indicating whether 95\% confidence intervals should
#' be calculated for the z value (the exponent in Taylor's power law).
#' @param na.rm A logical value indicating whether \code{NA} values should be
#' row-wise deleted. 
#'
#' @references
#' Taylor, L. 1961. Aggregation, Variance and the Mean. Nature 189:732-735. doi:
#' 10.1038/189732a0.
#'
#' Taylor, L., I. Woiwod, and J. Perry. 1978. The Density-Dependence of Spatial
#' Behaviour and the Rarity of Randomness. J. Anim. Ecol. 47:383-406.
#'
#' Taylor, L., and I. Woiwod. 1982. Comparative Synoptic Dynamics. I.
#' Relationships Between Inter- and Intra-Specific Spatial and Temporal
#' Variance/Mean Population Parameters. J. Anim. Ecol. 51:879-906.
#' @export
#' @examples
#' data(pinkbr)
#' fit_taylor(pinkbr[,-1])

fit_taylor <- function(x, ci = FALSE, na.rm = FALSE){
  if(na.rm) x <- na.omit(x)
  
  total_nas <- sum(is.na(x))
  ifelse(!na.rm & total_nas > 0, return_na <- TRUE, return_na <- FALSE)
  
  m <- apply(x, 2, mean)
  v <- apply(x, 2, var)
  log.m <- log(m)
  log.v <- log(v)
  fit <- lm(log.v ~ log.m)
  c.value <- as.numeric(coef(fit)[1])
  z.value <- as.numeric(coef(fit)[2])
  if(ci == TRUE) {
    z.se <- summary(fit)$coef[2,2]
    z.l <- z.value - 1.96 * z.se
    z.u <- z.value + 1.96 * z.se
    out <- list(c = c.value, z = z.value, z.l = z.l, z.u = z.u)
  }else{
    out <-  list(c = c.value, z = z.value) 
  }
  if(return_na){
    out[[1]] <- NA
    out[[2]] <- NA
  }
  return(out)
}

