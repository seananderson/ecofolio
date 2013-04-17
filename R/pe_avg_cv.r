#' Estimate the average-CV portfolio effect
#' 
#' Takes a matrix of abundance or biomass data and returns various estimates of
#' the average-CV portfolio effect. Options exist to detrend the time series
#' data.
#' 
#' @references Doak, D., D. Bigger, E. Harding, M. Marvier, R. O'Malley, and D. 
#' Thomson. 1998. The Statistical Inevitability of Stability-Diversity 
#' Relationships in Community Ecology. Amer. Nat. 151:264–276.
#' 
#' Tilman, D., C. Lehman, and C. Bristow. 1998. Diversity-Stability 
#' Relationships: Statistical Inevitability or Ecological Consequence? Amer. 
#' Nat. 151:277–282.
#' 
#' Schindler, D., R. Hilborn, B. Chasco, C. Boatright, T. Quinn, L. Rogers, and
#' M. Webster. 2010. Population diversity and the portfolio effect in an
#' exploited species. Nature 465:609-612. doi: 10.1038/nature09060.
#' 
#' @details This version of the portfolio effect consists of dividing the mean
#'   of the coefficient of variations (CV) of all individual subpopulations
#'   (assets) by the CV of the combined total population.
#'   
#' @param x A matrix of abundance or biomass data. The columns should represent 
#'   different subpopulations or species. The rows should represent the values 
#'   through time.
#' @param fit_type Character value describing if (and how) the time series
#'   should be detrended before estimating the portfolio effect. Defaults to not
#'   detrending.
#' @param ci Logical value (defaults to \code{FALSE}). Should a 95\% confidence
#'   interval should be calculated using a bootstrap procedure? Returns the
#'   bias-corrected (bca) version of the bootstrap confidence interval.
#' @param boot_reps Number of bootstrap replicates.
#'   
#' @return A numeric value representing the average-CV portfolio effect. If 
#' confidence intervals were requested then a list is returned with the
#' portfolio effect \code{pe} and 95\% bootstrapped confidence interval 
#' \code{ci}.
#' 
#' @examples 
#' dat = data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
#' pe_avg_cv(dat, ci = TRUE)
#' 
#' @export

pe_avg_cv <- function (x, fit_type = c("not_detrended", "linear_detrended", 
  "loess_detrended"),  ci = FALSE, boot_reps = 1000
)
{

if(!fit_type[1] %in% c("not_detrended", "linear_detrended", "loess_detrended")) 
  stop("not a valid fit_type")

if(fit_type[1] == "not_detrended") {
  cv_single_asset <- mean(apply(x, 2, cv))
  cv_portfolio <- cv(rowSums(x))
  pe <- cv_portfolio / cv_single_asset
}

if(fit_type[1] == "linear_detrended") {
# single assets:
  x_detrended <- x
  for(i in 1:ncol(x)) x_detrended[,i] <- residuals(lm(x[,i]~c(1:nrow(x))))
  single_asset_means <- apply(x, 2, mean)
  single_asset_sds <- apply(x_detrended, 2, sd)
  cv_single_asset <- mean(single_asset_sds / single_asset_means)
# portfolio:
  sd_portfolio <- sd(residuals(lm(rowSums(x)~c(1:nrow(x)))))
  mean_portfolio <- mean(rowSums(x))
  cv_portfolio <- sd_portfolio / mean_portfolio
  pe <- cv_portfolio / cv_single_asset
}

if(fit_type[1] == "loess_detrended") {
# single assets:
  x_detrended <- x
  for(i in 1:ncol(x)) x_detrended[,i] <- residuals(loess(x[,i]~c(1:nrow(x))))
  single_asset_means <- apply(x, 2, mean)
  single_asset_sds <- apply(x_detrended, 2, sd)
  cv_single_asset <- mean(single_asset_sds / single_asset_means)
# portfolio:
  sd_portfolio <- sd(residuals(loess(rowSums(x)~c(1:nrow(x)))))
  mean_portfolio <- mean(rowSums(x))
  cv_portfolio <- sd_portfolio / mean_portfolio
  pe <- cv_portfolio / cv_single_asset
}

  pe_avg_cv_for_boot <- function(x) {
    cv_single_asset <- mean(apply(x, 2, cv))
    cv_portfolio <- cv(rowSums(x))
    pe <- cv_portfolio / cv_single_asset
    #pe
  }
  if(ci) {
## confidence interval calculation
    require(boot)
    boot.out <- boot(t(x), function(y, i) pe_avg_cv_for_boot(t(y[i,])), R = boot_reps)
    pe_ci <- boot.ci(boot.out, type = "bca")$bca[c(4,5)]
    list(pe = pe, ci = pe_ci)
  } else pe
}
