#' Estimate the mean-variance portfolio effect
#' 
#' Takes a matrix of abundance or biomass data and returns various estimates of
#' the mean-variance portfolio effect. Options exist to fit various
#' mean-variance models and to detrend the time series data.
#' 
#' @details This version of the portfolio effect consists of dividing the CV of
#' a theoretical single population (single asset system) that has the same
#' overall mean but with the variance scaled according to the mean-variance
#' relationship by the CV of the combined total population. The calculation of
#' the portfolio CV is the same as in \code{\link{pe_avg_cv}} but the
#' calculation of the single asset system CV is different.
#'   
#' @param x A matrix or dataframe of abundance or biomass data. The columns
#' should represent different subpopulations or species. The rows should
#' represent the values through time.
#' @param type Type of model to fit to the log(variance)-log(mean) data.
#' Options are: \itemize{ 
#' \item \code{linear}: linear regression (the default), 
#' \item \code{linear_robust}: robust linear regression 
#' \item \code{quadratic}: quadratic regression 
#' \item \code{linear_quad_avg}: AICc-weighted model averaging of linear and
#' quadratic regression 
#' \item \code{linear_detrended}: detrend the time series with a linear model
#' before estimating z from a linear regression
#' \item \code{loess_detrended}: detrend the time series with a loess smoother
#' before estimating z from a linear regression
#' }
#' @param ci Logical value describing whether a 95\% confidence interval should
#' be calculated and returned (defaults to \code{TRUE}).
#' @param na.rm A logical value indicating whether \code{NA} values should be
#' row-wise deleted. 
#'
#' @return A numeric value representing the portfolio effect that
#' takes into account the mean-variance relationship. If confidence
#' intervals were requested then a list is returned with the portfolio
#' effect (\code{pe}) and 95\% confidence interval (\code{ci}).
#'   
#' @references 
#' Doak, D., D. Bigger, E. Harding, M. Marvier, R. O'Malley, and D. Thomson.
#' 1998. The Statistical Inevitability of Stability-Diversity Relationships in
#' Community Ecology. Amer. Nat. 151:264-276.
#' 
#' Tilman, D., C. Lehman, and C. Bristow. 1998. Diversity-Stability
#' Relationships: Statistical Inevitability or Ecological Consequence? Amer.
#' Nat. 151:277-282.
#' 
#' Tilman, D. 1999. The Ecological Consequences of Changes in Biodiversity: A
#' Search for General Principles. Ecology 80:1455-1474.
#' 
#' Taylor, L. 1961. Aggregation, Variance and the Mean. Nature
#' 189:732-735. doi: 10.1038/189732a0.
#' 
#' Taylor, L., I. Woiwod, and J. Perry. 1978. The Density-Dependence of Spatial
#' Behaviour and the Rarity of Randomness. J. Anim. Ecol. 47:383-406.
#' @export
#' @examples
#' data(pinkbr)
#' pe_mv(pinkbr[,-1], ci = TRUE)
#' pe_mv(pinkbr[,-1], type = "quadratic") # same as linear in this case
#' pe_mv(pinkbr[,-1], type = "linear_quad_avg")
#' pe_mv(pinkbr[,-1], type = "linear_robust")
#' pe_mv(pinkbr[,-1], type = "linear_detrended", ci = TRUE)
#' pe_mv(pinkbr[,-1], type = "loess_detrended", ci = TRUE)

# TODO move predictions within model type to deal with robustbase name space

pe_mv <- function(x, type = c("linear", "linear_robust", "quadratic",
  "linear_quad_avg",  "linear_detrended", "loess_detrended"), ci =
    FALSE, na.rm = FALSE) {
  
  type <- type[1]
  
  if(!type %in% c("linear", "linear_robust", "quadratic",
      "linear_quad_avg", "linear_detrended", "loess_detrended")){
    stop("not a valid type")
  }
  
  if(!type %in% c("linear", "linear_robust", "linear_detrended",
      "loess_detrended")){
    if(ci == TRUE){
      warning("Confidence intervals aren't supported for this type of
        mean-variance model. Setting ci = FALSE.")
    }
    ci <- FALSE
  }
  
  if(na.rm) x <- na.omit(x)
  
  total_nas <- sum(is.na(x))
  return_na <- ifelse(!na.rm & total_nas > 0, TRUE, FALSE)
  
  ## first get the means:
  m <- apply(x, 2, mean)
  single_asset_mean <- mean(rowSums(x))

  cv_portfolio <- cv(rowSums(x))
  
  ## now detrend if desired:
  if(type == "linear_detrended") {
    ## first get cv of detrended portfolio abundance:
    sd_portfolio <- sd(residuals(lm(rowSums(x)~c(1:nrow(x)))))
    mean_portfolio <- mean(rowSums(x))
    cv_portfolio <- sd_portfolio / mean_portfolio
    ## now detrend:
    x <- apply(x, 2, function(y) residuals(lm(y~c(1:length(y)))))
  }
  if(type == "loess_detrended") {
    ## first get CV of detrended portfolio abundance:
    sd_portfolio <- sd(residuals(loess(rowSums(x)~c(1:nrow(x)))))
    mean_portfolio <- mean(rowSums(x))
    cv_portfolio <- sd_portfolio / mean_portfolio
    ## now detrend:
    x <- apply(x, 2, function(y) residuals(loess(y~c(1:length(y)))))
  }
  
  ## and get the variances for the assets:
  v <- apply(x, 2, var)

  log.m <- log(m)
  log.v <- log(v)
  d <- data.frame(log.m = log.m, log.v = log.v, m = m, v = v)
  
  taylor_fit <- switch(type[1], 
    linear = {
      lm(log.v ~ log.m, data = d)
    },
    linear_robust = {
      robustbase::lmrob(log.v ~ log.m, data = d)
    },
    quadratic = {
      nls(log.v ~ B0 + B1 * log.m + B2 * I(log.m ^ 2),
        data = d, start = list(B0 = 0, B1 = 2, B2 = 0), lower =
        list(B0 = -1e9, B1 = 0, B2 = 0), algorithm = "port")
    },
    linear_detrended = {
      lm(log.v ~ log.m, data = d)
    }, 
    loess_detrended = {
      lm(log.v ~ log.m, data = d)
    },
    linear_quad_avg = {
      linear <- nls(log.v ~ B0 + B1 * log.m, data = d, start = list(B0 =
          0, B1 = 2), lower = list(B0 = -1e9, B1 = 0), algorithm =
          "port")
      quadratic <- nls(log.v ~ B0 + B1 * log.m + B2 * I(log.m ^ 2), data
        = d, start = list(B0 = 0, B1 = 2, B2 = 0), lower = list(B0 =
            -1e9, B1 = 0, B2 = 0), algorithm = "port")
      MuMIn::model.avg(list(linear=linear, quad=quadratic), rank = MuMIn::AICc)
    }
  )

  if(ci) {
    single_asset_variance_predict <- predict(taylor_fit, newdata =
      data.frame(log.m = log(single_asset_mean)), se = TRUE)
    single_asset_variance <- exp(single_asset_variance_predict$fit)
  } else {
    single_asset_variance_predict <- predict(taylor_fit, newdata =
      data.frame(log.m = log(single_asset_mean)), se = FALSE)
    single_asset_variance <- exp(single_asset_variance_predict)
  }
  
  cv_single_asset <- sqrt(single_asset_variance) / single_asset_mean
  pe <- as.numeric(cv_single_asset / cv_portfolio)
  
  if(ci == TRUE) {
    single_asset_variance_ci <- exp(single_asset_variance_predict$fit
      + c(-1.96, 1.96) * single_asset_variance_predict$se.fit)
    cv_single_asset_ci <- sqrt(single_asset_variance_ci) / single_asset_mean
    pe_ci <- as.numeric(cv_single_asset_ci / cv_portfolio)
    pe_ci <- pe_ci[order(pe_ci)] # make sure the lower value is first
    out <- list(pe = pe, ci = pe_ci)
  } else {
    out <- pe
  }

  if(return_na) out <- NA

  out
}


