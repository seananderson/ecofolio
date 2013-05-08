#' Plot mean-variance relationship
#'
#' Creates a scatter plot of the time series log(variance) vs.
#' log(mean). Shows a model fit to the mean-variance data and an
#' extrapolation to the size of the metapopulation. The linear version
#' of this model is referred to as Taylor's power law.
#'
#' @param x A matrix or dataframe of abundance or biomass data. The
#' columns should represent different subpopulations or species. The
#' rows should represent the values through time.
#' @param show A vector of character objects indicating which
#' mean-variance models to show.
#' @param col Colour for the mean-variance model fit. A vector of
#' length 3 with the three values corresponding to \code{linear},
#' \code{quadratic},
#' \code{robust}.                                                                                                                          
#' @param lty Line type for the mean-variance model fit. A vector of
#' length 3 with the three values corresponding to \code{linear},
#' \code{quadratic},
#' \code{robust}.                                                                                                                          
#' @param pch_sa Point type for the extrapolated
#' "single-asset" portfolio. A vector of length 3 with the three
#' values corresponding to \code{linear}, \code{quadratic},
#' \code{robust}.
#' @param ci Add a confidence interval around the model fit? Only
#' appears for the linear fit option.
#' @param pch_subpops Point type for the subpopulations.
#' @param pch_port Point type for the portfolio.
#' @param add_z Logical. Add Taylor's power law z value (based on a
#' linear model)?
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... Other values to pass to \code{plot}.
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
#' Taylor, L. 1961. Aggregation, Variance and the Mean. Nature 189:732-735. doi:
#' 10.1038/189732a0.
#' 
#' Taylor, L., I. Woiwod, and J. Perry. 1978. The Density-Dependence of Spatial
#' Behaviour and the Rarity of Randomness. J. Anim. Ecol. 47:383-406.
#' @export
#' @examples
#' data(pinkbr)
#' par(mfrow = c(1,3))
#' plot_mv(pinkbr[,-1], show = "linear")
#' mtext("Linear")
#' plot_mv(pinkbr[,-1], show = "quadratic", add_z = FALSE)
#' mtext("Quadratic")
#' plot_mv(pinkbr[,-1], show = "robust", add_z = FALSE)
#' mtext("Robust linear")

plot_mv <- function(x, show = c("linear", "quadratic", "robust"), col
  = c("#D95F02", "#1B9E77", "#E7298A"), lty = c(1, 1, 1),
  pch_sa = c(1, 5, 6), ci = FALSE, pch_subpops = 21,
  pch_port = 4, add_z = TRUE, xlab = "log(mean)", ylab =
  "log(variance)", ...) {

  require(plyr)
  require(reshape)

  ## get mean and variance of portfolio and assets:
  x.long <- melt(x, id.vars = NULL)
  overall.d <- apply(x[,-1], 1, sum)
  mv <- ddply(x.long, "variable", summarize, m = mean(value, na.rm =
      TRUE), v = var(value, na.rm = TRUE))
  overall.mean <- mean(overall.d, na.rm = TRUE)
  portfolio.var <- var(overall.d, na.rm = TRUE)
  m.t <- lm(log(v) ~ log(m), data = mv)
  overall.variance <- exp(as.numeric(predict(m.t, newdata =
          data.frame(m = overall.mean))))
  m.t.quad <- nls(log(v) ~ B0 + B1 * log(m) + B2 * I(log(m) ^ 2), data
      = mv, start = list(B0 = 0, B1 = 2, B2 = 0), lower = list(B0 =
        -1e9, B1 = 0, B2 = 0), algorithm = "port")
  overall.variance.quad <- exp(as.numeric(predict(m.t.quad, newdata =
          data.frame(m = overall.mean))))
    d1p <- seq(min(mv$m), max(mv$m), length.out = 200)
    d2p <- seq(max(mv$m), overall.mean, length.out = 200)

  ## set up plot:
  with(mv, plot((m), (v), xlim = c((min(m)), (overall.mean)*1.15),
      ylim = c(min((v)), max(overall.variance, portfolio.var,
          overall.variance.quad)*1.25), pch = pch_subpops, log = "xy", xlab =
      xlab, ylab =ylab, bg = "#00000020", col = "#00000070", ...))
  points((overall.mean), (portfolio.var), col = "black", pch = pch_port, lwd
    = 1.2, cex = 1.6)
  ## now add fits and extrapolations as requested:
  if("linear" %in% show) {
    m.t.p <- predict(m.t, newdata = data.frame(m = seq(min(mv$m),
          max(mv$m), length.out = 2)))
    m.t.p2 <- predict(m.t, newdata = data.frame(m = seq(max(mv$m),
          overall.mean, length.out = 2)))
    p1 <- predict(m.t, newdata = data.frame(m = d1p), se = TRUE)
    p2 <- predict(m.t, newdata = data.frame(m = d2p), se = TRUE)
    points((overall.mean), (overall.variance), col = col[1], pch = pch_sa[1],
      lwd = 1.5, cex = 1.4)
    segments((min(mv$m)), exp(m.t.p[1]), (max(mv$m)), exp(m.t.p[2]),
      col = col[1], lty = lty[1])
    segments((max(mv$m)), exp(m.t.p2[1]), (overall.mean),
      exp(m.t.p2[2]), lty = 2,  col = col[1])
    if(ci) {
      polygon(c(d1p, rev(d1p)), c(exp(p1$fit + 1.96*p1$se.fit),
          exp(rev(p1$fit - 1.96*p1$se.fit))), border = FALSE, col =
        "#00000020")
      polygon(c(d2p, rev(d2p)), c(exp(p2$fit + 1.96*p2$se.fit),
          exp(rev(p2$fit - 1.96*p2$se.fit))), border = FALSE, col =
        "#00000010")
    }
  }
  if("quadratic" %in% show){
    p.quad.1 <- predict(m.t.quad, newdata = data.frame(m = d1p), se = FALSE)
    p.quad.2 <- predict(m.t.quad, newdata = data.frame(m = d2p), se = FALSE)
    points((overall.mean), (overall.variance.quad), col =col[2], pch =
      pch_sa[2], lwd = 1.5, cex = 1.4)
    lines(d1p, exp(as.numeric(p.quad.1)), col = col[2], lwd = lty[2])
    lines(d2p, exp(as.numeric(p.quad.2)), col = col[2], lty = 2, lwd = 1.5)
  }
  if("robust" %in% show) {
    require(robustbase)
    m.t.rob <- lmrob(log(v) ~ log(m), data = mv)
    overall.variance.rob <- exp(as.numeric(predict(m.t.rob, newdata =
          data.frame(m = overall.mean))))
    p.rob.1 <- predict(m.t.rob, newdata = data.frame(m = d1p), se = FALSE)
    p.rob.2 <- predict(m.t.rob, newdata = data.frame(m = d2p), se = FALSE)
    points((overall.mean), (overall.variance.rob), col = col[3], pch =
      pch_sa[3], lwd = 1.5, cex = 1.1)
    lines(d1p, exp(p.rob.1), col = col[3], lty = lty[3])
    lines(d2p, exp(p.rob.2), col = col[3], lty = 2)
  }
  if(add_z) {
    mtext(paste("z =", formatC(round(coef(m.t)[2], 2), digits = 1,
          format = "f")), side = 1, adj = 0.9, line = -1.2)
  }
}
