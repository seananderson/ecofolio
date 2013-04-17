# ====================================================================
# Created by:    Sean Anderson, sean@seananderson.ca
# Created:       Feb 20, 2012, taken from the simulation file from
# June
# Last modified: Nov 20, 2012
# Purpose:       the functions to calculate theoretical PEs
# ====================================================================

odd.vector <- function(x) {
# return the odd indicies of a vector 
    vec <- 1:length(x)
    odd.vec <- vec[vec != as.integer(vec/2) * 2]
    x[odd.vec]
}

# a = c for purposes of object names
PE.tilman.2.assets <- function(mu1, mu2, z, a = 1) {
  CV.p <- sqrt(a * mu1 ^ z + a * mu2 ^ z) / (mu1 + mu2)
  CV.a <- sqrt(a * (mu1 + mu2) ^ z)/(mu1 + mu2)
  CV.p/CV.a
}

PE.tilman.2.assets.cor <- function(mu1, mu2, z, r, a = 1) {
  CV.p <- sqrt(mu1^z + mu2^z + 2 * (r * sqrt(mu1^z * mu2^z))) / (mu1 + mu2)
  CV.a <- sqrt(a * (mu1 + mu2) ^ z)/(mu1 + mu2)
  CV.p/CV.a
}

PE.schindler.2.assets <- function(mu1, mu2, z, a = 1) {
  CV.p <- sqrt(a * mu1 ^ z + a * mu2 ^ z) / (mu1 + mu2)
  CV.a <- (sqrt(a * mu1 ^ z) / mu1 + sqrt(a * mu2 ^ z)/ mu2) / 2
  CV.p/CV.a
}

PE.schindler.2.assets.cor <- function(mu1, mu2, z, r, a = 1) {
  CV.p <- sqrt(mu1^z + mu2^z + 2 * (r * sqrt(mu1^z * mu2^z))) / (mu1 + mu2)
  CV.a <- (sqrt(a * mu1 ^ z) / mu1 + sqrt(a * mu2 ^ z)/ mu2) / 2
  CV.p/CV.a
}

create.cor.dat.2.assets <- function(mu1, mu2, r, a, z, n = 1000) {
  require(mvtnorm)
  var1 <- a * mu1 ^ z
  var2 <- a * mu2 ^ z
  E <- matrix(c(var1, r * sqrt(var1 * var2), r * sqrt(var1 * var2), var2), 2, 2) 
  d <- rmvnorm(n, mean = c(mu1, mu2), E) 
  data.frame(mu1 = d[,1], mu2 = d[,2])
}

combinations <- function(x, m) nrow(combn(x, m)) * ncol(combn(x, m))

# TODO should these be multiplied by 2!!!!????? - all combinations
# answer - no, apparently not
sum.variance.cor <- function(k, mu, z, r) {
## assumes all mu the same
  v.s <- rep(mu ^ z, k)
  mu.s <- rep(mu, k)
  #browser()
  #combined.v <- sum(v.s) + 2 * sum(r * sqrt(rep(mu ^ z, combinations(k, 2))^2))
  combined.v <- sum(v.s) + sum(r * sqrt(rep(mu ^ z, combinations(k, 2))^2))
  ifelse(combined.v >= 0, combined.v, NA)
}

sum.variance.cor.generic <- function(mu.s, z, r) {
# for any k, r, mu, z
  k <- length(mu.s)
  v.s <- mu.s ^ z
  cor.matrix <- matrix(data = r, ncol = k, nrow = k)
  diag(cor.matrix) <- 1
  cov.matrix <- sqrt(v.s) %*% t(sqrt(v.s)) * cor.matrix
### weighted mean, not using for now:
  #X <- rep(1, length(mu.s))
  #X.matrix <- X %*% t(X)
  #combined.v <- sum(X.matrix * cov.matrix)
###
  combined.v <- sum(cov.matrix)
  ifelse(combined.v >= 0, combined.v, NA)
}

get_k_prob <- function(k) seq(1/(k+1), 1-(1/(k+1)), length.out = k)

# sim: 
#library(mvtnorm)
#z <- 1.5
#mu1 <- 3
##mu2 <- 5
#mu2 <- 7
#r.temp <- 0.9
#n <- 100000
#j <- rmvnorm(n, mean = c(mu1, mu2), matrix(c(mu1^z, r.temp, r.temp, mu2^z), ncol = 2))
## actual correlation = 0.0174
#cor(j[,1], j[,2])
#mean(apply(j, 1, sum))
#var(apply(j, 1, sum))

PE.schindler.k.cor <- function(k, mu, z, r) {
  mu.s <- rep(mu, k)
  v.s <- rep(mu ^ z, k)
  combined.v <- sum.variance.cor(k, mu, z, r)
  CV.p <- sqrt(combined.v) / sum(rep(mu, k))
  CV.a <- mean(sqrt(v.s) / mu.s)
  CV.p / CV.a
}

PE.tilman.k.cor <- function(k, mu, z, r) {
  mu.s <- rep(mu, k)
  combined.v <- sum.variance.cor(k, mu, z, r)
  if(!is.na(combined.v)) {
  CV.p <- sqrt(combined.v) / sum(rep(mu, k))
  CV.a <- sqrt(sum(mu.s) ^ z) / sum(mu.s)
  CV.p / CV.a
  }
  else
    NA
}



