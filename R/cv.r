#' Calculate the coefficient of variation
#'
#' @param x A numeric vector
#'
#' @examples 
#' cv(rnorm(20))
#'
#' @export

cv <- function(x){
  sd(x)/abs(mean(x))
}
