
#' @title Bias of the sample proportion
#' @description Calculates bias due to using sample proportion as estimator of the true
#' response rate.
#' @details For fixed designs the sample propotion is an unbiased (maximum likelihood)
#' estimator of the response rate, but in group sequential designs (e.g., Simon's)
#' it is biased.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @param p True success probability.
#' @return Bias.
#' @references Porcher, R. and Desseaux, K. What inference for two-stage phase II trials?
#' \emph{BMC Med Res Methodol},  2012, 12, 117.
#' @seealso \code{\link{sfms}} and \code{\link{pg}}.
#' @export
#' @examples
#' sbias(19, 4, 54, 0.4)
#' @author Arsenio Nhacolo
sbias <- function(n1, r1, n, p){
  sum1 <- 0
  sum2 <- 0
  for (x in 0:r1){
    sum1 <- sum1+x*sfms(x, n1, r1, n, p, m = 1)
  }
  for (x in (r1+1):n){
    sum2 <- sum2+x*sfms(x, n1, r1, n, p, m = 2)
  }
  return(sum1/n1+sum2/n-p)
}
