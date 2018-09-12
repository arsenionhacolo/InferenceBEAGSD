#' @title Bias-reduced  estimator
#' @description Calculates the bias-reduced estimator of the true
#' response rate as proposed by \emph{Guo and Liu (2005)}.
#' @details It uses bias subtraction, with bias calculated by \code{\link{sbias}}
#' and response rate estimated by \code{\link{pm}}.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @return Estimate of the response rate.
#' @references Guo, H. Y. and Liu, A. A simple and efficient bias-reduced estimator of
#' response probability following a group sequential phase II trial.
#' \emph{J Biopharm Stat},   2005, 15, 773-781.
#' @seealso \code{\link{sbias}}, \code{\link{pm}}, \code{\link{pu}}, \code{\link{pp}} and
#' \code{\link{pk}}.
#' @export
#' @examples
#' pg(21, 19, 4, 54)
#' @author Arsenio Nhacolo
pg <- function(s, n1, r1, n){
  p <- pm(s,n)
  return(p-sbias(n1,r1,n,p))
}
