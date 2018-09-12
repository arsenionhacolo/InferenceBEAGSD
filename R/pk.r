#' @title Median unbiased estimator
#' @description Calculates the median unbiased estimator of true response rate
#'  for for Simon-like designs.
#' @details Median unbiased estimator is the value of response rate such that
#' the p-value is 0.5 (\emph{Koyama and Chen, 2008}). The solution is found using
#' numerical search, with a precision of 0.000001.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @param p0 Response rate under the null hypothesis.
#' @return Estimate of the response rate.
#' @references Koyama, T. and Chen, H. Proper inference from Simon's two-stage designs.
#' \emph{Stat Med}, 2008, 27, 3145-3154.
#' @seealso \code{\link{pvaluek}}, \code{\link{pquantile}}, \code{\link{pm}},
#' \code{\link{pg}}, \code{\link{pu}} and \code{\link{pp}}.
#' @export
#' @examples
#' pk(21, 19, 4, 54, 0.2)
#' @author Arsenio Nhacolo
pk <- function(s, n1, r1, n, p0){
  return(pquantile(s, n1, r1, n, p0, 0.5))
}
