
#' @title Value of response rate to attain a given p-value
#' @description Finds, for Simon-like designs, the value of response probability that
#' would yield a given p-value.
#' @details The solution is found using numerical search, with a precision of 0.000001.
#' The p-value is as defined by \emph{Koyama and Chen (2008)}.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @param p0 Response rate under the null hypothesis.
#' @param pvalue The desired p-value.
#' @return Response probability.
#' @references Koyama, T. and Chen, H. Proper inference from Simon's two-stage designs.
#' \emph{Stat Med}, 2008, 27, 3145-3154.
#' @seealso \code{\link{pvaluek}} and \code{\link{pk}}.
#' @export
#' @examples
#' pquantile(21, 19, 4, 54, 0.2, 0.5)
#' @author Arsenio Nhacolo
pquantile <- function(s, n1, r1, n, p0, pvalue){
  if (s == 0) return(NA)
  pk <- -1; pkt <- -0.01
  while (pk != pkt){
    pkt <- pkt + 0.01
    if (pvaluek(s, n1, r1, n, pkt) > pvalue){
      pk <-  pkt
    }else if (pvaluek(s, n1, r1, n, pkt) == pvalue){
      return(pkt)
    }
  }
  pk <- -1; pkt <- pkt - 0.01
  while (pk != pkt){
    pkt <- pkt + 0.0001
    if (pvaluek(s, n1, r1, n, pkt) > pvalue){
      pk <-  pkt
    }else if (pvaluek(s, n1, r1, n, pkt) == pvalue){
      return(pkt)
    }
  }
  pk <- -1; pkt <- pkt - 0.0001
  while (pk != pkt){
    pkt <- pkt + 0.000001
    if (pvaluek(s, n1, r1, n, pkt) > pvalue){
      pk <-  pkt
    }else if (pvaluek(s, n1, r1, n, pkt) == pvalue){
      return(pkt)
    }
  }
  return(pk)
}
