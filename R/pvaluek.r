
#' @title P-value
#' @description Calculates p-value for Simon-like designs.
#' @details It is based on the definition of p-value by \emph{Koyama and Chen (2008)}.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @param p0 Response rate under the null hypothesis.
#' @return P-value.
#' @references Koyama, T. and Chen, H. Proper inference from Simon's two-stage designs.
#' \emph{Stat Med}, 2008, 27, 3145-3154.
#' @seealso \code{\link{pquantile}} and \code{\link{pk}}.
#' @export
#' @examples
#' pvaluek(21, 19, 4, 54, 0.2)
#' @author Arsenio Nhacolo
pvaluek <- function(s, n1, r1, n, p0){
  n2 <- n - n1
  if (s <= r1){
    return(pbinom(s - 1, n1, p0, lower.tail = F))
  }else{
    x1 <- (r1+1):n1
    return(sum(dbinom(x1, n1, p0) * (pbinom(s - x1 - 1, n2, p0, lower.tail = F))))
  }
}
