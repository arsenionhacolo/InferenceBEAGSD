
#' @title Probability mass function of \code{(M, S)}
#' @description Probability mass function of \code{M} (stage) and \code{S} (number of successes).
#' @details Probability mass function of the statistic \code{(M, S)} for Simon-like
#' designs (allowing early stopping for futility only).
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @param p True success probability.
#' @param m Stage number (1 or 2). It is automatically determined based on \code{s} and \code{r1},
#' therefore it shouldn't be provided, unless there are reasons to do so.
#' @return Density.
#' @references Jung, S.-H. and Kim, K. M. On the estimation of the binomial probability in
#' multistage clinical trials. \emph{Stat Med}, 2004, 23, 881-896.
#' @seealso \code{\link{sbias}} and \code{\link{pg}}.
#' @export
#' @examples
#' sfms(21, 19, 4, 54, 0.4)
#' @author Arsenio Nhacolo
sfms <- function(s, n1, r1, n, p, m = NA){
  if (is.na(m)) m <- ifelse(s <= r1,1,2)
  if (m==1){
    dens <- choose(n1,s) * p^s * (1-p)^(n1-s)
  }else if (m==2){
    n2 <- n-n1
    x1 <- max(r1+1,s-n2):min(s,n1)
    dens <- sum(choose(n1,x1)*choose(n2,s-x1))*p^s*(1-p)^(n-s)
  }else{
    stop("Wrong stage number (m). Allowed values: 1 or 2.")
  }
  return(dens)
}

