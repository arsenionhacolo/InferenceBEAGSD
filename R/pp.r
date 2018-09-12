#' @title UMVCUE
#' @description Calculates the uniformly minimum variance conditionally unbiased estimator
#' (UMVCUE) of the true response probability.
#' @details The UMVCUE (\emph{Pepe et al., 2009}) is conditional on on proceeding to the
#' second stage. he sample proportion is used when the trial stopped at first stage.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @return Estimate of the response rate.
#' @references Pepe, M. S.; Feng, Z.; Longton, G. and Koopmeiners, J. Conditional estimation
#' of sensitivity and specificity from a phase 2 biomarker study allowing early termination
#' for futility. \emph{Stat Med}, 2009, 28, 762-779.
#' @seealso \code{\link{pm}}, \code{\link{pg}}, \code{\link{pu}} and
#' \code{\link{pk}}.
#' @export
#' @examples
#' pp(21, 19, 4, 54)
#' @author Arsenio Nhacolo
pp <- function(s, n1, r1, n){
  if (s<=r1){
    return(pm(s,n1))
  }else{
    n2 <- n-n1
    x1 <- max(r1+1,s-n2):min(s,n1)
    return(sum(choose(n1,x1)*choose(n2-1,s-x1-1))/sum(choose(n1,x1)*choose(n2,s-x1)))
  }
}