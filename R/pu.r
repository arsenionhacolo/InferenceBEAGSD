#' @title UMVUE
#' @description Calculates the uniformly minimum variance unbiased estimator (UMVUE) of the true
#' response probability.
#' @details The UMVUE is based on approach by \emph{Grishick et al. (1946)}. It was first
#' considered by \emph{Chang et al. (1989)} and further studied by \emph{Jung et al. (2004)}.
#' @param s Total number of successes.
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 critical value (trial is stopped at stage 1 if the number of successes
#' is at most \code{r1}).
#' @param n Total sample size.
#' @return Estimate of the response rate.
#' @references Jung, S.-H. and Kim, K. M. On the estimation of the binomial probability in
#' multistage clinical trials. \emph{Stat Med}, 2004, 23, 881-896.
#' @seealso \code{\link{pm}}, \code{\link{pg}}, \code{\link{pp}} and
#' \code{\link{pk}}.
#' @export
#' @examples
#' pu(21, 19, 4, 54)
#' @author Arsenio Nhacolo
pu <- function(s, n1, r1, n){
  if (s<=r1){
    return(pm(s,n1))
  }else{
    n2 <- n-n1
    x1 <- max(r1+1,s-n2):min(s,n1)
    return(sum(choose(n1-1,x1-1)*choose(n2,s-x1))/sum(choose(n1,x1)*choose(n2,s-x1)))
  }
}
