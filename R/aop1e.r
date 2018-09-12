#' @title Overall p-value for CI (Method 1 of Nhacolo and Brannath, 2018).
#' @description \code{aop1e} is a modified version of \code{\link{aop1}} used
#' for getting the confidence interval.
#' @details This is one of the four methods proposed by Nhacolo and Brannath (2018) primarily
#' for single-arm adaptive two-stage group sequential designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @param newpi0 New response probability that replaces the one under the null hypothesis.
#' @return p-value.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{aop2e}}, \code{\link{aop2ev2}}, \code{\link{aop3e}}, \code{\link{aop1}}.
#' @export
#' @author Arsenio Nhacolo
aop1e <- function(dsgn, x1o, xo, newpi0){
  stopifnot(xo >= x1o)

  l1 <- dsgn$l1[1]
  u1 <- dsgn$u1[1]
  n1 <- dsgn$n1[1]
  pi0 <- newpi0
  if (x1o<=l1 | x1o>=u1){
    pvalue <- 1-pbinom(x1o-1,n1,pi0)
  }else{
    lx1o <- dsgn$l[dsgn$x1==x1o]
    d <- dsgn[dsgn$x1>l1 & dsgn$x1<u1,]
    pvalue <- 1-pbinom(u1-1,n1,pi0)+
      sum(dbinom(d$x1,n1,pi0)*(1-pbinom(xo-lx1o+d$l-d$x1-1,d$n2,pi0)))
  }
  return(pvalue)
}

