#' @title Overall p-value (Method 3 of Nhacolo and Brannath, 2018).
#' @description \code{aop3e} calculates the overall p-value for adaptive two-stage designs
#' with binary endpoint using the Method 3 (see Nhacolo and Brannath, 2018).
#' @details This is one of the four methods proposed by Nhacolo and Brannath (2018) primarily
#' for single-arm adaptive two-stage group sequential designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @param newpi0 New response probability that replaces the one under the null hypothesis.
#' Omit it if the intention is only to calculate the overall p-value.
#' @return p-value.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{aop1}}, \code{\link{aop1e}}, \code{\link{aop2}}, \code{\link{aop2e}}, \code{\link{aop2v2}}, \code{\link{aop2ev2}}
#' @export
#' @author Arsenio Nhacolo
aop3e <- function(dsgn, x1o, xo, newpi0=NULL){
  stopifnot(xo >= x1o)
  l1 <- dsgn$l1[1]
  u1 <- dsgn$u1[1]
  n1 <- dsgn$n1[1]
  if (is.null(newpi0)) pi0 <- dsgn$pi0[1] else pi0 <- newpi0
  if (x1o<=l1 | x1o>=u1){
    pvalue <- 1-pbinom(x1o-1,n1,pi0)
  }else{
    d <- dsgn[dsgn$x1>l1 & dsgn$x1<u1,]
    #n2o <- d$n2[d$x1==x1o]
    p2o <- 1-pbinom(xo-x1o-1,d$n2[d$x1==x1o],pi0)
    #w1o <- dsgn$w1[dsgn$x1==x1o]
    #w2o <- dsgn$w2[dsgn$x1==x1o]
    z <- (d$w1[d$x1==x1o]*qnorm(1-d$p1B[d$x1==x1o]) + d$w2[d$x1==x1o]*qnorm(1-p2o) - d$w1*qnorm(1-d$p1B))/d$w2
    zp <- 1-pnorm(z)
    pvalue <- 1-pbinom(u1-1,n1,pi0) + sum(dbinom(d$x1,n1,pi0)*(zp))
  }
  return(pvalue)
}

