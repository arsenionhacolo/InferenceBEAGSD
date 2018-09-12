#' @title Overall p-value (Method 2 of Nhacolo and Brannath, 2018).
#' @description \code{aop2} calculates the overall p-value for adaptive two-stage designs
#' with binary endpoint using the Method 2 (see Nhacolo and Brannath, 2018).
#' @details This is one of the four methods proposed by Nhacolo and Brannath (2018) primarily
#' for single-arm adaptive two-stage group sequential designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @param verbose If \code{TRUE} (default) messages will be printed.
#' @return p-value.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{aop1}}, \code{\link{aop2v2}}, \code{\link{aop3e}}
#' @export
#' @author Arsenio Nhacolo
aop2 <- function(dsgn, x1o, xo, verbose = TRUE){
  stopifnot(xo >= x1o)
  l1 <- dsgn$l1[1]
  u1 <- dsgn$u1[1]
  n1 <- dsgn$n1[1]
  pi0 <- dsgn$pi0[1]
  if (x1o<=l1 | x1o>=u1){
    pvalue <- 1-pbinom(x1o-1,n1,pi0)
  }else{
    lx1o <- dsgn$l[dsgn$x1==x1o]
    n2o <- dsgn$n2[dsgn$x1==x1o]
    p2o <- 1-pbinom(xo-x1o-1,n2o,pi0)
    Do <- dsgn$D[dsgn$x1==x1o]
    d <- dsgn[dsgn$x1>l1 & dsgn$x1<u1,]
    DDo <- d$D-Do+p2o
    DDo[DDo<0] <- 0
    DDo[DDo>1] <- 1
    pvalue <- 1-pbinom(u1-1,n1,pi0)+
      sum(dbinom(d$x1,n1,pi0)*(DDo))
  }
  if (verbose){
    cat("Design: (pi0, pi1, alpha, beta, n1) = (",pi0,", ",dsgn$pi1[1],", ",
        dsgn$alpha[1],", ",dsgn$beta[1],", ",dsgn$n1[1],")\nP-value: ",pvalue,"\n", sep = "")
  }
  return(pvalue)
}

