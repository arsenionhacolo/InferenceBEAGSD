# A (1-alpha)100% CI is collection of pi0 such that the coresponding p-value is
# within [alpha,1-alpha] (Koyama and Chen, 2008). To allow consistence with
# one-sided hypothesis test, for alpha=a we calculate (1-2a)100% two-sided CI.
#' @title Confidence interval (using Method 2 of Nhacolo and Brannath, 2018).
#' @description \code{ci2} computes confidence interval.
#' @details This CI is obtained using the Method 2, one of the four methods proposed by
#' Nhacolo and Brannath (2018) primarily for single-arm adaptive two-stage group sequential
#' designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @param alpha The significance level.
#' @param twosided If \code{FALSE} (default) a one-sided CI is produced.
#' @return CI is a list with lower and upper bounds.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{ci1}}, \code{\link{ci2v2}}, \code{\link{ci3}}, \code{\link{aop2}}, \code{\link{aop2e}}, \code{\link{pipv2}},
#'  \code{\link{mue2}}.
#' @export
#' @author Arsenio Nhacolo
ci2 <- function(dsgn, x1o, xo, alpha=0.05, twosided=FALSE){
  lower <- pipv2(dsgn=dsgn, x1o=x1o, xo=xo, pv=alpha)
  if (twosided){
    upper <- pipv2(dsgn=dsgn, x1o=x1o, xo=xo, pv=1-alpha)
    #msg <- paste("Two-sided ", (1-2*alpha)*100, "% CI for pi: [", sep = "")
  }else{
    upper <- 1
    #msg <- paste("One-sided ", (1-alpha)*100, "% CI for pi: [", sep = "")
  }
  #cat(msg, lower, ",", upper, "]\n\n", sep = "")
  return(list(lower = lower, upper = upper))
}