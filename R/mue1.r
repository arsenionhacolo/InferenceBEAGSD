
#Median unbiased estimate of response rate pi
#' @title Median estimate (using Method 1 of Nhacolo and Brannath, 2018).
#' @description \code{mue1} calculates the median estimate of the response rate.
#' @details This estimate is obtained using the Method 1, one of the four methods proposed by
#' Nhacolo and Brannath (2018) primarily for single-arm adaptive two-stage group sequential
#' designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @return Median estimate of response probability.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{mue2}}, \code{\link{mue2v2}}, \code{\link{mue3}}, \code{\link{aop1}}, \code{\link{aop1e}}, \code{\link{pipv1}}.
#' @export
#' @author Arsenio Nhacolo
mue1 <- function(dsgn, x1o, xo){
  if (x1o == 0) return(0)
  pipv1(dsgn=dsgn, x1o=x1o, xo=xo, pv=0.5)
}
