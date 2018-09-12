
#' @title Power for single-stage parallel-group RCT.
#' @description \code{Pwr} calculates the power in an equal-size group
#' two-arm randomized clinical trial with a binary response.
#' @param pc Response probability in control group.
#' @param pt Response probability in treatment group.
#' @param Nc Sample size per group.
#' @param alp Significance level (default: 0.05).
#' @return Sample size for one group.
#' @references Ahn, C., Heo, M. and Zhang, S. \emph{Sample Size Calculations for Clustered and Longitudinal Outcomes in Clinical Research}.
#' CRC Press, 2014.
#' @seealso \code{\link{Nct}}.
#' @export
#' @examples
#' Pwr(0.2,0.3,389,0.05)
#' @author Arsenio Nhacolo
Pwr <- function(pc,pt,Nc,alp=0.05){
  Nt <- Nc
  pnorm(abs(pt-pc)/(sqrt(pt*(1-pt)/Nt+pc*(1-pc)/Nc))-qnorm(1-alp/2))
}

