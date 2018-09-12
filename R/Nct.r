#' @title Sample size per group for single-stage parallel-group RCT.
#' @description \code{Nct} calculates sample size for one group in an equal-size group
#' two-arm randomized clinical trial with a binary response.
#' @details The sample size is for one group (arm), double the number to get the total.
#' @param pc Response probability in control group.
#' @param pt Response probability in treatment group.
#' @param alp Significance level (default: 0.05).
#' @param pow Power (default: 0.8)
#' @return Sample size for one group.
#' @references Ahn, C., Heo, M. and Zhang, S. \emph{Sample Size Calculations for Clustered and Longitudinal Outcomes in Clinical Research}.
#' CRC Press, 2014.
#' @seealso \code{\link{Pwr}}.
#' @export
#' @examples
#' Nct(0.2,0.3,0.05,0.9)
#' @author Arsenio Nhacolo
Nct <- function(pc,pt,alp=0.05,pow=0.8){
  return((qnorm(1-alp/2)+qnorm(pow))^2*(pt*(1-pt)+pc*(1-pc))/(pt-pc)^2)
}
