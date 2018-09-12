
#' @title Sample proportion
#' @description Calculates the sample proportion.
#' @details For fixed designs the sample propotion is an unbiased (maximum likelihood)
#' estimator of the response rate, but in group sequential designs (e.g., Simon's)
#' it is biased.
#' @param s Total number of successes.
#' @param n Total sample size.
#' @return Estimate of the response rate.
#' @seealso \code{\link{pg}}, \code{\link{pu}}, \code{\link{pp}} and \code{\link{pk}}.
#' @export
#' @examples
#' pm(21, 54)
#' @author Arsenio Nhacolo
pm <- function(s, n){
  return(s/n)
}
