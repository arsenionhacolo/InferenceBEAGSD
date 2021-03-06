#' Englert and Kieser (2013)'s optimal adaptive designs.
#'
#' A dataframe containing all optimal adaptive two-stage designs for phase II
#'cancer clinical trials present in Englert and Kieser (2013).
#'
#' @format A dataframe with 709 rows and 11 variables:
#' \describe{
#'   \item{id}{Identifier of the designs}
#'   \item{x1}{Number of successes (responses) at stage 1}
#'   \item{n2}{Stage 2 sample size}
#'   \item{D}{Discrete conditional error function}
#'   \item{l}{Stage 2 decision boundary}
#'   \item{pi0}{Response probability under the null hypothesis}
#'   \item{pi1}{Response probability under the alternative hypothesis}
#'   \item{alpha}{Type I error rate}
#'   \item{beta}{Type II error rate}
#'   \item{n1}{Stage 1 sample size}
#'   \item{n2max}{Maximum stage 2 sample size}
#' }
#' @source \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201200220}
"EKOptAdaptDesigns"
