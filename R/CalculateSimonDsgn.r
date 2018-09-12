
# Code for simulating and analysing Simon's two-stage optimal and minimax phase II designs
#require(OneArmPhaseTwoStudy) # for functions to get Simon's designs

#' @title Simon's designs 
#' @description \code{CalculateSimonDsgn} finds Simon's optimal and minimax designs.
#' @details Simon's designs are two-stage single-arm for phase II clinical trials.
#' They consist in first stage and overall sample sizes and critical values,
#' \code{n1} and \code{n}, and \code{r1} and \code{r}, respectively.
#' @param p0 The response rate under the null hypothesis.
#' @param p1 The response rate under the alternative hypothesis.
#' @param alpha Type I error rate.
#' @param beta Type II error rate.
#' @param verbose If TRUE (default) the designs are printed (gives messy
#' printout when the function is run without assigment).
#' @return A two-row dataframe containing the optimal and the minimax designs.
#' @references Simon, R. Optimal two-stage designs for phase II clinical trials.
#' \emph{Control Clin Trials}, 1989, 10, 1-10.
#' @seealso \code{\link{SimulateSimonDsgn}} and \code{\link{SimulateSimonDsgnAdaptN}}.
#' @export
#' @importFrom OneArmPhaseTwoStudy setupSimon
#' @importFrom OneArmPhaseTwoStudy getSolutions
#' @examples
#' d <- CalculateSimonDsgn(0.2, 0.4, 0.05, 0.1)
#' @author Arsenio Nhacolo
CalculateSimonDsgn <- function(p0, p1, alpha, beta, verbose = TRUE){
  setUp <- setupSimon(alpha, beta, p0, p1)
  designs <- getSolutions(setUp)$Solutions
  designs <- designs[designs$Type %in% c("MiniMax","Optimal"),]
  designs$targetAlpha <- alpha
  designs$targetBeta <- beta
  if (verbose) print(designs[,c(17,2:13,18,19)])
  return(designs)
}

