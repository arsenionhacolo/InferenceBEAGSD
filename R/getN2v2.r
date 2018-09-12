#' @title Number of patients to be enrolled in the second stage
#' @description Calculates the number of patients which should be enrolled in the
#' second stage if the conditional power should be altert to "cp". It's a version
#' of \code{\link[OneArmPhaseTwoStudy]{getN2}}.
#' @details This functon is the same as \code{\link[OneArmPhaseTwoStudy]{getN2}}
#' (\code{OneArmPhaseTwoStudy} package), with some changes in arguments' validation.
#' It's is a helper to
#' \code{\link{SimulateSimonDsgnAdaptN}}.
#' @inheritParams OneArmPhaseTwoStudy::getN2
#' @references Englert S., Kieser M. Adaptive designs for single-arm phase II trials in oncology.
#' \emph{Pharm Stat}, 2012, 11, 241-249.
#' @seealso \code{\link[OneArmPhaseTwoStudy]{getN2}}, \code{\link{SimulateSimonDsgnAdaptN}}.
#' @importFrom OneArmPhaseTwoStudy getCP
#' @examples
#' designParam <- CalculateSimonDsgn(0.2, 0.4, 0.05, 0.1)
#' dsgn <- designParam[designParam$Type == "Optimal",]
#' getN2v2(0.9, dsgn$p1, dsgn, 7)
getN2v2 <- function(cp, p1, design, k, mode = 0, alpha = 0.05){
  stopifnot(class(cp) == "numeric", cp > 0, cp <= 1)
  stopifnot(class(p1) == "numeric", p1 > 0, p1 <= 1)
  stopifnot(class(design) == "data.frame")
  dgn_names = names(design)
  stopifnot("r1" %in% dgn_names, "r" %in% dgn_names, "p0" %in%
              dgn_names, "n1" %in% dgn_names, "n" %in% dgn_names)
  stopifnot(class(k) == "numeric" | class(k) == "integer",
            k >= 0)
  stopifnot(class(mode) == "numeric", mode %in% 0:3)
  stopifnot(class(alpha) == "numeric", alpha > 0, alpha <= 1)
  if (cp == 1) {
    return(0)
  }
  cpInt = 0
  n2 = 0
  while (cpInt < cp) {
    n2 = n2 + 1
    cpInt = getCP(n2, p1, design, k, mode, alpha)
  }
  n2
}

