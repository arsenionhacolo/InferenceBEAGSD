#' @title Response rate to attain a specified p-value (using Method 2 of Nhacolo and Brannath, 2018).
#' @description \code{pipv2} finds the response probability under the null hypotheisis that,
#' given the observed data, would yield a desired overall p-value.
#' @details The p-value is obtained using the Method 2, one of the four methods proposed by
#' Nhacolo and Brannath (2018) primarily for single-arm adaptive two-stage group sequential
#' designs with a binary endpoint.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param x1o The observed stage 1 number of responses.
#' @param xo The total observed number of responses.
#' @param pv The desired p-value.
#' @return Response probability.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{pipv1}}, \code{\link{pipv2v2}}, \code{\link{pipv3}}, \code{\link{aop2}}, \code{\link{aop2e}}.
#' @export
#' @author Arsenio Nhacolo
pipv2 <- function(dsgn, x1o, xo, pv){
  if (x1o == 0) return(NA)
  pi <- -1; pit <- -0.01
  pitv <- c(); cpvv <- c()
  while (pi != pit){
    pit <- pit + 0.01
    if (pit > 0.995) pit <- round(pit)
    cpv <- aop2e(dsgn=dsgn, x1o=x1o, xo=xo, newpi0=pit)
    pitv <- c(pitv, pit)
    cpvv <- c(cpvv, cpv)
    if (cpv > pv){
      pi <-  pit
    }else if (cpv == pv){
      return(pit)
    }else if (pit == 1){
      maxcpv <- max(cpvv)
      cat("WARNING: No pi0 was found that gives pvalue of ", pv,
          ".\n\t The closest pvalue is ", maxcpv, " for pi0 = ", pitv[cpvv==maxcpv], "\n", sep = "")
      return(NA)
    }
  }
  pi <- -1; pit <- pit - 0.01
  while (pi != pit){
    pit <- pit + 0.0001
    cpv <- aop2e(dsgn=dsgn, x1o=x1o, xo=xo, newpi0=pit)
    if (cpv > pv){
      pi <-  pit
    }else if (cpv == pv){
      return(pit)
    }
  }
  pi <- -1; pit <- pit - 0.0001
  while (pi != pit){
    pit <- pit + 0.000001
    cpv <- aop2e(dsgn=dsgn, x1o=x1o, xo=xo, newpi0=pit)
    if (cpv > pv){
      pi <-  pit
    }else if (cpv == pv){
      return(pit)
    }
  }
  return(pi)
}

