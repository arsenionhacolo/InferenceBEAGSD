#MONOTONICITY CHECK OF DELTA AS FUNCTION OF X1 AND X2, FOR METHOD 3 (INVERSE NORMAL COMB. FUNCTION)
#' @title Check the monotonicity of the sample space ordering.
#' @description \code{checkMonoDCF} checks the monotonicity of the sample space ordering defined
#' based on inverse normal combination function (see Nhacolo and Brannath, 2018).
#' @details The monotonicity is with respect to the stage 2 number of successes.
#' @param d Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param verbose If \code{TRUE} (default) messages about monotonicity will be printed.
#' @return A list containing a dataframe (\code{mono}) with detailed info, and a logical variable
#' \code{notmono} indicating wether non-monotonicity was concluded.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @export
#' @examples
#' \dontrun{
#' #Check for all Englert and Kieser designs
#' notmov <- c()
#' for (i in 1:max(EKOADwn$id)){
#'   notmov <- c(notmov,checkMonoDCF(EKOADwn[EKOADwn$id==1,],verbose=FALSE)[[2]])
#' }
#' isMonotone <- !any(notmov);isMonotone
#' }
#' @author Arsenio Nhacolo
checkMonoDCF <- function(d, verbose=TRUE){
  nm <- names(d)
  if (!('w1'%in%nm|'w2'%in%nm|'p1B'%in%nm)) stop('A preprocessed design (from EKOADwn.RData).')
  id <- d$id[1]
  pi0 <- d$pi0[1]
  d <- d[d$n2!=0,c('x1','n2','w1','w2','p1B')]
  nrow <- nrow(d)
  trow <- sum(d$n2)+nrow
  cna <- rep(NA,trow)
  mono <- data.frame(did=rep(id,trow),x1=cna,x2=cna,delta=cna, dif=cna)
  t <- 0
  for (i in 1:nrow){
    ct <- 0
    n2 <- d$n2[i]
    for (x2 in 0:n2){
      ct <- ct+1
      t <- t+1
      p2 <- 1-pbinom(x2-1,n2,pi0)
      mono$delta[t] <- pnorm(d$w1[i]*qnorm(1-d$p1B[i])+d$w2[i]*qnorm(1-p2))
      mono$x1[t] <- d$x1[i]
      mono$x2[t] <- x2
    }
    mono$dif[(t-ct+2):t] <- mono$delta[(t-ct+2):t]-mono$delta[(t-ct+1):(t-1)]
  }
  notmo <- any(mono$dif<0,na.rm = TRUE)
  if (verbose){
    if (notmo){
      cat('Delta is NOT monotone in x2!\n')
    }else{
      cat('Delta is monotone in x2.\n')
    }
  }
  list('monodf'=mono,'notmono'=notmo)
}
