
#' @title Pre-process the Englert and Kieser (2013) optimal adaptive designs.
#' @description \code{dsgnPrep} takes Englert and Kieser's optimal adaptive design and adds
#' information that is needed by other functions.
#' @details The function adds, to each x1 leading to 2nd stage, the corresponding p-value (p1)
#' and its back-wards image (p1B), the stage-wise weights w1 and w2 and other information used
#' in inference methods proposed by  Nhacolo and Brannath (2018).
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOptAdaptDesigns}}.
#' @param w1,w2 Stage 1 and 2 weights. If \code{w1="n"} (default), weights a calculated
#' based on stage-wise sample sizes as described in Nhacolo and Brannath (2018).
#' If \code{w1="sr2"}, then \code{w1=w2=1/sqrt(2)}.
#' @return Dataframe containing the input dataframe with added information.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @export
#' @importFrom adaptTest tsT
#' @examples
#' \dontrun{
#' #Designs with w1a and w2 calculated based on sample sizes
#' EKOADwn <- data.frame()
#' for (j in 1:max(EKOptAdaptDesigns$id)){
#'   EKOADwn <- rbind(EKOADwn, dsgnPrep(dsgn = EKOptAdaptDesigns[EKOptAdaptDesigns$id==j,],w1 = "n"))
#' }
#' save(EKOADwn,file = "EKOADwn.RData")
#' }
#' @author Arsenio Nhacolo
dsgnPrep <- function(dsgn=NULL, w1="n", w2=NULL){
  if (w1=="sr2"){
    dsgn$w1 <- dsgn$w2 <- 1/sqrt(2)
  }else if (w1=="n"){
    dsgn$w1 <- sqrt(dsgn$n1/(dsgn$n1+dsgn$n2))
    dsgn$w2 <- sqrt(dsgn$n2/(dsgn$n1+dsgn$n2))
  }else if (mode(w1)!="numeric" | mode(w2)!="numeric"){
    stop("Invalid w1 or w2")
  }else if (abs(w1^2+w2^2-1) > 0.0001){
    stop("w1^2+w2^2 must be equal to 1.")
  }
  l1 <- min(dsgn$x1[dsgn$D>0])-1
  u1 <- max(dsgn$x1[dsgn$D<1])+1
  dsgn$alpha0 <- 1-pbinom(l1-1,dsgn$n1[1], dsgn$pi0[1])#RE-CHECK HERE (l1 or l1-1 ?)
  dsgn$alpha1 <- 1-pbinom(u1-1,dsgn$n1[1],dsgn$pi0[1])
  dsgn$p1 <- 1-pbinom(dsgn$x1-1,dsgn$n1,dsgn$pi0)
  alpha2 <- adaptTest::tsT(typ = "l",a = dsgn$alpha[1],a0 = dsgn$alpha0[1],a1 = dsgn$alpha1[1],a2 = NA)
  dsgn$p1B <- dsgn$p1
  ds <- dsgn[dsgn$x1>l1 & dsgn$x1<u1,]
  dsgn$p1B[dsgn$x1>l1 & dsgn$x1<u1] <- pnorm((qnorm(alpha2)-ds$w2*qnorm(ds$D))/ds$w1)
  dsgn$l1 <- l1
  dsgn$u1 <- u1
  dsgn$alpha2 <- alpha2
  return(dsgn)
}
