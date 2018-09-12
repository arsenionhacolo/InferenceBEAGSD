#' @title Phase III sample size adjustment factor (Method 2).
#' @description \code{adjustMet2} calculates the multiplicative ajustment  factor \eqn{\rho} to be applied to Phase III
#' sample size estimate using Method 2 proposed by Nhacolo and Brannath (in press).
#' @details The aim of the adjustment is to get an adequately powered Phase III trial based
#' on Phase II data. \eqn{\rho} is found using numerical rearch. See the documentation of the function
#' \code{\link{AnIItoIIIRe}} for more details about the designs.
#' @param p2d Dataframe with Phase II design, with similar as in \code{\link{EKOptAdaptDesigns}}.
#' @param p2r Dataframe containing results of Phase II trials following the design \code{p2d}. It
#' is the output of the function \code{\link{AnalyzeEKOAD}}.
#' @param p2e Phase II estimate to consider among the estimates used by code{\link{AnalyzeEKOAD}}. It
#' can be \code{"pip"} (naive MLE) or one of the four estimates from methods proposed by
#' Nhacolo and Brannath (2018): \code{"pim1"}, \code{"pim2"}, \code{"pim2v2"} or \code{"pim3"}.
#' @param p2p0 Phase II response rate under \eqn{H_0}. If \code{NULL} (default), the value is taken \code{p2d}.
#' @param p2p1 Phase II response rate under \eqn{H_1}. If \code{NULL} (default), the value is taken \code{p2d}.
#' @param p2a Phase II type I error rate. If \code{NULL} (default), the value is taken \code{p2d}.
#' @param p2b Phase II type II error rate. If \code{NULL} (default), the value is taken \code{p2d}.
#' @param p3p0 Phase III response rate of the control group. If \code{NULL} (default), the value is set to \code{p2p0}.
#' @param p3p1 Phase III response rate of the treatment group. If \code{NULL} (default), the value is set to \code{p2p1}.
#' @param p3a Phase III type I error rate. If \code{NULL} (default), the value is set \code{p2a}.
#' @param p3b hase III type II error rate. If \code{NULL} (default), the value is set \code{p2b}.
#' @param nsimul Number of (parametric) bootstrap samples (default 5000).
#' @param seed Seed for random number generator. If \code{NULL} (default), no seed is set.
#' @param rhorange A vector specifying a range to serach for \eqn{\rho}. The default is \code{c(0.5,5)}.
#' @param p3mpt Tolerated error margin for the power, i.e., maximum allowed absolute diference between the estimated
#' expected power and the target. The default is 0.001.
#' @param rhot Search for \eqn{\rho} is interrupted and deem unsuccessful if the absolute difference between
#' current and the previous is less than or equal to \code{rhot}.
#' @return A list containing two dataframes \code{final} and \code{intermed}. \code{final} contains the final
#' measures for the adjustment factor (\eqn{\rho}), and for the unadjusted and adjusted power. \code{intermed}
#' holds the intermediate results (of each bootstrap sample).
#' @references
#' Nhacolo, A. and Brannath, W. Using Estimates from Adaptive Phase II Oncology Trials to Plan Phase III Trials.
#' \emph{In press}.
#'
#' Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#'
#' Ahn, C., Heo, M. and Zhang, S. \emph{Sample Size Calculations for Clustered and Longitudinal Outcomes in Clinical Research}.
#' CRC Press, 2014.
#'
#' @seealso \code{\link{adjustMet2}}, \code{\link{SimulateEKOAD}}, \code{\link{AnalyzeEKOAD}}.
#' @export
#' @examples
#' \dontrun{
#' vdid <- c(6,10) # design ids
#' vp2est <- c("pip","pim1","pim2","pim2v2","pim3")
#' nse <- 1000#number of simulations for each phase
#' cur <- 1; tot <- length(vdid)*length(vp2est)
#' for (did in vdid){
#'  for (p2est in vp2est){
#'    cat('Processing ',cur,' of ',tot,' (',100*round(cur/tot,1),'%)\n',sep = '')
#'     load(paste0("p2r",did,".rdata")) # output of the function AnalyzeEKOAD
#'     out <- adjustMet2(p2d = EKOADwn[EKOADwn$id==did,], p2r = rslt[1:nse,], p2e = p2est, nsimul = nse, seed = 3343)
#'     write.csv(out$final,file = paste0("final",did,p2est,".csv"),row.names = FALSE)
#'    write.csv(out$intermed,file = paste0("intermed",did,p2est,".csv"),row.names = FALSE)
#'    cur <- cur+1
#'  }
#' }
#'
#'
#' vdid <- c(6,10)
#' vp2est <- c("pip","pim1","pim2","pim2v2","pim3")
#' fa <- data.frame()
#' for (did in vdid)
#' {
#'  for (p2est in vp2est){
#'     f <- read.csv(paste0("final",did,p2est,".csv"))
#'     fn <- names(f)
#'     f$dsgn <- did
#'     f <- f[,c('dsgn',fn)]
#'     fa <- rbind(fa,f)
#'   }
#' }
#' write.csv(fa,file = "final_all.csv",row.names = FALSE)
#' }
#' @author Arsenio Nhacolo
adjustMet2 <- function(p2d,p2r,p2e,p2p0=NULL,p2p1=NULL,p2a=NULL,p2b=NULL,p3p0=NULL,p3p1=NULL,
                           p3a=NULL,p3b=NULL, nsimul=5000, seed=NULL,rhorange=c(0.5,5),p3mpt=0.001,
                           rhot=0.0001){ #old name: rhoSimulBootRe
  stopifnot(all(p2r$pi1==p2r$spi1))
  stopifnot(p2d$pi0[1]==p2r$pi0[1],p2d$pi1[1]==p2r$pi1[1],p2d$alpha[1]==p2r$alpha[1],p2d$beta[1]==p2r$beta[1])
  set.seed(seed)

  if (is.null(p2p0)) p2p0 <- p2d$pi0[1]
  if (is.null(p2a)) p2a <- p2d$alpha[1]
  if (is.null(p2b)) p2b <- p2d$beta[1]

  if (is.null(p3p0)) p3p0 <-p2p0
  if (is.null(p3p1)) p3p1 <- p2d$pi1[1]
  if (is.null(p3a)) p3a <- p2a
  if (is.null(p3b)) p3b <-  p2b

  l1 <- min(p2d$x1[p2d$D>0])-1
  u1 <- max(p2d$x1[p2d$D<1])+1
  n1 <- p2d$n1[1]

  nr <- nrow(p2r)
  cna <- rep(NA,nr)
  iout <- data.frame(rho=cna,n=cna,nrho=cna,pwr=cna,pwrrho=cna)#intermediary output


  if (p2e=="pip"){
    for (j in 1:nr){
      r <- data.frame(x1=rep(NA,nsimul))
      phat <- p2r[j,p2e]
      r$x1 <- rbinom(n=nsimul,size=n1,prob=phat)
      r <- merge(r,p2d[,c("x1","n2","D","l")],by="x1",all.x=TRUE,all.y=FALSE)
      r$x2 <- rbinom(n=nrow(r),size=r$n2,prob=phat)
      r$x <- r$x1+r$x2; r$n <- n1+r$n2
      r <- r[r$x1>l1 & (r$x1>=u1|r$x>r$l),]#exclude unseccessful trials
      r$pstar <- r$x/r$n #  Naive MLE
      r$p3nstar <- NA
      rnr <- nrow(r)
      pwrhat <- rep(NA,rnr)
      p3pwr=1-p3b
      for (i in 1:rnr){
        r$p3nstar[i] <- Nct(pc=p3p0,pt=r$pstar[i],alp=p3a,pow=p3pwr)
        pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=r$p3nstar[i],alp=p3a)
      }
      mph <- mean(pwrhat)
      pdif <- p3pwr-mph
      if (abs(pdif)<=p3mpt){
        iout$rho[j] <- 1
      }else if (pdif>0){
        nstar <- r$p3nstar*rhorange[2]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[2]
        }else if (pdif>0){
          warning("Need to increase the upper bound of rho range!")
        }else{
          rl <- 1
          ru <- rhorange[2]
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }else{
        nstar <- r$p3nstar*rhorange[1]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[1]
        }else if (pdif<0){
          warning("Need to decrease the lower bound of rho range!")
        }else{
          rl <- rhorange[1]
          ru <- 1
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }
    }
  }else if (p2e=="pim1"){
    for (j in 1:nr){
      r <- data.frame(x1=rep(NA,nsimul))
      phat <- p2r[j,p2e]
      r$x1 <- rbinom(n=nsimul,size=n1,prob=phat)
      r <- merge(r,p2d[,c("x1","n2","D","l")],by="x1",all.x=TRUE,all.y=FALSE)
      r$x2 <- rbinom(n=nrow(r),size=r$n2,prob=phat)
      r$x <- r$x1+r$x2; r$n <- n1+r$n2
      r <- r[r$x1>l1 & (r$x1>=u1|r$x>r$l),]#exclude unseccessful trials
      r$pstar <- NA #
      r$p3nstar <- NA
      rnr <- nrow(r)
      pwrhat <- rep(NA,rnr)
      p3pwr=1-p3b
      for (i in 1:rnr){
        r$pstar[i] <- mue1(dsgn=p2d, x1o = r$x1[i],xo = r$x[i])
        r$p3nstar[i] <- Nct(pc=p3p0,pt=r$pstar[i],alp=p3a,pow=p3pwr)
        pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=r$p3nstar[i],alp=p3a)
      }
      mph <- mean(pwrhat)
      pdif <- p3pwr-mph
      if (abs(pdif)<=p3mpt){
        iout$rho[j] <- 1
      }else if (pdif>0){
        nstar <- r$p3nstar*rhorange[2]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[2]
        }else if (pdif>0){
          warning("Need to increase the upper bound of rho range!")
        }else{
          rl <- 1
          ru <- rhorange[2]
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }else{
        nstar <- r$p3nstar*rhorange[1]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[1]
        }else if (pdif<0){
          warning("Need to decrease the lower bound of rho range!")
        }else{
          rl <- rhorange[1]
          ru <- 1
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }
    }
  }else if (p2e=="pim2"){
    for (j in 1:nr){
      r <- data.frame(x1=rep(NA,nsimul))
      phat <- p2r[j,p2e]
      r$x1 <- rbinom(n=nsimul,size=n1,prob=phat)
      r <- merge(r,p2d[,c("x1","n2","D","l")],by="x1",all.x=TRUE,all.y=FALSE)
      r$x2 <- rbinom(n=nrow(r),size=r$n2,prob=phat)
      r$x <- r$x1+r$x2; r$n <- n1+r$n2
      r <- r[r$x1>l1 & (r$x1>=u1|r$x>r$l),]#exclude unseccessful trials
      r$pstar <- NA #
      r$p3nstar <- NA
      rnr <- nrow(r)
      pwrhat <- rep(NA,rnr)
      p3pwr=1-p3b
      for (i in 1:rnr){
        r$pstar[i] <- mue2(dsgn=p2d, x1o = r$x1[i],xo = r$x[i])
        r$p3nstar[i] <- Nct(pc=p3p0,pt=r$pstar[i],alp=p3a,pow=p3pwr)
        pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=r$p3nstar[i],alp=p3a)
      }
      mph <- mean(pwrhat)
      pdif <- p3pwr-mph
      if (abs(pdif)<=p3mpt){
        iout$rho[j] <- 1
      }else if (pdif>0){
        nstar <- r$p3nstar*rhorange[2]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[2]
        }else if (pdif>0){
          warning("Need to increase the upper bound of rho range!")
        }else{
          rl <- 1
          ru <- rhorange[2]
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }else{
        nstar <- r$p3nstar*rhorange[1]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[1]
        }else if (pdif<0){
          warning("Need to decrease the lower bound of rho range!")
        }else{
          rl <- rhorange[1]
          ru <- 1
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }
    }
  }else if (p2e=="pim2v2"){
    for (j in 1:nr){
      r <- data.frame(x1=rep(NA,nsimul))
      phat <- p2r[j,p2e]
      r$x1 <- rbinom(n=nsimul,size=n1,prob=phat)
      r <- merge(r,p2d[,c("x1","n2","D","l")],by="x1",all.x=TRUE,all.y=FALSE)
      r$x2 <- rbinom(n=nrow(r),size=r$n2,prob=phat)
      r$x <- r$x1+r$x2; r$n <- n1+r$n2
      r <- r[r$x1>l1 & (r$x1>=u1|r$x>r$l),]#exclude unseccessful trials
      r$pstar <- NA #  Naive MLE
      r$p3nstar <- NA
      rnr <- nrow(r)
      pwrhat <- rep(NA,rnr)
      p3pwr=1-p3b
      for (i in 1:rnr){
        r$pstar[i] <- mue2v2(dsgn=p2d, x1o = r$x1[i],xo = r$x[i])
        r$p3nstar[i] <- Nct(pc=p3p0,pt=r$pstar[i],alp=p3a,pow=p3pwr)
        pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=r$p3nstar[i],alp=p3a)
      }
      mph <- mean(pwrhat)
      pdif <- p3pwr-mph
      if (abs(pdif)<=p3mpt){
        iout$rho[j] <- 1
      }else if (pdif>0){
        nstar <- r$p3nstar*rhorange[2]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[2]
        }else if (pdif>0){
          warning("Need to increase the upper bound of rho range!")
        }else{
          rl <- 1
          ru <- rhorange[2]
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }else{
        nstar <- r$p3nstar*rhorange[1]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[1]
        }else if (pdif<0){
          warning("Need to decrease the lower bound of rho range!")
        }else{
          rl <- rhorange[1]
          ru <- 1
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }
    }
  }else if (p2e=="pim3"){
    for (j in 1:nr){
      r <- data.frame(x1=rep(NA,nsimul))
      phat <- p2r[j,p2e]
      r$x1 <- rbinom(n=nsimul,size=n1,prob=phat)
      r <- merge(r,p2d[,c("x1","n2","D","l")],by="x1",all.x=TRUE,all.y=FALSE)
      r$x2 <- rbinom(n=nrow(r),size=r$n2,prob=phat)
      r$x <- r$x1+r$x2; r$n <- n1+r$n2
      r <- r[r$x1>l1 & (r$x1>=u1|r$x>r$l),]#exclude unsuccessful trials
      r$pstar <- NA #
      r$p3nstar <- NA
      rnr <- nrow(r)
      pwrhat <- rep(NA,rnr)
      p3pwr=1-p3b
      for (i in 1:rnr){
        r$pstar[i] <- mue3(dsgn=p2d, x1o = r$x1[i],xo = r$x[i])
        r$p3nstar[i] <- Nct(pc=p3p0,pt=r$pstar[i],alp=p3a,pow=p3pwr)
        pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=r$p3nstar[i],alp=p3a)
      }
      mph <- mean(pwrhat)
      pdif <- p3pwr-mph
      if (abs(pdif)<=p3mpt){
        iout$rho[j] <- 1
      }else if (pdif>0){
        nstar <- r$p3nstar*rhorange[2]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[2]
        }else if (pdif>0){
          warning("Need to increase the upper bound of rho range!")
        }else{
          rl <- 1
          ru <- rhorange[2]
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }else{
        nstar <- r$p3nstar*rhorange[1]
        for (i in 1:rnr){
          pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
        }
        mph <- mean(pwrhat)
        pdif <- p3pwr-mph
        if (abs(pdif)<=p3mpt){
          iout$rho[j] <- rhorange[1]
        }else if (pdif<0){
          warning("Need to decrease the lower bound of rho range!")
        }else{
          rl <- rhorange[1]
          ru <- 1
          rdf <- ru-rl
          while (abs(pdif)>p3mpt & rdf>rhot) {
            rmd <- rl+(ru-rl)/2
            nstar <- r$p3nstar*rmd
            for (i in 1:rnr){
              pwrhat[i] <- Pwr(pc=p3p0,pt=phat,Nc=nstar[i],alp=p3a)
            }
            mph <- mean(pwrhat)
            pdif <- p3pwr-mph
            if (pdif>0) rl <- rmd else ru <- rmd
            rdf <- ru-rl
          }
          if (rdf<=rhot){
            warning("Maximum iterations reached!")
          }else{
            iout$rho[j] <- rmd
          }
        }
      }
    }
  }else{
    stop("Invalid estimation method (p2e). Valid values are pip, pim1, pim2, pim2v2 and pim3")
  }

  p3pwr=1-p3b
  for (j in 1:nr){
    phat <- p2r[j,p2e]
    iout$n[j] <- Nct(pc=p3p0,pt=phat,alp=p3a,pow=p3pwr)
    iout$nrho[j] <- iout$rho[j]*iout$n[j]
    iout$pwr[j] <- Pwr(pc=p3p0,pt=p3p1,Nc=iout$n[j],alp=p3a)
    iout$pwrrho[j] <- Pwr(pc=p3p0,pt=p3p1,Nc=iout$nrho[j],alp=p3a)
  }

  #a - average, m - median
  fout <- data.frame(est=p2e,rhomean=mean(iout$rho,na.rm=TRUE),rhomedian=median(iout$rho,na.rm=TRUE),rhosd=sd(iout$rho,na.rm=TRUE),
                     nmean=mean(iout$n,na.rm=TRUE),nmedian=median(iout$n,na.rm=TRUE),nsd=sd(iout$n,na.rm=TRUE),
                     nrhomean=mean(iout$nrho,na.rm=TRUE),nrhomedian=median(iout$nrho,na.rm=TRUE),nrhosd=sd(iout$nrho,na.rm=TRUE),
                     pwrmean=mean(iout$pwr,na.rm=TRUE),pwrmedian=median(iout$pwr,na.rm=TRUE),pwrsd=sd(iout$pwr,na.rm=TRUE),
                     pwrrhomean=mean(iout$pwrrho,na.rm=TRUE),pwrrhomedian=median(iout$pwrrho,na.rm=TRUE),pwrrhosd=sd(iout$pwrrho,na.rm=TRUE)
  )
  list("final"=fout, "intermed"=iout)
}

