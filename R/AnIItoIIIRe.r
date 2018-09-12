#' @title Use of Phase II estimates to plan Phase III sample size.
#' @description \code{AnIItoIIIRe} calculates the power in a Phase III equal-size group
#' two-arm randomized clinical trial with a binary response planned using estimates
#' from Phase II adaptive two-stage trial.
#' @details The sample size (N) of the Phase III trial is based on the estimates
#' naive MLE and estimators proposed by Nhacolo and Brannath (2018). Different values of retention
#' factor f proposed by Kirby et al. (2012) are applied.The control group response rate is considered
#' to be equal to that under the null hypothesis of the Phase II design, and the hypothesized treatment
#' group response rate considered to be equal to that estimated from the Phase II trial. The target type I
#' error and power are the same as of the Phase II design.Two-sided hypothesis test is assumed.N is a sample
#' size per group, and equal size groups are assume. Hence, N total is 2*N.
#' When calculating the power, the true response rate (in treatment group) is considered to be the one under
#' which the Phase II trial was simulated (spi1).
#' @param rslt Dataframe containing the output from the function \code{\link{AnalyzeEKOAD}}, but with
#' only successful trials (rslt$suco==1), i.e., trials in which \eqn{H_0} was rejected.
#' @param f Vector of length 5 containing multiplicative ajustment factors to be applied
#' to Phase II estimates. The default is \code{f = c(.95,.96,.97,.98,.99)}.
#' @return The input dataframe with corresponding Phase III sample size and power.
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
#' @seealso \code{\link{AnalyzeEKOAD}}, \code{\link{SimulateEKOAD}}, \code{\link{PerforIItoIIIRe}}.
#' @export
#' @author Arsenio Nhacolo
AnIItoIIIRe <- function(rslt, f = c(.95,.96,.97,.98,.99)){
  stopifnot(length(f)==5)

  rslt$Np <- NA
  rslt$Npf1 <- NA
  rslt$Npf2 <- NA
  rslt$Npf3 <- NA
  rslt$Npf4 <- NA
  rslt$Npf5 <- NA

  rslt$Nm1 <- NA
  rslt$Nm1f1 <- NA
  rslt$Nm1f2 <- NA
  rslt$Nm1f3 <- NA
  rslt$Nm1f4 <- NA
  rslt$Nm1f5 <- NA

  rslt$Nm2 <- NA
  rslt$Nm2f1 <- NA
  rslt$Nm2f2 <- NA
  rslt$Nm2f3 <- NA
  rslt$Nm2f4 <- NA
  rslt$Nm2f5 <- NA

  rslt$Nm2v2 <- NA
  rslt$Nm2v2f1 <- NA
  rslt$Nm2v2f2 <- NA
  rslt$Nm2v2f3 <- NA
  rslt$Nm2v2f4 <- NA
  rslt$Nm2v2f5 <- NA

  rslt$Nm3 <- NA
  rslt$Nm3f1 <- NA
  rslt$Nm3f2 <- NA
  rslt$Nm3f3 <- NA
  rslt$Nm3f4 <- NA
  rslt$Nm3f5 <- NA

  rslt$Pp <- NA
  rslt$Ppf1 <- NA
  rslt$Ppf2 <- NA
  rslt$Ppf3 <- NA
  rslt$Ppf4 <- NA
  rslt$Ppf5 <- NA

  rslt$Pm1 <- NA
  rslt$Pm1f1 <- NA
  rslt$Pm1f2 <- NA
  rslt$Pm1f3 <- NA
  rslt$Pm1f4 <- NA
  rslt$Pm1f5 <- NA

  rslt$Pm2 <- NA
  rslt$Pm2f1 <- NA
  rslt$Pm2f2 <- NA
  rslt$Pm2f3 <- NA
  rslt$Pm2f4 <- NA
  rslt$Pm2f5 <- NA

  rslt$Pm2v2 <- NA
  rslt$Pm2v2f1 <- NA
  rslt$Pm2v2f2 <- NA
  rslt$Pm2v2f3 <- NA
  rslt$Pm2v2f4 <- NA
  rslt$Pm2v2f5 <- NA

  rslt$Pm3 <- NA
  rslt$Pm3f1 <- NA
  rslt$Pm3f2 <- NA
  rslt$Pm3f3 <- NA
  rslt$Pm3f4 <- NA
  rslt$Pm3f5 <- NA

  nrows <- nrow(rslt)
  for (i in 1:nrows){
    #cat("Trial ",i," of ",nrows,sep="")

    rslt$Np[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pip[i],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Npf1[i] <-  Nct(pc=rslt$pi0[i],pt=rslt$pip[i]*f[1],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Npf2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pip[i]*f[2],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Npf3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pip[i]*f[3],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Npf4[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pip[i]*f[4],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Npf5[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pip[i]*f[5],alp=rslt$alpha[i],pow=1-rslt$beta[i])

    rslt$Nm1[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim1[i],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm1f1[i] <-  Nct(pc=rslt$pi0[i],pt=rslt$pim1[i]*f[1],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm1f2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim1[i]*f[2],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm1f3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim1[i]*f[3],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm1f4[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim1[i]*f[4],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm1f5[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim1[i]*f[5],alp=rslt$alpha[i],pow=1-rslt$beta[i])

    rslt$Nm2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2[i],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2f1[i] <-  Nct(pc=rslt$pi0[i],pt=rslt$pim2[i]*f[1],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2f2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2[i]*f[2],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2f3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2[i]*f[3],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2f4[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2[i]*f[4],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2f5[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2[i]*f[5],alp=rslt$alpha[i],pow=1-rslt$beta[i])

    rslt$Nm2v2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2v2f1[i] <-  Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i]*f[1],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2v2f2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i]*f[2],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2v2f3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i]*f[3],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2v2f4[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i]*f[4],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm2v2f5[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim2v2[i]*f[5],alp=rslt$alpha[i],pow=1-rslt$beta[i])

    rslt$Nm3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim3[i],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm3f1[i] <-  Nct(pc=rslt$pi0[i],pt=rslt$pim3[i]*f[1],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm3f2[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim3[i]*f[2],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm3f3[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim3[i]*f[3],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm3f4[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim3[i]*f[4],alp=rslt$alpha[i],pow=1-rslt$beta[i])
    rslt$Nm3f5[i] <- Nct(pc=rslt$pi0[i],pt=rslt$pim3[i]*f[5],alp=rslt$alpha[i],pow=1-rslt$beta[i])

    #Power
    rslt$Pp[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Np[i])
    rslt$Ppf1[i] <-  Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Npf1[i])
    rslt$Ppf2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Npf2[i])
    rslt$Ppf3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Npf3[i])
    rslt$Ppf4[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Npf4[i])
    rslt$Ppf5[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Npf5[i])

    rslt$Pm1[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1[i])
    rslt$Pm1f1[i] <-  Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1f1[i])
    rslt$Pm1f2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1f2[i])
    rslt$Pm1f3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1f3[i])
    rslt$Pm1f4[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1f4[i])
    rslt$Pm1f5[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm1f5[i])

    rslt$Pm2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2[i])
    rslt$Pm2f1[i] <-  Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2f1[i])
    rslt$Pm2f2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2f2[i])
    rslt$Pm2f3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2f3[i])
    rslt$Pm2f4[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2f4[i])
    rslt$Pm2f5[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2f5[i])

    rslt$Pm2v2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2[i])
    rslt$Pm2v2f1[i] <-  Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2f1[i])
    rslt$Pm2v2f2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2f2[i])
    rslt$Pm2v2f3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2f3[i])
    rslt$Pm2v2f4[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2f4[i])
    rslt$Pm2v2f5[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm2v2f5[i])

    rslt$Pm3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3[i])
    rslt$Pm3f1[i] <-  Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3f1[i])
    rslt$Pm3f2[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3f2[i])
    rslt$Pm3f3[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3f3[i])
    rslt$Pm3f4[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3f4[i])
    rslt$Pm3f5[i] <- Pwr(pc=rslt$pi0[i],pt=rslt$spi1[i],alp=rslt$alpha[i],Nc=rslt$Nm3f5[i])

  }
  return(rslt)
}

