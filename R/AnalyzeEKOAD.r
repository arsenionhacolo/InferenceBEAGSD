
#' @title Analyze simulated adaptive trials.
#' @description \code{AnalyzeEKOAD} performs inference on trials simulated by the function
#' \code{\link{SimulateEKOAD}} using the methods proposed by Nhacolo and Brannath (2018)
#' and naive maximum likelihood.
#' @details Overall p-values, point estimates and confidence intervals are calculated.
#' @param replicates Number of simulated trials to be analysed. If \code{NULL} (default),
#' all trials found in \code{./basedir/SimulatedTrials} are analysed.
#' @param basedir The base directory containing the sub-directory \code{SimulatedTrials}
#' with the simulated trials. If \code{NULL} (default), the current working directory
#' is uded.
#' @return A dataframe with the results. A copy is saved in the file \code{Results.csv}
#' in the \code{basedir}.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{SimulateEKOAD}, \code{\link{mue1}}, \code{\link{mue2}}, \code{\link{mue2v2}}, \code{\link{mue3}}}.
#' @export
#' @author Arsenio Nhacolo
AnalyzeEKOAD <- function(replicates = NULL, basedir = NULL){
  if (!is.null(basedir)) setwd(basedir)
  dsgn <- read.csv("Design.csv")
  files <- length(list.files("./SimulatedTrials"))
  if (is.null(replicates)){
    replicates <- files
  }else if (replicates > files){
    replicates <- files
    cat("The specified number of replicates is greater than the available (",
        files, "), ... using the available!", sep = "")
  }
  cat("Analyzing ", replicates, " trials... ", sep = "")
  trial <- read.csv("./SimulatedTrials/trial1.csv")
  s1 = sum(trial$resp[trial$interim==1])
  rslt <- with(trial, data.frame(trial = 1, s1 = s1, s = sum(resp),
                                 n1 = length(resp[interim==1]),
                                 n = length(resp), stop = stop[1],
                                 l = dsgn$l[dsgn$x1==s1],
                                 D = dsgn$D[dsgn$x1==s1]))
  for (i in 2:replicates){
    trial <- read.csv(paste("./SimulatedTrials/trial", i, ".csv", sep = ""))
    s1 = sum(trial$resp[trial$interim==1])
    rslt <- with(trial, rbind(rslt, data.frame(trial = i, s1 = s1, s = sum(resp),
                                               n1 = length(resp[interim==1]),
                                               n = length(resp), stop = stop[1],
                                               l = dsgn$l[dsgn$x1==s1],
                                               D = dsgn$D[dsgn$x1==s1])))
  }
  rslt$s2 <- rslt$s - rslt$s1
  rslt$n2 <- rslt$n - rslt$n1
  rslt$l1 <- dsgn$l1[1]
  rslt$u1 <- dsgn$u1[1]
  rslt$pi0 <- dsgn$pi0[1]
  rslt$pi1 <- dsgn$pi1[1]
  rslt$spi1 <- dsgn$spi1[1]
  rslt$dsgnn1 <- dsgn$n1[1]
  rslt$alpha <- dsgn$alpha[1]
  rslt$beta <- dsgn$beta[1]
  #Success (H0 rejected)
  rslt$suco <- as.numeric(rslt$s1 > rslt$l1 & (rslt$s1 >= rslt$u1 | rslt$s > rslt$l)) #original decision rule
  # Naive estimation of the response rate p
  rslt$pip1 <- rslt$s1/rslt$n1 #  Sample proportion (stage 1 only)
  rslt$pip2 <- NA
  rslt$pip2[rslt$n2 > 0] <- rslt$s2[rslt$n2 > 0]/rslt$n2[rslt$n2 > 0] #Sample proportion (stage 2 only)
  rslt$pip <- rslt$s/rslt$n #  Sample proportion

  #Overall p-value and estimation based on methods 1, 2, 2v2 and 3
  rslt$pvm1 <- NA  # Overall p-value
  rslt$pim1 <- NA # Median estimate
  rslt$cilm1 <- NA # Lower bound of CI
  rslt$cium1 <- NA # Upper bound of CI

  rslt$pvm2 <- NA  # Overall p-value
  rslt$pim2 <- NA # Median estimate
  rslt$cilm2 <- NA # Lower bound of CI
  rslt$cium2 <- NA # Upper bound of CI

  rslt$pvm2v2 <- NA  # Overall p-value
  rslt$pim2v2 <- NA # Median estimate
  rslt$cilm2v2 <- NA # Lower bound of CI
  rslt$cium2v2 <- NA # Upper bound of CI

  rslt$pvm3 <- NA  # Overall p-value
  rslt$pim3 <- NA # Median estimate
  rslt$cilm3 <- NA # Lower bound of CI
  rslt$cium3 <- NA # Upper bound of CI

  for (i in 1:nrow(rslt)){
    rslt$pvm1[i] <- aop1(dsgn = dsgn, x1o = rslt$s1[i],xo = rslt$s[i], verbose = FALSE)
    rslt$pim1[i] <- mue1(dsgn=dsgn, x1o = rslt$s1[i],xo = rslt$s[i])
    cim1 <- ci1(dsgn = dsgn, x1o = rslt$s1[i], xo =  rslt$s[i], alpha = rslt$alpha[1], twosided = FALSE)
    rslt$cilm1[i] <- cim1$lower
    rslt$cium1[i] <- cim1$upper

    rslt$pvm2[i] <- aop2(dsgn = dsgn, x1o = rslt$s1[i],xo = rslt$s[i], verbose = FALSE)
    rslt$pim2[i] <- mue2(dsgn=dsgn, x1o = rslt$s1[i],xo = rslt$s[i])
    cim2 <- ci2(dsgn = dsgn, x1o = rslt$s1[i], xo =  rslt$s[i], alpha = rslt$alpha[1], twosided = FALSE)
    rslt$cilm2[i] <- cim2$lower
    rslt$cium2[i] <- cim2$upper

    rslt$pvm2v2[i] <- aop2v2(dsgn = dsgn, x1o = rslt$s1[i],xo = rslt$s[i], verbose = FALSE)
    rslt$pim2v2[i] <- mue2v2(dsgn=dsgn, x1o = rslt$s1[i],xo = rslt$s[i])
    cim2v2 <- ci2v2(dsgn = dsgn, x1o = rslt$s1[i], xo =  rslt$s[i], alpha = rslt$alpha[1], twosided = FALSE)
    rslt$cilm2v2[i] <- cim2v2$lower
    rslt$cium2v2[i] <- cim2v2$upper

    rslt$pvm3[i] <- aop3e(dsgn = dsgn, x1o = rslt$s1[i],xo = rslt$s[i])
    rslt$pim3[i] <- mue3(dsgn=dsgn, x1o = rslt$s1[i],xo = rslt$s[i])
    cim3 <- ci3(dsgn = dsgn, x1o = rslt$s1[i], xo =  rslt$s[i], alpha = rslt$alpha[1], twosided = FALSE)
    rslt$cilm3[i] <- cim3$lower
    rslt$cium3[i] <- cim3$upper
  }
  rslt$sucm1 <- as.numeric(rslt$pvm1 <= rslt$alpha) #Success (H0 rejected)
  rslt$covm1 <- as.numeric(rslt$cilm1<=rslt$spi1 & rslt$spi1<=rslt$cium1)# Is the true response rate covered by the CI?

  rslt$sucm2 <- as.numeric(rslt$pvm2 <= rslt$alpha) #Success (H0 rejected)
  rslt$covm2 <- as.numeric(rslt$cilm2<=rslt$spi1 & rslt$spi1<=rslt$cium2)# Is the true response rate covered by the CI?

  rslt$sucm2v2 <- as.numeric(rslt$pvm2v2 <= rslt$alpha) #Success (H0 rejected)
  rslt$covm2v2 <- as.numeric(rslt$cilm2v2<=rslt$spi1 & rslt$spi1<=rslt$cium2v2)# Is the true response rate covered by the CI?

  rslt$sucm3 <- as.numeric(rslt$pvm3 <= rslt$alpha) #Success (H0 rejected)
  rslt$covm3 <- as.numeric(rslt$cilm3<=rslt$spi1 & rslt$spi1<=rslt$cium3)# Is the true response rate covered by the CI?

  write.csv(rslt, "Results.csv", row.names = FALSE)
  cat("Done. Results saved in Results.csv\n", sep = "")
  return(rslt)
}
