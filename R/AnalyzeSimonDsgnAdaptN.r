#' @title Analysis of simulated adaptive Simon's design trials
#' @description Analyses the trials simulated by \code{\link{SimulateSimonDsgnAdaptN}}.
#' @details In addition to hypothesis testing, the response rate is estimated
#' using different estimators: \code{\link{pm}}, \code{\link{pg}}, \code{\link{pu}},
#' \code{\link{pp}} and \code{\link{pk}}. The overall critical value, \code{r}, is
#' recalculated using conditional type I error (\emph{Englert and Kieser, 2012}).
#' @param replicates Number of trials to be analysed. By default all simulated trials are
#' analysed.
#' @param basedir The root directory in which simulations were performed. The current
#' working directory is assumed by default. It must contain all the files and folders
#' created by \code{\link{SimulateSimonDsgnAdaptN}}.
#' @return Creates two data files in \code{basedir} containing results for optimal
#' (\emph{ResultsOptimalDesignAdapt.csv}) and minimax (\emph{ResultsMinimaxDesignAdapt.csv}).
#' The files contain a \code{trial} ID, stage 1, stage 2 and overall number of
#' successful responses, \code{s1}, \code{s2} and \code{s}, sample sizes
#' (equal to those pre-specified by design),  \code{n1}, \code{n2} and \code{n},
#' and critical values, \code{r1} and \code{r}. \code{p0} the response rate
#' assumed under \eqn{H_0} and \code{dsgnp1} under \eqn{H_1}. \code{p1} is the
#' true response rate (used for generating trial data). \code{pm1} and \code{pm2}
#' are, respectively, \code{\link{pm}} based only of stage 1 and stage 2 data.
#' \code{stop} indicates whether the trial stopped at first stage (\code{stop = 1}),
#' and \code{success} indicates whether \eqn{H_0} was rejected (\code{success = 1}).
#' @seealso \code{\link{CalculateSimonDsgn}}, \code{\link{SimulateSimonDsgnAdaptN}},
#' \code{\link{AnalyzePerformanceSimon}} and \code{\link{AnalyzeSimonDsgn}}.
#' @export
#' @importFrom OneArmPhaseTwoStudy getCE
#' @importFrom OneArmPhaseTwoStudy getD_proportionally
#' @importFrom OneArmPhaseTwoStudy getD_equally
#' @importFrom OneArmPhaseTwoStudy getD_distributeToOne
#' @importFrom OneArmPhaseTwoStudy get_r2_flex
#' @examples
#' AnalyzeSimonDsgnAdaptN()
#' @author Arsenio Nhacolo
AnalyzeSimonDsgnAdaptN <- function(replicates = NA, basedir = NA){
  if (!is.na(basedir)) setwd(basedir)
  designParam <- read.csv("DesignParametersAdapt.csv")
  # Optimal design
  dsgn <- designParam[designParam$Type == "Optimal",]
  files <- length(list.files("OptimalAdapt/SimulatedTrials"))
  if (is.na(replicates)){
    replicates <- files
  }else if (replicates > files){
    replicates <- files
    cat("The specified number of replicates is greater than the available (",
        files, "), ... using the available!", sep = "")
  }
  trial <- read.csv("OptimalAdapt/SimulatedTrials/trial1.csv")
  result <- with(trial, data.frame(trial = 1, s1 = sum(resp[interim==1]),
                                   s = sum(resp),
                                   r1 = dsgn$r1,
                                   r = dsgn$r,
                                   n1 = length(resp[interim==1]),
                                   n = length(resp), stop = stop[1]))
  for (i in 2:replicates){
    trial <- read.csv(paste("OptimalAdapt/SimulatedTrials/trial", i, ".csv", sep = ""))
    result <- with(trial, rbind(result, data.frame(trial = i, s1 = sum(resp[interim==1]),
                                                   s = sum(resp),
                                                   r1 = dsgn$r1,
                                                   r = dsgn$r,
                                                   n1 = length(resp[interim==1]),
                                                   n = length(resp), stop = stop[1])))
  }
  result$p0 <- dsgn$p0
  result$p1 <- dsgn$p1
  result$dsgnp1 <- dsgn$dsgnp1
  result$alpha <- dsgn$Alpha
  result$targetAlpha <- dsgn$targetAlpha
  result$beta <- dsgn$Beta
  result$targetBeta <- dsgn$targetBeta
  # Calculate new r
  result$newr2 <- NA
  result$newr <- result$r
  for (i in 1:nrow(result)){
    if (result$stop[i] == 0){
      ce <- switch(dsgn$restAlphaMet + 1,
                   getCE(dsgn, result$s1[i]),
                   getD_proportionally(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1],
                   getD_equally(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1],
                   getD_distributeToOne(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1])
      n2 <- as.numeric(result$n[i] - result$n1[i])
      result$newr2[i] <- get_r2_flex(ce, dsgn$p0, n2)
      result$newr[i] <- result$r1[i] + result$newr2[i]
    }

  }

  result$success <- as.numeric(result$s > result$newr)
  result$s2 <- result$s - result$s1
  result$n2 <- result$n - result$n1
  # Estimation of the response rate p
  result$pm1 <- result$s1/result$n1 #  Maximum likelihood (stage 1 only)
  result$pm2 <- NA
  result$pm2[result$n2 > 0] <- result$s2[result$n2 > 0]/result$n2[result$n2 > 0] #  Maximum likelihood (stage 2 only)
  result$pm <- result$s/result$n #  Maximum likelihood
  result$pg <- NA
  result$pu <- NA
  result$pp <- NA
  result$pk <- NA
  for (i in 1:nrow(result)){
    result$pg[i] <- pg(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pu[i] <- pu(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pp[i] <- pp(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pk[i] <- pk(result$s[i], result$n1[i], result$r1[i], result$n[i], result$p0[i])
  }
  write.csv(result, "ResultsOptimalDesignAdapt.csv", row.names = FALSE)

  # Minimax design
  dsgn <- designParam[designParam$Type == "MiniMax",]
  files <- length(list.files("MinimaxAdapt/SimulatedTrials"))
  if (is.na(replicates)){
    replicates <- files
  }else if (replicates > files){
    replicates <- files
    cat("The specified number of replicates is greater than the available (",
        files, "), ... using the available!", sep = "")
  }
  trial <- read.csv("MinimaxAdapt/SimulatedTrials/trial1.csv")
  result <- with(trial, data.frame(trial = 1, s1 = sum(resp[interim==1]),
                                   s = sum(resp),
                                   r1 = dsgn$r1,
                                   r = dsgn$r,
                                   n1 = length(resp[interim==1]),
                                   n = length(resp), stop = stop[1]))
  for (i in 2:replicates){
    trial <- read.csv(paste("MinimaxAdapt/SimulatedTrials/trial", i, ".csv", sep = ""))
    result <- with(trial, rbind(result, data.frame(trial = i, s1 = sum(resp[interim==1]),
                                                   s = sum(resp),
                                                   r1 = dsgn$r1,
                                                   r = dsgn$r,
                                                   n1 = length(resp[interim==1]),
                                                   n = length(resp), stop = stop[1])))
  }
  result$p0 <- dsgn$p0
  result$p1 <- dsgn$p1
  result$dsgnp1 <- dsgn$dsgnp1
  result$alpha <- dsgn$Alpha
  result$targetAlpha <- dsgn$targetAlpha
  result$beta <- dsgn$Beta
  result$targetBeta <- dsgn$targetBeta
  # Calculate new r
  result$newr2 <- NA
  result$newr <- result$r
  for (i in 1:nrow(result)){
    if (result$stop[i] == 0){
      ce <- switch(dsgn$restAlphaMet + 1,
                   getCE(dsgn, result$s1[i]),
                   getD_proportionally(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1],
                   getD_equally(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1],
                   getD_distributeToOne(dsgn, dsgn$targetAlpha)$ce[result$s1[i]+1])
      n2 <- as.numeric(result$n[i] - result$n1[i])
      result$newr2[i] <- get_r2_flex(ce, dsgn$p0, n2)
      result$newr[i] <- result$r1[i] + result$newr2[i]
    }

  }

  result$success <- as.numeric(result$s > result$newr)
  result$s2 <- result$s - result$s1
  result$n2 <- result$n - result$n1
  # Estimation of the response rate p
  result$pm1 <- result$s1/result$n1 #  Maximum likelihood (stage 1 only)
  result$pm2 <- NA
  result$pm2[result$n2 > 0] <- result$s2[result$n2 > 0]/result$n2[result$n2 > 0] #  Maximum likelihood (stage 2 only)
  result$pm <- result$s/result$n #  Maximum likelihood
  result$pg <- NA
  result$pu <- NA
  result$pp <- NA
  result$pk <- NA
  for (i in 1:nrow(result)){
    result$pg[i] <- pg(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pu[i] <- pu(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pp[i] <- pp(result$s[i], result$n1[i], result$r1[i], result$n[i])
    result$pk[i] <- pk(result$s[i], result$n1[i], result$r1[i], result$n[i], result$p0[i])
  }
  write.csv(result, "ResultsMinimaxDesignAdapt.csv", row.names = FALSE)
}

