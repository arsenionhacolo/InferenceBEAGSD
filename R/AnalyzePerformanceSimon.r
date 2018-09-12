#' @title Performance of estimation methods
#' @description It takes the results produced by \code{\link{AnalyzeSimonDsgn}} and
#' \code{\link{AnalyzeSimonDsgnAdaptN}} and produces a dataframe containing bias,
#' mean square error and variance of the estimators. It also calculates the power
#' and the expected sample size (EN) where applicable.
#' @details Computations are done for different combinations of values of
#' \code{stop}, (0,1), and \code{success}, (0,1). See \code{\link{AnalyzeSimonDsgn}} or
#' \code{\link{AnalyzeSimonDsgnAdaptN}}. For instance, computations done on all
#' simulated trials are marked with \code{"both"} in the columns \code{stop} and
#' \code{success}, while the ones done only on trials that continued to the final
#' stage have \code{stop = "no"} and \code{success = "both"}.
#' @param designs Taking values \code{"fixed"}, \code{"adaptive"} or \code{"all"},
#' indicating whether only classical, adaptive or all designs should be included.
#' The default is \code{"all"}.
#' @param basedir The root directory in which simulations were performed. The current
#' working directory is assumed by default. It must contain all the files and folders
#' created by \code{\link{SimulateSimonDsgn}} and/or \code{\link{SimulateSimonDsgnAdaptN}}.
#' @return Dataframe containing bias, mean square error and variance of the estimators,
#' power, expected sample size, and design information.
#' @seealso \code{\link{AnalyzeSimonDsgn}}, \code{\link{AnalyzeSimonDsgnAdaptN}},
#' \code{\link{pdata}} and \code{\link{AnalyzePerformanceSimon2}}.
#' @export
#' @examples
#' \dontrun{
#' AnalyzePerformanceSimon()
#' }
#' @author Arsenio Nhacolo
AnalyzePerformanceSimon <- function(designs = "all", basedir = NA){
  if (!is.na(basedir)) setwd(basedir)
  #FIXED DESIGNS
  if (designs %in%  c("all","fixed")){
    # Optimal design
    rslt <- read.csv("ResultsOptimalDesign.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- pdata(t, "Optimal", "both", "both", nrep)
    t <- rslt[rslt$stop == 1,]
    presult <- rbind(presult, pdata(t, "Optimal", "yes", "no", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Optimal", "no", "both", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 1,]
    presult <- rbind(presult, pdata(t, "Optimal", "no", "yes", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 0,]
    presult <- rbind(presult, pdata(t, "Optimal", "no", "no", nrep))
    # Minimax design
    rslt <- read.csv("ResultsMinimaxDesign.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- rbind(presult, pdata(t, "Minimax", "both", "both", nrep))
    t <- rslt[rslt$stop == 1,]
    presult <- rbind(presult, pdata(t, "Minimax", "yes", "no", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Minimax", "no", "both", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 1,]
    presult <- rbind(presult, pdata(t, "Minimax", "no", "yes", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 0,]
    presult <- rbind(presult, pdata(t, "Minimax", "no", "no", nrep))

  }
  # ADAPTIVE DESIGNS
  if (designs %in%  c("all","adaptive")){
    # Adaptive optimal design
    rslt <- read.csv("ResultsOptimalDesignAdapt.csv")
    nrep <- nrow(rslt)
    t <- rslt
    if (!exists("presult")){
      presult <- pdata(t, "Adaptive optimal", "both", "both", nrep)
    }else{
      presult <- rbind(presult, pdata(t, "Adaptive optimal", "both", "both", nrep))
    }
    t <- rslt[rslt$stop == 1,]
    presult <- rbind(presult, pdata(t, "Adaptive optimal", "yes", "no", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive optimal", "no", "both", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 1,]
    presult <- rbind(presult, pdata(t, "Adaptive optimal", "no", "yes", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive optimal", "no", "no", nrep))

    # Adaptive minimax design
    rslt <- read.csv("ResultsMinimaxDesignAdapt.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "both", "both", nrep))
    t <- rslt[rslt$stop == 1,]
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "yes", "no", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "no", "both", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 1,]
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "no", "yes", nrep))
    t <- rslt[rslt$stop == 0 & rslt$success == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "no", "no", nrep))
  }
  write.csv(presult, file = "PerformanceResults.csv", row.names = F)
  presult
}
