
#' @title Performance of estimation methods
#' @description It takes the results produced by \code{\link{AnalyzeSimonDsgn}} and
#' \code{\link{AnalyzeSimonDsgnAdaptN}} and produces a dataframe containing bias,
#' mean square error and variance of the estimators. It also calculates the power
#' and the expected sample size (EN) where applicable.
#' @details It is the same as \code{\link{AnalyzePerformanceSimon}}, but here the
#' estimation is done only for two sets: all trials (unconditional), and only trials
#' that continued to final stage (conditional).
#' @param designs Taking values \code{"fixed"}, \code{"adaptive"} or \code{"all"},
#' indicating whether only classical, adaptive or all designs should be included.
#' The default is \code{"all"}.
#' @param basedir The root directory in which simulations were performed. The current
#' working directory is assumed by default. It must contain all the files and folders
#' created by \code{\link{SimulateSimonDsgn}} and/or \code{\link{SimulateSimonDsgnAdaptN}}.
#' @return Dataframe containing bias, mean square error and variance of the estimators,
#' power, expected sample size, and design information.
#' @seealso \code{\link{AnalyzeSimonDsgn}}, \code{\link{AnalyzeSimonDsgnAdaptN}},
#' \code{\link{pdata}} and \code{\link{AnalyzePerformanceSimon}}.
#' @export
#' @examples
#' \dontrun{
#' AnalyzePerformanceSimon2()
#'
#' # Simulation example
#' seed = 1986
#' p0 <- 0.1
#' alpha <- 0.05
#' beta <- 0.1
#' repl <- 100 # number of replicated trials for each p
#' if (file.exists("PerforAll.csv")) unlink("PerforAll.csv")
#' coln <- TRUE
#' while (p0 < 0.5){
#'   pv <- seq(p0+0.2,p0+0.4,0.1) # p to simulate data
#'   p1v <- seq(p0+0.2,p0+0.3,0.1) # p to get design
#'   for (p1 in p1v){
#'     designParam <- CalculateSimonDsgn(p0, p1, alpha, beta)
#'     pstart <- p0+0.1
#'     SimulateSimonDsgn(repl, designParam, pstart, seed = seed)
#'     SimulateSimonDsgnAdaptN(repl, designParam, pstart, seed = seed)
#'     AnalyzeSimonDsgn()
#'     AnalyzeSimonDsgnAdaptN()
#'     perf <- AnalyzePerformanceSimon2()
#'     for (p in pv){
#'       SimulateSimonDsgn(repl, designParam, p, seed = seed)
#'       SimulateSimonDsgnAdaptN(repl, designParam, p, seed = seed)
#'       AnalyzeSimonDsgn()
#'       AnalyzeSimonDsgnAdaptN()
#'       perf <- rbind(perf, AnalyzePerformanceSimon2())
#'     }
#'     write.csv(perf, file = paste("PerforAll_a",alpha,"b",beta,"p0",p0,"p1",
#'                                  p1,".csv", sep = ""), row.names = F)
#'     write.table(perf, file ="PerforAll.csv", append = T, sep = ",", row.names = F, col.names = coln)
#'     coln <- FALSE
#'   }
#'   p0 <- p0+0.1
#' }
#' }
#' @author Arsenio Nhacolo
AnalyzePerformanceSimon2 <- function(designs = "all", basedir = NA){
  if (!is.na(basedir)) setwd(basedir)
  #FIXED DESIGNS
  if (designs %in%  c("all","fixed")){
    # Optimal design
    rslt <- read.csv("ResultsOptimalDesign.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- pdata(t, "Optimal", "both", "both", nrep)
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Optimal", "no", "both", nrep))

    # Minimax design
    rslt <- read.csv("ResultsMinimaxDesign.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- rbind(presult, pdata(t, "Minimax", "both", "both", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Minimax", "no", "both", nrep))
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
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive optimal", "no", "both", nrep))
    # Adaptive minimax design
    rslt <- read.csv("ResultsMinimaxDesignAdapt.csv")
    nrep <- nrow(rslt)
    t <- rslt
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "both", "both", nrep))
    t <- rslt[rslt$stop == 0,]
    presult <- rbind(presult, pdata(t, "Adaptive minimax", "no", "both", nrep))
  }
  write.csv(presult, file = "PerformanceResults.csv", row.names = F)
  presult
}

