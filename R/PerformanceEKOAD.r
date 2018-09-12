
#' @title Performance of estimation methods
#' @description \code{PerformanceEKOAD} calculates performance measures (bias, mean
#' square error, coverage probability) of the estimation methods based on the
#' results produced by \code{\link{AnalyzeEKOAD}}.
#' @param basedir The base directory containing the file with the results (\code{Results.csv}).
#'  If \code{NULL} (default), the current working directory is uded.
#' @return A dataframe with the performance results. A copy is saved in the file the
#'  \code{PerformanceResults.csv} in the \code{basedir}.
#' @references Nhacolo, A. and Brannath, W. Interval and point estimation in adaptive Phase II trials with binary endpoint.
#' \emph{Stat Methods Med Res}, 2018.
#' @seealso \code{\link{SimulateEKOAD}}, \code{\link{AnalyzeEKOAD}}.
#' @export
#' @examples
#' \dontrun{
#' #SIMULATIONS
#' for (did in c(6,10)){#Design ID
#' cat("============================ Design ",did," ============================\n")
#' repl <- 50000 # number of replicated trials for each p
#' dir.create(as.character(did))
#' setwd(as.character(did))
#' design <- EKOADwn[EKOADwn$id==did,]
#' seed = 3343
#' if (file.exists("PerforAll.csv")) unlink("PerforAll.csv")
#' piv <- seq(0,1,0.025) # p to simulate data
#' resul <- data.frame()
#' perf <- data.frame()
#' k <- 0
#' pl <- length(piv)
#' for (pi in piv){
#'   k <- k+1
#'   cat("_________________________ pi = ",pi, " (",k," of ",pl,") _________________________\n",sep = "")
#'   SimulateEKOAD(replicates = repl, dsgn = design, newpi1 = pi, seed = seed)
#'   resul <- rbind(resul, AnalyzeEKOAD())
#'   perf <- rbind(perf, PerformanceEKOAD())
#' }
#' write.table(resul, file ="ResultsAll.csv", sep = ",", row.names = F, col.names = TRUE)
#' write.table(perf, file ="PerforAll.csv", sep = ",", row.names = F, col.names = TRUE)
#' cat("Design ID: ", design$id[1], "\nReplicates: ", repl, "\nSeed: ", seed,
#'     "\nDate last run: ", date(),file = "info.txt", sep = "", append = FALSE)
#' }
#' }
#' @author Arsenio Nhacolo
PerformanceEKOAD <- function(basedir = NULL){
  if (!is.null(basedir)) setwd(basedir)
  cat("Analyzing performance of estimation methods... ", sep = "")
  rslt <- read.csv("Results.csv")
  nrep <- nrow(rslt)
  t <- rslt
  presult <- pdata2(t, "both", "both", nrep)
  #t <- rslt[rslt$stop == 0,]
  #presult <- rbind(presult, pdata(t, "no", "both", nrep))

  write.csv(presult, file = "PerformanceResults.csv", row.names = F)
  cat("Done. \nPerformance results saved in PerformanceResults.csv\n", sep = "")
  return(presult)
}
