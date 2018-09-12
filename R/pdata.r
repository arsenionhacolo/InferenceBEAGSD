
#' @title Helper function for analyzing the performance of estimators
#' @description It takes the results produced by \code{\link{AnalyzeSimonDsgn}} or
#' \code{\link{AnalyzeSimonDsgnAdaptN}} and produces a dataframe containing bias,
#' mean square error and variance of the estimators.
#' @details It is a helper function for \code{\link{AnalyzePerformanceSimon}}.
#' It also calculates the power and the expected sample size (EN) where applicable.
#' @param t Dataframe containing results produced by \code{\link{AnalyzeSimonDsgn}}
#' or \code{\link{AnalyzeSimonDsgnAdaptN}}.
#' @param stop Taking value \code{"yes"}, \code{"no"} or \code{"both"}, indicating
#' that only trials that stopped, continued or both were analyzed.
#' @param success Taking value \code{"yes"}, \code{"no"} or \code{"both"}, indicating
#' that only trials that were successful, unsuccessful or both were analyzed.
#' @param replicates Number of trials analysed. It is equal to the number of
#' rows in \code{t}.
#' @return Dataframe containing bias, mean square error and variance of the estimators.
#' @seealso \code{\link{AnalyzeSimonDsgn}}, \code{\link{AnalyzeSimonDsgnAdaptN}} and
#' \code{\link{AnalyzePerformanceSimon}}.
#' @export
#' @examples
#' \dontrun{
#' rslt <- read.csv("ResultsOptimalDesign.csv")
#' nrep <- nrow(rslt)
#' t <- rslt
#' presult <- pdata(t, "Optimal", "both", "both", nrep)
#' t <- rslt[rslt$stop == 0,]
#' presult <- rbind(presult, pdata(t, "Optimal", "no", "both", nrep))
#' }
#' @author Arsenio Nhacolo
pdata <- function(t, design, stop, success, replicates){
  return(data.frame(design, p0 = t$p0[1], p1 = t$p1[1], dsgnp1 = t$dsgnp1[1],
                    targetAlpha = t$targetAlpha[1], targetBeta = t$targetBeta[1],
                    stop, success,
                    trials = nrow(t),
                    ptrials = 100*nrow(t)/replicates,
                    biaspm1 = mean(t$pm1 - t$p1),
                    msepm1 = mean((t$pm1 - t$p1)^2),
                    varpm1 = var(t$pm1),
                    #
                    biaspm2 = mean(t$pm2 - t$p1),
                    msepm2 = mean((t$pm2 - t$p1)^2),
                    varpm2 = var(t$pm2),
                    #
                    biaspm = mean(t$pm - t$p1),
                    msepm = mean((t$pm - t$p1)^2),
                    varpm = var(t$pm),
                    biaspg = mean(t$pg - t$p1),
                    msepg = mean((t$pg - t$p1)^2),
                    varpg = var(t$pg),
                    biaspu = mean(t$pu - t$p1),
                    msepu = mean((t$pu - t$p1)^2),
                    varpu = var(t$pu),
                    biaspp = mean(t$pp - t$p1),
                    msepp = mean((t$pp - t$p1)^2),
                    varpp = var(t$pp),
                    biaspk = mean(t$pk - t$p1),
                    msepk = mean((t$pk - t$p1)^2),
                    varpk = var(t$pk),
                    EN = ifelse(stop == "both" & success == "both", mean(t$n), NA),
                    power = ifelse(stop == "both" & success == "both", mean(t$success), NA)))
}

