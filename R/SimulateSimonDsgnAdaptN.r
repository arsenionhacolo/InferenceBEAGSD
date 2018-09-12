
# Function to simulate data from adaptive Simon's optimal and minimax phase II designs
#' @title Simon's adaptive designs data simulation
#' @description Simulates data from adaptive versions of Simon's optimal and minimax designs,
#' propposed by \emph{Englert and Kieser (2012)}. Adaptation consists in recalculating the
#' second stage sample size n2 in order to achieve a desired conditional power given the
#' number of successes at first stage.
#' @details The simulated trials are stored in the sub-directories \code{/OptimalAdapt/SimulatedTrials}
#' and \code{/MinimaxAdapt/SimulatedTrials} for optimal and minimax designs, repectively, under the
#' current working directory. The sub-directories are automatically created. Individual trial
#' data are stored in a CSV file named \code{trial#}, where # is the replicate number.
#' @param replicates Number of trials to be generated.
#' @param designParam A dataframe containing Simon's optimal and minimax designs, as returned
#' by the function \code{\link{CalculateSimonDsgn}}.
#' @param newp1 If \code{NA} (default) data are generated assuming the same response probability
#' under alternative hypothesis, \code{p1}, used to get the designs
#' (see \code{\link{CalculateSimonDsgn}}). One may provide different values of \code{newp1} if
#' there is interest in studying the effect of departure from the design's assumed \code{p1}.
#' @param condPwr The desired conditional power. The default is \code{1-beta}.
#' @param restAlphaMet The method for spending the "rest alpha" (difference between nominal alpha
#' level and actual alpha level for the given design).
#' \itemize{
#'  \item 0: "rest alpha" is not used (default);
#'  \item 1: "rest alpha" is spent proportionally;
#'  \item 2: "rest alpha" is spent equally;
#'  \item 3: "rest alpha" is spent only to the worst case scenario (minimal number of responses at
#'            the interim analysis so that the study can proceed to the second stage).
#' }
#' @param seed Initial value (any integer) of random-number seed. It  is useful for creating
#' simulations that can be reproduced. The default is \code{NA}, meaning no reproducibility.
#' @param deleteOld If TRUE (default) the sub-directories \code{/OptimalAdapt/SimulatedTrials} and
#' \code{/MinimaxAdapt/SimulatedTrials} are deleted, if they exist, before simulation starts. The
#' old data files are still replaced by the new ones even if \code{deleteOld} is set to
#' \code{FALSE}, but some old files remain in cases where the previous \code{replicates} was
#' greater that the current one.
#' @return The function is not intended to return an R object, instead it creates files
#' (in CSV format) containing simulated trials data. See \emph{Details}. It also saves
#' in the current working directory the \code{designParam} argument (\emph{DesignParametersAdapt.csv}).
#' @references Englert S., Kieser M. Adaptive designs for single-arm phase II trials in oncology.
#' \emph{Pharm Stat}, 2012, 11, 241-249.
#' @seealso \code{\link{CalculateSimonDsgn}}, \code{\link[OneArmPhaseTwoStudy]{getN2}},
#' \code{\link{SimulateSimonDsgn}} and \code{\link{AnalyzeSimonDsgnAdaptN}}.
#' @export
#' @examples
#' d <- CalculateSimonDsgn(0.2, 0.4, 0.05, 0.1)
#' SimulateSimonDsgnAdaptN(100, d, seed = 1986)
#' @author Arsenio Nhacolo
SimulateSimonDsgnAdaptN <- function(replicates, designParam, newp1 = NA, condPwr = NA,
                                    restAlphaMet = 0, seed = NA, deleteOld = TRUE){
  designParam$dsgnp1 <- designParam$p1
  if (!is.na(newp1)) designParam$p1 <- newp1
  if (is.na(condPwr)) condPwr <- 1 - designParam$targetBeta
  designParam$condPwr <- condPwr
  designParam$restAlphaMet <- restAlphaMet
  # Save design parameters
  write.csv(designParam, "DesignParametersAdapt.csv", row.names = FALSE)
  # Optimal design trials
  dsgn <- designParam[designParam$Type == "Optimal",]
  if (!is.na(seed)) set.seed(seed)
  if (deleteOld && dir.exists("./OptimalAdapt/SimulatedTrials")){
    unlink("./OptimalAdapt/SimulatedTrials", recursive = TRUE)
  }
  dir.create("./OptimalAdapt/SimulatedTrials", recursive = TRUE)
  for (i in 1:replicates){
    data1 <- data.frame(id = 1:dsgn$n1,
                        resp = sample(c(0,1), dsgn$n1, TRUE,
                                      c(1-dsgn$p1,dsgn$p1)))
    data1$interim <- 1
    s1 <- sum(data1$resp)
    if (s1 <= dsgn$r1){
      data1$stop <- 1
      write.csv(data1, paste("./OptimalAdapt/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }else{
      n2adpt <- getN2v2(dsgn$condPwr, dsgn$p1, dsgn, s1, restAlphaMet, dsgn$targetAlpha)
      data2 <- data.frame(id = (dsgn$n1+1):(dsgn$n1+n2adpt),
                          resp = sample(c(0,1), n2adpt, TRUE, c(1-dsgn$p1,dsgn$p1)))
      data2$interim <- 2
      data <- rbind(data1, data2)
      data$stop <- 0
      write.csv(data, paste("./OptimalAdapt/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }
  }
  # Minimax design trials
  dsgn <- designParam[designParam$Type == "MiniMax",]
  if (!is.na(seed)) set.seed(seed)
  if (deleteOld && dir.exists("./MinimaxAdapt/SimulatedTrials")){
    unlink("./MinimaxAdapt/SimulatedTrials", recursive = TRUE)
  }
  dir.create("./MinimaxAdapt/SimulatedTrials", recursive = TRUE)
  for (i in 1:replicates){
    data1 <- data.frame(id = 1:dsgn$n1,
                        resp = sample(c(0,1), dsgn$n1, TRUE,
                                      c(1-dsgn$p1,dsgn$p1)))
    data1$interim <- 1
    s1 <- sum(data1$resp)
    if (s1 <= dsgn$r1){
      data1$stop <- 1
      write.csv(data1, paste("./MinimaxAdapt/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }else{
      n2adpt <- getN2v2(dsgn$condPwr, dsgn$p1, dsgn, s1, restAlphaMet, dsgn$targetAlpha)
      data2 <- data.frame(id = (dsgn$n1+1):(dsgn$n1+n2adpt),
                          resp = sample(c(0,1), n2adpt, TRUE, c(1-dsgn$p1,dsgn$p1)))
      data2$interim <- 2
      data <- rbind(data1, data2)
      data$stop <- 0
      write.csv(data, paste("./MinimaxAdapt/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }
  }
}

