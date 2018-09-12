#' @title Simon's designs data simulation
#' @description \code{SimulateSimonDsgn} simulates data from Simon's optimal and minimax designs.
#' @details The simulated trials are stored in the sub-directories \code{/Optimal/SimulatedTrials}
#' and \code{/Minimax/SimulatedTrials} for optimal and minimax designs, repectively, under the
#' current working directory. The sub-directories are automatically created. Individual trial
#' data are stored in a CSV file named \code{trial#}, where # is the replicate number.
#' @param replicates Number of trials to be generated.
#' @param designParam A dataframe containing Simon's optimal and minimax designs, as returned
#' by the function \code{\link{CalculateSimonDsgn}}.
#' @param newp1 If \code{NA} (default) data are generated assuming the same response probability
#' under alternative hypothesis, \code{p1}, used to get the designs
#' (see \code{\link{CalculateSimonDsgn}}). One may provide different values of \code{newp1} if
#' there is interest in studying the effect of departure from the design's assumed \code{p1}.
#' @param seed Initial value (any integer) of random-number seed. It  is useful for creating
#' simulations that can be reproduced. The default is \code{NA}, meaning no reproducibility.
#' @param deleteOld If TRUE (default) the sub-directories \code{/Optimal/SimulatedTrials} and
#' \code{/Minimax/SimulatedTrials} are deleted, if they exist, before simulation starts. The
#' old data files are still replaced by the new ones even if \code{deleteOld} is set to
#' \code{FALSE}, but some old files remain in cases where the previous \code{replicates} was
#' greater that the current one.
#' @return The function is not intended to return an R object, instead it creates files
#' (in CSV format) containing simulated trials data. See \emph{Details}. It also saves
#' in the current working directory the \code{designParam} argument (\emph{DesignParameters.csv}).
#' @seealso \code{\link{CalculateSimonDsgn}}, \code{\link{SimulateSimonDsgnAdaptN}}
#' and \code{\link{AnalyzeSimonDsgn}}.
#' @export
#' @examples
#' d <- CalculateSimonDsgn(0.2, 0.4, 0.05, 0.1)
#' SimulateSimonDsgn(100, d, seed = 1986)
#' @author Arsenio Nhacolo
SimulateSimonDsgn <- function(replicates, designParam, newp1 = NA, seed = NA, deleteOld = TRUE){
  designParam$dsgnp1 <- designParam$p1
  if (!is.na(newp1)) designParam$p1 <- newp1
  # Save design parameters
  write.csv(designParam, "DesignParameters.csv", row.names = FALSE)
  # Optimal design trials
  if (!is.na(seed)) set.seed(seed)
  if (deleteOld && dir.exists("./Optimal/SimulatedTrials")){
    unlink("./Optimal/SimulatedTrials", recursive = TRUE)
  }
  dir.create("./Optimal/SimulatedTrials", recursive = TRUE)
  dsgn <- designParam[designParam$Type == "Optimal",]
  for (i in 1:replicates){
    data1 <- data.frame(id = 1:dsgn$n1,
                        resp = sample(c(0,1), dsgn$n1, TRUE,
                                      c(1-dsgn$p1,dsgn$p1)))
    data1$interim <- 1
    if (sum(data1$resp) <= dsgn$r1){
      data1$stop <- 1
      write.csv(data1, paste("./Optimal/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }else{
      data2 <- data.frame(id = (dsgn$n1+1):dsgn$n,
                          resp = sample(c(0,1), dsgn$n - dsgn$n1,
                                        TRUE, c(1-dsgn$p1,dsgn$p1)))
      data2$interim <- 2
      data <- rbind(data1, data2)
      data$stop <- 0
      write.csv(data, paste("./Optimal/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }
  }
  # Minimax design trials
  if (!is.na(seed)) set.seed(seed)
  if (deleteOld && dir.exists("./Minimax/SimulatedTrials")){
    unlink("./Minimax/SimulatedTrials", recursive = TRUE)
  }
  dir.create("./Minimax/SimulatedTrials", recursive = TRUE)
  dsgn <- designParam[designParam$Type == "MiniMax",]
  for (i in 1:replicates){
    data1 <- data.frame(id = 1:dsgn$n1,
                        resp = sample(c(0,1), dsgn$n1, TRUE,
                                      c(1-dsgn$p1,dsgn$p1)))
    data1$interim <- 1
    if (sum(data1$resp) <= dsgn$r){
      data1$stop <- 1
      write.csv(data1, paste("./Minimax/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }else{
      data2 <- data.frame(id = (dsgn$n1+1):dsgn$n,
                          resp = sample(c(0,1), dsgn$n - dsgn$n1,
                                        TRUE, c(1-dsgn$p1,dsgn$p1)))
      data2$interim <- 2
      data <- rbind(data1, data2)
      data$stop <- 0
      write.csv(data, paste("./Minimax/SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }
  }
}
