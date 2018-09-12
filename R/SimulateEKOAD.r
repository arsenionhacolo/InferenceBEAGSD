#' @title Simulate single-arm binary endpoint two-stage adaptive designs.
#' @description \code{SimulateEKOAD} Simulate trials following designs similar to that of Englert and Kieser(2013)'s.
#' @details The original designs (like the ones in \code{\link{EKOptAdaptDesigns}}) must be pre-processed
#' using the function \code{\link{dsgnPrep}} to get extra information like the designs in \code{\link{EKOADwn}}.
#' @param replicates Number of trials to be simulated.
#' @param dsgn Dataframe containing one of the designs in \code{\link{EKOADwn}}.
#' @param newpi1 New response rate under the alternative hypothesis used to simulate trials. If \code{NULL} (default),
#' the one from the design is used.
#' @param seed The seed for random number generator. If \code{NULL} (default), no seed is set and , hence, results
#' are not reproducible.
#' @param deleteOld If \code{TRUE} (default), the simulation sub-directory is cleared before simulations start.
#' @return Simulated trials are saved in the sub-directory ./SimulatedTrials.
#' @references Englert, S. and Kieser, M. Optimal adaptive two-stage designs for phase {II} cancer clinical trials.
#' \emph{Biometrical Journal}, 2013.
#' @seealso \code{\link{EKOptAdaptDesigns}}, \code{\link{EKOADwn}}.
#' @export
#' @author Arsenio Nhacolo
SimulateEKOAD <- function(replicates, dsgn,
                          newpi1 = NULL, seed = NULL, deleteOld = TRUE){

  l1 <- dsgn$l1[1]
  u1 <- dsgn$u1[1]
  n1 <- dsgn$n1[1]
  pi0 <- dsgn$pi0[1]

  if (is.null(newpi1)) spi1 <- dsgn$pi1[1] else spi1 <- newpi1 #simulation pi1
  dsgn$spi1 <- spi1
  write.csv(dsgn, "Design.csv", row.names = FALSE)# Save design parameters
  if (!is.null(seed)) set.seed(seed)
  if (deleteOld && dir.exists("./SimulatedTrials")){
    unlink("./SimulatedTrials", recursive = TRUE)
  }
  dir.create("./SimulatedTrials", recursive = TRUE)
  cat("Simulating ", replicates, " trials... ", sep = "")
  for (i in 1:replicates){
    data1 <- data.frame(id = 1:n1, resp = sample(c(0,1), n1, TRUE, c(1-spi1,spi1)))
    data1$interim <- 1
    s1 <- sum(data1$resp)
    if (s1 <= l1 | s1 >= u1){
      data1$stop <- 1
      write.csv(data1, paste("./SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }else{
      n2 <- dsgn$n2[dsgn$x1 == s1]
      data2 <- data.frame(id = (n1+1):(n1+n2), resp = sample(c(0,1), n2,
                                                             TRUE, c(1-spi1,spi1)))
      data2$interim <- 2
      data <- rbind(data1, data2)
      data$stop <- 0
      write.csv(data, paste("./SimulatedTrials/trial", i, ".csv", sep = ""),
                row.names = FALSE)
    }
  }
  cat("Done. Trials saved in ./SimulatedTrials/\n", sep = "")
}
