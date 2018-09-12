
#' @title Performance, with respect to Phase III power, of phase II estimates.
#' @description \code{PerforIItoIIIRe} calculates the mean and median power in a Phase III
#' trials from the output of \code{\link{AnIItoIIIRe}}.
#' @param t Dataframe containing the output from the function \code{\link{AnIItoIIIRe}}.
#' @return Dataframe.
#' @references
#' Nhacolo, A. and Brannath, W. Using Estimates from Adaptive Phase II Oncology Trials to Plan Phase III Trials.
#' \emph{In press}.
#' @seealso \code{\link{AnIItoIIIRe}}.
#' @export
#' @examples
#' \dontrun{
#' rslt <- read.csv("ResultsAll.csv")
#' rsltfull <- rslt
#' rslt <- rslt[rslt$suco==1,]
#' rslt <- rslt[,c("pi0", "pi1", "spi1", "alpha", "beta", "suco", "pip",
#'                 "pim1", "pim2", "pim2v2", "pim3")]
#' rslt <- rslt[rslt$spi1>=rslt$pi0+0.1 & rslt$spi1<=rslt$pi1+0.3,]
#' rslt$spi1f <- factor(rslt$spi1)
#' cats <- levels(rslt$spi1f)
#' ncats <- length(cats)
#' setwd(paste0("C:/Users/arsenio/Documents/PhD/Simulations/Paper2/Reuse/pi01by0.01/50000/",did))
#' save(ncats,file = "ncats.rdata")
#' for (i in 1:ncats){
#'   sr <- rslt[rslt$spi1f==cats[i],]#Single result (result of a specific spi1)
#'   save(sr,file = paste0("sr",i,".rdata"))
#' }
#' load("ncats.rdata")
#' PerfAll <- data.frame()
#' for (k in 1:ncats){
#'  load(paste0("sr",k,".rdata"))
#'  sre <- AnIItoIIIRe(rslt = sr,f = c(.95,.96,.97,.98,.99))
#'   PerfAll <- rbind(PerfAll,PerforIItoIIIRe(sre))
#'   rm(sr)
#' }
#' write.csv(PerfAll, file = "PerfAllIItoIII.csv", row.names = F)
#' }
#' @author Arsenio Nhacolo
PerforIItoIIIRe <- function(t){
  pd <- data.frame(pi0 = t$pi0[1], pi1 = t$pi1[1], spi1 = t$spi1[1],
                   alpha = t$alpha[1], beta = t$beta[1])
  #Mean
  pd$meaPp <- mean(t$Pp)
  pd$meaPpf1 <- mean(t$Ppf1)
  pd$meaPpf2 <- mean(t$Ppf2)
  pd$meaPpf3 <- mean(t$Ppf3)
  pd$meaPpf4 <- mean(t$Ppf4)
  pd$meaPpf5 <- mean(t$Ppf5)

  pd$meaPm1 <- mean(t$Pm1)
  pd$meaPm1f1 <- mean(t$Pm1f1)
  pd$meaPm1f2 <- mean(t$Pm1f2)
  pd$meaPm1f3 <- mean(t$Pm1f3)
  pd$meaPm1f4 <- mean(t$Pm1f4)
  pd$meaPm1f5 <- mean(t$Pm1f5)

  pd$meaPm2 <- mean(t$Pm2)
  pd$meaPm2f1 <- mean(t$Pm2f1)
  pd$meaPm2f2 <- mean(t$Pm2f2)
  pd$meaPm2f3 <- mean(t$Pm2f3)
  pd$meaPm2f4 <- mean(t$Pm2f4)
  pd$meaPm2f5 <- mean(t$Pm2f5)

  pd$meaPm2v2 <- mean(t$Pm2v2)
  pd$meaPm2v2f1 <- mean(t$Pm2v2f1)
  pd$meaPm2v2f2 <- mean(t$Pm2v2f2)
  pd$meaPm2v2f3 <- mean(t$Pm2v2f3)
  pd$meaPm2v2f4 <- mean(t$Pm2v2f4)
  pd$meaPm2v2f5 <- mean(t$Pm2v2f5)

  pd$meaPm3 <- mean(t$Pm3)
  pd$meaPm3f1 <- mean(t$Pm3f1)
  pd$meaPm3f2 <- mean(t$Pm3f2)
  pd$meaPm3f3 <- mean(t$Pm3f3)
  pd$meaPm3f4 <- mean(t$Pm3f4)
  pd$meaPm3f5 <- mean(t$Pm3f5)

  #Median
  pd$medPp <- median(t$Pp)
  pd$medPpf1 <- median(t$Ppf1)
  pd$medPpf2 <- median(t$Ppf2)
  pd$medPpf3 <- median(t$Ppf3)
  pd$medPpf4 <- median(t$Ppf4)
  pd$medPpf5 <- median(t$Ppf5)

  pd$medPm1 <- median(t$Pm1)
  pd$medPm1f1 <- median(t$Pm1f1)
  pd$medPm1f2 <- median(t$Pm1f2)
  pd$medPm1f3 <- median(t$Pm1f3)
  pd$medPm1f4 <- median(t$Pm1f4)
  pd$medPm1f5 <- median(t$Pm1f5)

  pd$medPm2 <- median(t$Pm2)
  pd$medPm2f1 <- median(t$Pm2f1)
  pd$medPm2f2 <- median(t$Pm2f2)
  pd$medPm2f3 <- median(t$Pm2f3)
  pd$medPm2f4 <- median(t$Pm2f4)
  pd$medPm2f5 <- median(t$Pm2f5)

  pd$medPm2v2 <- median(t$Pm2v2)
  pd$medPm2v2f1 <- median(t$Pm2v2f1)
  pd$medPm2v2f2 <- median(t$Pm2v2f2)
  pd$medPm2v2f3 <- median(t$Pm2v2f3)
  pd$medPm2v2f4 <- median(t$Pm2v2f4)
  pd$medPm2v2f5 <- median(t$Pm2v2f5)

  pd$medPm3 <- median(t$Pm3)
  pd$medPm3f1 <- median(t$Pm3f1)
  pd$medPm3f2 <- median(t$Pm3f2)
  pd$medPm3f3 <- median(t$Pm3f3)
  pd$medPm3f4 <- median(t$Pm3f4)
  pd$medPm3f5 <- median(t$Pm3f5)

  return(pd)
}

