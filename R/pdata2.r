
#' @title Helper function for analyzing the performance of estimators
#' @description \code{pdata2} is a helper function used by function \code{\link{PerformanceEKOAD}}.
#' @details Not to be used directly.
#' @return Dataframe
#' @export
#' @author Arsenio Nhacolo
pdata2 <- function(t, stop, success, replicates){
  pd <- data.frame(pi0 = t$pi0[1], pi1 = t$pi1[1], spi1 = t$spi1[1],
                   alpha = t$alpha[1], beta = t$beta[1],
                   stop, success,
                   trials = nrow(t),
                   ptrials = 100*nrow(t)/replicates,

                   meanbpip1 = mean(t$pip1 - t$spi1),
                   medianbpip1 = median(t$pip1 - t$spi1),
                   msepip1 = mean((t$pip1 - t$spi1)^2),
                   varpip1 = var(t$pip1),

                   meanbpip = mean(t$pip - t$spi1),
                   medianbpip = median(t$pip - t$spi1),
                   msepip = mean((t$pip - t$spi1)^2),
                   varpip = var(t$pip),
                   EN = ifelse(stop == "both" & success == "both", mean(t$n), NA),
                   powero = ifelse(stop == "both" & success == "both", mean(t$suco), NA))

  pd$meanbpim1 = mean(t$pim1 - t$spi1)
  pd$medianbpim1 = median(t$pim1 - t$spi1)
  pd$msepim1 = mean((t$pim1 - t$spi1)^2)
  pd$varpim1 = var(t$pim1)
  pd$powerm1 = ifelse(stop == "both" & success == "both", mean(t$sucm1), NA)
  pd$covpm1 = mean(t$covm1) # coverage probability

  pd$meanbpim2 = mean(t$pim2 - t$spi1)
  pd$medianbpim2 = median(t$pim2 - t$spi1)
  pd$msepim2 = mean((t$pim2 - t$spi1)^2)
  pd$varpim2 = var(t$pim2)
  pd$powerm2 = ifelse(stop == "both" & success == "both", mean(t$sucm2), NA)
  pd$covpm2 = mean(t$covm2) # coverage probability

  pd$meanbpim2v2 = mean(t$pim2v2 - t$spi1)
  pd$medianbpim2v2 = median(t$pim2v2 - t$spi1)
  pd$msepim2v2 = mean((t$pim2v2 - t$spi1)^2)
  pd$varpim2v2 = var(t$pim2v2)
  pd$powerm2v2 = ifelse(stop == "both" & success == "both", mean(t$sucm2v2), NA)
  pd$covpm2v2 = mean(t$covm2v2) # coverage probability

  pd$meanbpim3 = mean(t$pim3 - t$spi1)
  pd$medianbpim3 = median(t$pim3 - t$spi1)
  pd$msepim3 = mean((t$pim3 - t$spi1)^2)
  pd$varpim3 = var(t$pim3)
  pd$powerm3 = ifelse(stop == "both" & success == "both", mean(t$sucm3), NA)
  pd$covpm3 = mean(t$covm3) # coverage probability

  return(pd)
}
