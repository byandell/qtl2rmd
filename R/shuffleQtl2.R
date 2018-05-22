shuffleQtl2 <- function(driver_tar, target, kinship, covar) {
  tmp <- target * NA
  fit <- qtl2::fit1(driver_tar, target, kinship, covar)
  # Pay attention to missing values
  tmp[names(fit$fitted),] <- fit$fitted
  # Shuffle residuals from fit
  m <- which(!is.na(tmp[,]))
  ms <- sample(m, length(m))
  res <- (target - tmp)[ms,]
  tmp[m,] <- tmp[m,] + res

  tmp
}
shuffleQtl2M <- function(driver_med, mediator, kinship, covar, index) {
  for(i in seq_len(ncol(mediator))) {
    mediator[,i] <- shuffleQtl2(driver_med[,, index[i]], mediator[,i, drop = FALSE],
                                kinship, covar)
  }
  mediator
}
