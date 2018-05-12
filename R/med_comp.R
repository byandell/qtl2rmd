med_comp <- function(target, mediator, annotation, covar_tar, covar_med, kinship, intcovar, driver, driver_med) {
  # one driver
  med_test1 <- intermediate::mediation_test(
    target   = target,
    mediator = mediator,
    annotation = annotation,
    covar_tar = covar_tar,
    covar_med = covar_med,
    kinship = kinship,
    intcovar = intcovar,
    driver = driver,
    driver_med = NULL)
  (sum_med1 <- 
      dplyr::arrange(
        summary(med_test1),
        pvalue))
  
  med_testp <- intermediate::mediation_test(
    target   = target,
    mediator = mediator,
    annotation = annotation,
    covar_tar = covar_tar,
    covar_med = covar_med,
    kinship = kinship,
    intcovar = intcovar,
    driver = driver,
    driver_med = driver_med)
  (sum_medp <- 
      dplyr::arrange(
        summary(med_testp),
        pvalue))
  
  # two drivers assumed uncorrelated
  med_test2 <- intermediate::mediation_test(
    target   = target,
    mediator = mediator,
    annotation = annotation,
    covar_tar = covar_tar,
    covar_med = covar_med,
    kinship = kinship,
    intcovar = intcovar,
    driver = driver,
    driver_med = driver_med, frobenius = 1)
  (sum_med2 <- 
      dplyr::arrange(
        summary(med_test2),
        pvalue))

  # one driver at mediator
  med_testm <- intermediate::mediation_test(
    target   = target,
    mediator = mediator,
    annotation = annotation,
    covar_tar = covar_tar,
    covar_med = covar_med,
    kinship = kinship,
    intcovar = intcovar,
    driver = NULL,
    driver_med = driver_med)
  (sum_medm <- 
      dplyr::arrange(
        summary(med_testm),
        pvalue))
  
  out <- list(
    best = dplyr::bind_rows(
      target = sum_med1,
      two = sum_med2,
      depend = sum_medp,
      mediator = sum_medm,
      .id = "fit"),
    test = dplyr::bind_rows(
      target = med_test1$test,
      two = med_test2$test,
      depend = med_testp$test,
      mediator = med_testm$test,
      .id = "fit"))
  out$best <- dplyr::mutate(out$best, fit = factor(fit, c("mediator","target","two","depend")))
  out$test <- dplyr::mutate(out$test, fit = factor(fit, c("mediator","target","two","depend")))
  out
}
