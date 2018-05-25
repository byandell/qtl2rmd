mediator_quants <- function(med_tar_test, medfn = quantiles, ...){
  med_quants <-
    dplyr::mutate(
      dplyr::bind_rows(
        purrr::map(
          split(targets, targets$file),
          function(x) {
            med_tar_test <- readRDS(file.path(params$batch, x$file))
            out <- dplyr::bind_rows(
              purrr::map(
                split(
                  med_tar_test,
                  med_tar_test$id),
                medfn, ...),
              .id = "id")
            if(nrow(out)) {
              data.frame(
                target = x$target,
                chr = x$chr,
                out,
                stringsAsFactors = FALSE)
            } else {
              NULL
            }
          })),
      chr = factor(chr, c(1:19,"X")))
}

quantiles <- function(x, quants = c(.1,.2,.5,.8,.9), ...) {
  data.frame(pct = quants * 100,
             tar = quantile(x$pvalue.t, quants),
             med = quantile(x$pvalue.m, quants))
}
counts <- function(x, cutoff = 0.1, ...) {
  if(is.null(x))
    return(NULL)
  if(!all(c("triad.t","triad.m") %in% names(x)))
    return(NULL)
  tidyr::spread(
    dplyr::ungroup(
      dplyr::count(
        dplyr::group_by(
          dplyr::mutate(
            x,
            triad.t = ifelse(pvalue.t <= cutoff, triad.t, "n.s."),
            triad.m = ifelse(pvalue.m <= cutoff, triad.m, "n.s.")),
          triad.t, triad.m))),
    triad.m, n)
}
alts <- function(x, cutoff = 0.1, ...) {
  if(is.null(x))
    return(NULL)
  if(!all(c("triad.t","triad.m") %in% names(x)))
    return(NULL)
  tidyr::spread(
    dplyr::ungroup(
      dplyr::count(
        dplyr::group_by(
          dplyr::mutate(
            x,
            triad.t = ifelse(pvalue.t <= cutoff, triad.t, paste0(triad.t, "_ns"))),
          triad.t, alt.t))),
    alt.t, n)
}
logp <- function(x, quants = c(.1,.2,.5,.8,.9), ...) {
  tmpfn <- function(x, model_val) {
    x <- dplyr::mutate(
      dplyr::filter(
        x, 
        model == model_val),
      diffp = log10(pvalue.t) - log10(pvalue.m))
    data.frame(pct = quants * 100,
               diffp = quantile(x$diffp, quants),
               pvalue = quantile(x$pvalue.t, quants))
  }
  out <- dplyr::inner_join(
    tmpfn(x, "causal"),
    tmpfn(x, "independent"),
    suffix = c(".c",".i"),
    by = "pct")
  out2 <- 
    tidyr::spread(
      dplyr::select(
        dplyr::filter(
          x,
          model %in% c("causal","independent")),
        simnum, model, pvalue.t),
      model, pvalue.t)
  out2 <- data.frame(pct = quants * 100,
                     causal = quantile(out2$causal, quants),
                     independent = quantile(out2$independent, quants),
                     caus_indep = quantile(out2$causal + out2$independent, quants))
  dplyr::inner_join(out, out2, by = "pct")
}

summaryplot_quants <- function(med_quants, driver = "tar") {
  sums <- as.data.frame(t(as.data.frame(
    purrr::map(
      split(med_quants, med_quants$pct),
      function(x) c(summary(x[[driver]]))),
    check.names = FALSE)),
    check.names = FALSE)
  nsums <- names(sums)
  sums$quantile <- as.numeric(rownames(sums))
  sums <- sums %>%
    tidyr::gather("summary", "pvalue", -quantile) %>%
    dplyr::mutate(summary = factor(summary, rev(nsums)))
  ggplot(sums) +
    aes(pvalue, quantile, col = summary) +
    geom_abline(slope = 100, intercept = 0, col = "gray") +
    geom_line()
}
