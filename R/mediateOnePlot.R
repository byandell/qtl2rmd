mediateOnePlot <- function(topsnps, type = c("pvalue","IC")) {
  type <- match.arg(type)
  switch(
    type,
    pvalue = {
      ggplot2::ggplot(best_snp) +
        ggplot2::aes(pos, -log10(pvalue), col = pattern, group = triad) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~ triad)
    },
    IC = {
      ggplot2::ggplot(best_snp) +
        ggplot2::aes(pos, IC, col = pattern, group = triad) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~ triad) +
        ggplot2::ylab("BIC on log10 scale")
    })
}