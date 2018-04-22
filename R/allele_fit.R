allele_fit <- function(out) {
  out1 <- 
    dplyr::arrange(
      dplyr::ungroup(
        dplyr::mutate(
          dplyr::group_by(
            dplyr::mutate(
              tidyr::gather(
                dplyr::select(
                  out,
                  pheno, QTL, symbol, pvalue, A_m:H_m, A:H, A_p:H_p),
                allele, effect, -pheno, -QTL, -symbol, -pvalue),
              fitType = stringr::str_remove(allele, "[A-H]_*"),
              allele = stringr::str_remove(allele, "_[a-z]"),
              fitType = ifelse(fitType == "", "a", fitType),
              fitType = c(m="Mediator", p="Target", a="Adjusted")[fitType],
              fitType = factor(fitType, c("Mediator", "Target", "Adjusted"))),
            pheno, QTL, symbol, fitType),
          effect = effect - mean(effect))),
      pvalue)
  class(out1) <- c("allele_fit", class(out1))
  out1
}
ggplot_allele_fit <- function(out1,
                              colors = qtl2::CCcolors) {
  allele_names <- names(colors)
  names(allele_names) <- sort(unique(out1$allele))
  out1 <- 
    dplyr::mutate(
      dplyr::filter(
        dplyr::mutate(
          out1,
          symbol = paste0(symbol, " (", signif(pvalue, 2), ")")),
        symbol %in% unique(symbol)[1:12]),
      allele = allele_names[allele],
      symbol = factor(symbol, unique(symbol)))
  ggplot2::ggplot(out1) +
    ggplot2::aes(fitType, effect, color = allele, group = allele) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size=2) +
    ggplot2::facet_wrap(~ symbol) +
    ggplot2::ggtitle(paste(out1$pheno[1], "for QTL", out1$QTL[1])) +
    ggplot2::scale_color_manual(name = "allele",
                                values = colors)
  
}
autoplot.allele_fit <- function(...) ggplot_allele_fit(...)
