peaks %>%
  select(pheno, chr, pos, lod) %>%
  filter(chr != "X") %>%
  group_by(pheno) %>%
  filter(lod == max(lod)) %>%
  arrange(desc(lod)) %>%
  as.data.frame %>%
  head(params$showPeaks) %>%
  print
