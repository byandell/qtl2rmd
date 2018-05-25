extract_tables <- function(batch = "sims",
                           target_name = "AA_G83_ins_secrete_gm_adj_IPI",
                           chr_id = "1", model = "tar", result = "test", driver = "snp") {
  dirname <- file.path(batch, target_name)
  filenames <- list.files(dirname, paste("", chr_id, model, result, driver, sep = "_"))
  index <- 
    as.numeric(
      stringr::str_extract(
        stringr::str_extract(
          filenames, "_snp_[0-9]+"), "[0-9]+"))
  filenames <- as.list(file.path(dirname, filenames))
  dplyr::bind_rows(
    purrr::map(
      filenames,
      function(x) 
        read.csv(x, stringsAsFactors = FALSE, row.names = 1)),
    .id = "simnum")
}
  