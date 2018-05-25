# must supply target_name, chr_id
args <- commandArgs(TRUE)
target_name <- args[[1]]
chr_id <- args[[2]]

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

tar_test <- extract_tables(batch = "sims",
                           target_name = target_name,
                           chr_id = chr_id, model = "tar", result = "test", driver = "snp")
med_test <- extract_tables(batch = "sims",
                           target_name = target_name,
                           chr_id = chr_id, model = "med", result = "test", driver = "snp")
med_tar_test <-
  dplyr::inner_join(
    tar_test,
    med_test,
    suffix = c(".t",".m"),
    by = c("simnum","id","model"))
saveRDS(med_tar_test, file.path("sims", target_name,
                                paste(target_name, chr_id, "med_tar_test.rds", sep = "_")))

tar_best <- extract_tables(batch = "sims",
                           target_name = target_name,
                           chr_id = chr_id, model = "tar", result = "best", driver = "snp")
med_best <- extract_tables(batch = "sims",
                           target_name = target_name,
                           chr_id = chr_id, model = "med", result = "best", driver = "snp")
med_tar_best <-
  dplyr::inner_join(
    tar_best,
    med_best,
    suffix = c(".t",".m"),
    by = c("simnum","id","symbol","chr"))
saveRDS(med_tar_best, file.path("sims", target_name,
                                paste(target_name, chr_id, "med_tar_best.rds", sep = "_")))
