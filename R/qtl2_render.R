#' Run Rmd file through purl and execute
#'
#' Use `purl` to extract R code from Rmd and source it.
#'
#' @param rmd_name base name of Rmd file
#' @param package package containing \code{rmd_name}
#' @param envir environment for temporary file
#'
#' @return nothing, but evaluates code in \code{rmd_name}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{rmd_run(rmd_name, ...)}
#'
#' @export
rmd_run <- function(rmd_name, package="doqtl2", envir=globalenv()){
  rmd_name <- system.file(rmd_name, package=package)
  tempR <- tempfile(tmpdir = ".", fileext = ".R")
  on.exit(unlink(tempR))
  knitr::purl(rmd_name, output=tempR)
  sys.source(tempR, envir=envir)
}


#' R/qtl2 analysis pipeline for one phenotype
#'
#' Run genome scan then SNP analysis for each significant peak.
#'
#' @param phename_filter character string to select phenotype (partial match)
#' @param phename_drop space-separated character string of partial matches to drop
#' @param datapath character string path to directory with data
#' @param resultpath character string path to directory with results
#' @param window_Mbp size of window around peak_Mbp (in Mbp)
#' @param extend_bp how many bp to extend SNP search beyond feature (default 5000)
#' @param analysis_type which analysis type to use as character string, typically "anal1" or "anal2" (default "")
#' @param lod_min minimum lod peak to process (default 5.5)
#'
#' @return phename phenotype name = directory name for results
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(qtl2_onepheno(phename, pheno_set, dirpath, outfile)
#'
#' @export
qtl2_onepheno <- function(phename_filter, phename_drop="",
                          datapath="",
                          resultpath="~",
                          window_Mbp = 3, # Mbp window around peak
                          extend_bp = 5000, # base pairs beyond genes
                          analysis_type = "", # analysis setup if more than one done
                          lod_min = 5.5) { # minimum lod peak

  ## Create output directory if it does not exist.
  if(!dir.exists(file.path(resultpath)))
    dir.create(file.path(resultpath))

  ## Do genome scan to find major peaks.
  ## For now, this uses file.path(dirpath, "peaks.csv") as guide.
  outtype <- "html"
  qtl2_render(phename_filter, phename_drop,
              datapath=datapath, resultpath=resultpath,
              analysis_type=analysis_type,
              outfile=file.path(resultpath,
                                paste("Scan", outtype, sep=".")))

  ## Read peaks from scan.
  peaks <- read_csv(file.path(resultpath, "peaks.csv"))

  ## Run analysis for each peak.
  for(peak in seq_len(nrow(peaks))) {
    if(peaks$lod[peak] >= lod_min)
      qtl2_render_snp(phename_filter, phename_drop,
                      peaks$chr[peak], peaks$pos[peak],
                      window_Mbp, extend_bp,
                      datapath=datapath, resultpath=resultpath,
                      analysis_type=analysis_type,
                      outfile=file.path(resultpath,
                                        paste0("SNP_",
                                               peaks$chr[peak], ".",
                                               outtype)))
  }
  file.path(resultpath)
}

#' Render R/qtl2 analysis
#'
#' Use R/rmarkdown::render() to customize Rmd for trait.
#'
#' @param phename_filter character string to select phenotype (partial match)
#' @param phename_drop space-separated character string of partial matches to drop
#' @param datapath character string path to directory with data
#' @param resultpath character string path to directory with results
#' @param analysis_type which analysis type to use as character string, typically "anal1" or "anal2" (default "")
#' @param outfile character string of output HTML file
#'
#' @return
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(qtl2_render(phename, pheno_set, dirpath, outfile)
#'
#' @export
qtl2_render <- function(phename_filter = "",
                        pheno_drop = "",
                        datapath = NULL,
                        resultpath = datapath,
                        analysis_type = "",
                        outfile = NULL) {
  if(is.null(datapath))
    stop("must supply datapath")
  if(is.null(outfile))
    stop("must supply outfile")

  param_vals <- list(phename_filter=phename_filter,
                 phename_drop=phename_drop,
                 datapath=datapath,
                 resultpath=resultpath,
                 analysis_type=analysis_type)
  rmarkdown::render(system.file("doc/analysis_onepheno.Rmd",
                                package="doqtl2"),
                    output_file=outfile,
                    params=param_vals)
}

#' Render R/qtl2 SNP study
#'
#' Use R/rmarkdown::render() to customize Rmd for SNP study.
#'
#' @param phename_filter character string to select phenotype (partial match)
#' @param phename_drop space-separated character string of partial matches to drop
#' @param chr chromosome name
#' @param peak_Mbp chromosome position in Mbp
#' @param window_Mbp size of window around peak_Mbp (in Mbp)
#' @param extend_bp how many bp to extend SNP search beyond feature (default 5000)
#' @param datapath character string path to directory with data
#' @param resultpath character string path to directory with results
#' @param analysis_type which analysis type to use as character string, typically "anal1" or "anal2" (default "")
#' @param outfile character string of output HTML file
#'
#' @return
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(qtl2_render_snp(phename, pheno_set, dirpath, outfile)
#'
#' @export
qtl2_render_snp <- function(phename_filter,
                            phename_drop="",
                            chr = "9",
                            peak_Mbp = 105.5,
                            window_Mbp = 3, # Mbp window around peak
                            extend_bp = 5000, # base pairs beyond genes
                            datapath = file.path("~/Documents/Research/attie_alan/DO", "data"),
                            resultpath = datapath,
                            analysis_type = "",
                            outfile = NULL) {
  if(is.null(outfile))
    stop("must supply outfile")

  param_vals <- list(phename_filter=phename_filter,
                 phename_drop=phename_drop,
                 chr_id=as.character(chr),
                 peak_Mbp=peak_Mbp,
                 window_Mbp=window_Mbp,
                 extend_bp=extend_bp,
                 datapath=datapath,
                 resultpath=resultpath,
                 analysis_type=analysis_type)
  rmarkdown::render(system.file("doc/analysis_onepeak.Rmd",
                                package="doqtl2"),
                    output_file=outfile,
                    params=param_vals)
}

#' Size of data files
#'
#' Find size and types of data files.
#'
#' @param datapath character string path to folder with files
#'
#' @return tbl summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(qtl2_datafiles("."))
#'
#' @export
qtl2_datafiles <- function(datapath) {
  datafiles <- system(paste("ls -sh", datapath), intern=TRUE)
  datafiles <- lapply(strsplit(datafiles, " ", fixed=TRUE),
                      function(x) x[x!=""])
  datafiles[[1]] <- NULL
  datafiles <- data.frame(size=sapply(datafiles, function(x) as.numeric(x[1])),
                          name=sapply(datafiles, function(x) x[2]))
  datafiles <- datafiles %>%
    filter(size>0) %>%
    mutate(name=as.character(name)) %>%
    mutate(type=sapply(strsplit(name, ".", fixed=TRUE),
                       function(x)x[length(x)])) %>%
    select(name,size,type)
}
