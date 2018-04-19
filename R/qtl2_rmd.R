#' R/qtl2 analysis pipeline for one locus
#'
#' Run analysis for multiple traits at one locus.
#'
#' @param rmd_name base name of Rmd file
#' @param datapath path to data files
#' @param resultpath path to result files
#' @param output_name name of output file to be created in \code{resultpath} (without MIME ending)
#' @param ... named list of parameters to pass to rmarkdown file
#'
#' @return phename phenotype name = directory name for results
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dontrun(qtl2_rmd(rmd_name, ...)
#'
#' @export
qtl2_rmd <- function(rmd_name,
                     datapath,
                     resultpath,
                     output_name,
                     ...) {
  
  if(is.null(rmd_name))
    stop("must supply rmd_name")
  if(is.null(output_name))
    stop("must supply output_name")
  param_vals <- list(datapath=datapath,
                     resultpath=resultpath,
                     ...)
  if(interactive())
    param_vals <- rmarkdown::knit_params_ask(system.file(rmd_name, package="doqtl2"),
                                             params = param_vals)
  if(is.null(param_vals))
    stop(paste("render of", rmd_name, "cancelled"))
  
  set_output <- param_vals$set_output
  if(is.null(set_output))
    set_output <- "html"
  switch(set_output,
         pdf = {
           out_format <- "pdf_document"
           out_filename <- file.path(resultpath, paste0(output_name, ".pdf"))
         },
         word = {
           out_format <- "word_document"
           out_filename <- file.path(resultpath, paste0(output_name, ".docx"))
         },
         {
           out_format <- "html_document"
           out_filename <- file.path(resultpath, paste0(output_name, ".html"))
         })
  rmarkdown::render(system.file(rmd_name, package="doqtl2"),
                    output_format=out_format,
                    output_file=out_filename,
                    params=param_vals)
}
