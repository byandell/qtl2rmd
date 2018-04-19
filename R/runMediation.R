# Interactive setting of other parameters
if(!exists("doBatch"))
  doBatch <- FALSE
if(interactive() & !doBatch) {
  param_vals <- rmarkdown::knit_params_ask(
    RmdFilename,
    params = param_vals)
}

# Preset output format and filename.
# Only HTML enables plotly use.
out_format <- "html_document"

# Render document
rmarkdown::render(
  RmdFilename,
  output_format=out_format,
  output_file=out_filename,
  params=param_vals)                                           
