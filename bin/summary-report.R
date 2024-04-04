install.packages("rmarkdown")
library(rmarkdown)

args <- commandArgs(trailingOnly = TRUE)

# Define parameters
params <- list(bammix_plot = args[1],
               freyja_plot = args[2],
               aafplot_mut = args[3],
               aafplot_amp = args[4])

summary_report_rmd = args[5]

# Render the R Markdown file with parameters
render(summary_report_rmd,
       output_format = 'pdf_document',
       params = params)