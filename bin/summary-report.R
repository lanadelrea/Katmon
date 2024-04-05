install.packages("rmarkdown")
install.packages("knitr")
library(rmarkdown)
library(knitr)

args <- commandArgs(trailingOnly = TRUE)

# Define parameters
params <- list(lineage_assignment = args[1],
               bammix_plot = args[2],
               freyja_plot = args[3],
               aafplot_mut = args[4],
               aafplot_amp = args[5])

summary_report_rmd = args[6]
#sample = args[7]

# Render the R Markdown file with parameters
render(summary_report_rmd,
       output_format = 'pdf_document',
       params = params)