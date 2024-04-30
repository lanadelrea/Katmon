library(rmarkdown)
library(knitr)
library(kableExtra)
library(readr)


args <- commandArgs(trailingOnly = TRUE)

lineage_table_tsv <- args[1]
lineage_table <- read_tsv(lineage_table_tsv)

# Define parameters
params <- list(lineage_assignment = lineage_table,
               bammix_plot = args[2],
               freyja_plot = args[3],
               aafplot_mut = args[4],
               aafplot_amp = args[5])

summary_report_rmd = args[6]
#sample = args[7]

# Render the R Markdown file with parameters
render(summary_report_rmd,
       output_format = 'html_document',
       params = params)