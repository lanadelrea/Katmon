library(rmarkdown)
library(knitr)
library(kableExtra)
library(readr)


args <- commandArgs(trailingOnly = TRUE)


# Input table and plots

lineage_table_tsv <- args[1]
lineage_table <- read_tsv(lineage_table_tsv)

freyja_plot <- args[2]

virstrain_table_tsv <- args[3]
virstrain_table <- read_tsv(virstrain_table_tsv)
summary_report_rmd = args[4]


# Defining parameters
params <- list(lineage_assignment = lineage_table,
               virstrain = virstrain_table,
               freyja_plot = freyja_plot)


# Render the R Markdown file with parameters
render(summary_report_rmd,
       output_format = 'html_document',
       params = params)