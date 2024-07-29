library(rmarkdown)
library(knitr)
library(kableExtra)
library(readr)


args <- commandArgs(trailingOnly = TRUE)


# Input table and plots

lineage_table_tsv <- args[1]
lineage_table <- read_tsv(lineage_table_tsv)
bammix_plots_str <- args[2]
freyja_plot <- args[3]
aafplots_mut_str <- args[4]
aafplots_amp_str <- args[5]
virstrain_table_tsv <- args[6]
virstrain_table <- read_tsv(virstrain_table_tsv)
summary_report_rmd = args[7]


bammix_plots <- strsplit(bammix_plots_str, ",")[[1]]
aafplots_mut <- strsplit(aafplots_mut_str, ",")[[1]]
aafplots_amp <- strsplit(aafplots_amp_str, ",")[[1]]


# Defining parameters
params <- list(lineage_assignment = lineage_table,
               virstrain = virstrain_table,
               bammix_plots = bammix_plots,
               freyja_plot = freyja_plot,
               aafplots_mut = aafplots_mut,
               aafplots_amp = aafplots_amp)


# Render the R Markdown file with parameters
render(summary_report_rmd,
       output_format = 'html_document',
       params = params)