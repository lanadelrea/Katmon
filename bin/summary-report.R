library(rmarkdown)
library(knitr)
library(kableExtra)
library(readr)


args <- commandArgs(trailingOnly = TRUE)


# Input table and plots

lineage_table_tsv <- args[1]
lineage_table <- read_tsv(lineage_table_tsv, locale(encoding = "UTF-8"))

bammix_plots <- readLines(args[2], encoding = "UTF-8") %>%
                strsplit("\t") %>%
                unlist()

freyja_plot <- args[3]

aafplots_mut <- readLines(args[4], encoding = "UTF-8") %>%
                strsplit("\t") %>%
                unlist()

aafplots_amp <- readLines(args[5], encoding = "UTF-8") %>%
                strsplit("\t") %>%
                unlist()

virstrain_table_tsv <- args[6]
virstrain_table <- read_tsv(virstrain_table_tsv, locale(encoding = "UTF-8"))

summary_report_rmd = args[7]

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