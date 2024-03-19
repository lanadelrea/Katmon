args <- commandArgs(trailingOnly = TRUE)

rmarkdown::render("summary-report.Rmd",
                  params = list(bammix_plot = args[1],
                                freyja_plot = args[2],
                                aafplot_mut = args[3],
                                aafplot_amp = args[4]))
