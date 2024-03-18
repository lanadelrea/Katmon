library(rmarkdown)

---
title: "SARS-CoV-2 Variant Coinfection"
author: "A.M. Realingo"
date: "`r format(Sys.Date(), '%B %d, %Y')`"
output: pdf_document
---
  

## Section 1: Bammix
figure1 <- args[1]
![Figure 1 Bammix plot](figure1)

## Section 2: Freyja
figure2 <- args[2]
![Figure 2 Freyja](figure2)

## Section 3: Alternative Allele Fraction per Mutation
figure3 <- args[3]
![Figure 3 Alternative Allele Fraction per Mutation](figure3)

## Section 4: Alternative Allele Fraction pero Amplicon
figure4 <- args[4]
![Figure 4 Alternative Allele Fraction per Amplicon](figure4)