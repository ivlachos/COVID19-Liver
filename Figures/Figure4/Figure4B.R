
rm(list=ls())
options(stringsAsFactors = FALSE)

library(rafalib)
library(tidyverse)
library(pheatmap)
library(parallel)
library(ggplot2)
library(cowplot)


# ==== DSP COVID Probe Enrichment (Donor) =====
res_covid= readRDS(
  "LiverDSPCovidEnrichmentScores.RDS"
)


ggplot(
  res_covid,
  aes(x=Patient,y=COVIDScore)
) + 
  geom_boxplot() +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=11),
    legend.position = "none"
  ) +
  xlab("") + ylab("Enrichment Score")



