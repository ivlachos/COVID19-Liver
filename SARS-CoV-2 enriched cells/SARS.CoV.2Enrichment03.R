rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(Matrix)
library(parallel)
library(ggplot2)
library(cowplot)
library(tidyverse)

# ==== Liver Data ====
liver_metadata = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsMetadata6.RDS"
)

covid_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsCOVIDcounts.RDS"
)

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)

# ===== Abundance Tests =====
bootDF = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversCOVIDBootTestFull.RDS"
)
rownames(bootDF) = bootDF$Nuclei

# subset
bootDF = bootDF[rownames(liver_metadata),]


# ==== Assign Status ====
covid_zero = colnames(covid_counts)[c(covid_counts) == 0]
covid_nonzero = setdiff(colnames(covid_counts),covid_zero)

boot_ambient = covid_nonzero[bootDF[covid_nonzero,"FDR"] >= 0.01]
boot_plus = covid_nonzero[bootDF[covid_nonzero,"FDR"] < 0.01]
length(covid_zero) + length(boot_ambient) + length(boot_plus)

liver_metadata$Boot.SARS.CoV.2 = "Ambient.RNA"
liver_metadata[covid_zero,"Boot.SARS.CoV.2"] = "RNA.Minus"
liver_metadata[boot_ambient,"Boot.SARS.CoV.2"] = "RNA.Ambient"
liver_metadata[boot_plus,"Boot.SARS.CoV.2"] = "RNA.Plus"
liver_metadata$Boot.SARS.CoV.2 = factor(liver_metadata$Boot.SARS.CoV.2 ,levels=c("RNA.Plus","RNA.Ambient","RNA.Minus"))

# ==== Plots ====
pbboot = liver_metadata[,c("SubCluster","Boot.SARS.CoV.2")] %>%
  table %>%
  prop.table(margin=1) %>%
  reshape2::melt() %>%
  dplyr::filter(Boot.SARS.CoV.2 == "RNA.Plus") %>%
  ggplot(
    aes(x=SubCluster,y=100*value,fill=Boot.SARS.CoV.2)
  )  +
  geom_bar(stat="identity") +
  ylab("%") + xlab("") +
  theme(
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
    legend.position = "none"
  )

pbboot

# ===== Save Results ====
saveRDS(
  object = liver_metadata,
  file = "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_BootCovidTest.RDS"
)

