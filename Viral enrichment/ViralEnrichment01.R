rm(list=ls())
options(stringsAsFactors = FALSE)


library(parallel)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(rafalib)
library(matrixStats)
library(googledrive)
library(readxl)


# ==== Cluster Manual Annotation ====
subcluster_annot = read_xlsx(
  path="/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/subcluster_annotation.xlsx",
  sheet=1
) %>% as.data.frame()

# cluster names used in manuscript
short_names = subcluster_annot$Label
names(short_names) = subcluster_annot$`Clustering v.1.1`

compartment_names = c(
  "0" = "Hepatocytes",
  "1" = "Immune",
  "2" = "Endothelial",
  "3" = "Mesenchymal",
  "4" = "BECs"
)



# ==== Functions ====
EnrichScore = function(clust_vcells,total_vcells,clust_prop,epsilon = 0.0001){
  log(
    (clust_vcells + epsilon)/(total_vcells*clust_prop + epsilon)
  )
}

# ==== COVID Liver =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1_DonorIDs.RDS"
)

tmp00 = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_BootCovidTest.RDS"
)
liver_metadata$Boot.SARS.CoV.2 = tmp00[rownames(liver_metadata),"Boot.SARS.CoV.2"]
rm(tmp00)

liver_metadata$Batch2 = gsub("ncLab","",liver_metadata$Batch)
liver_metadata$Batch2 = gsub("Broad","",liver_metadata$Batch2)


covid_umi_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsCOVIDcounts.RDS"
)

liver_metadata$SARS.CoV.2.UMI = covid_umi_counts[,rownames(liver_metadata)]
rm(covid_umi_counts)


