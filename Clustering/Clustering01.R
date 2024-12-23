rm(list=ls())
options(stringsAsFactors = FALSE)


library(tidyverse)
library(rafalib)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(Seurat)
library(harmony)
library(cluster)
library(clValid)
library(clustree)
library(scales)
library(parallel)



# ==== Liver Data =====
liver_all = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObj.RDS"
)



# ==== Doublet Prediction (Scrublet) =====
liver_scrub = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAllScrubletScores.RDS"
)

setkey(liver_scrub,barcodekey)

# add Scrublet Score to liver meta data
liver_all$ScrubletScore = liver_scrub[colnames(liver_all)]$doublet_scores
liver_all$ScrubletPredicted = liver_scrub[colnames(liver_all)]$predicted_doublets



# ==== Doublet Prediction (Pegasus) =====
liver_pegasus = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAllPegasusScores.RDS"
)

setkey(liver_pegasus,barcodekey)



# add Pegasus Score to liver meta data
liver_all$PegasusScore = liver_pegasus[colnames(liver_all)]$doublet_score
liver_all$PegasusPredicted = liver_pegasus[colnames(liver_all)]$pred_dbl

liver_all$BothPredicted = liver_all$PegasusPredicted & liver_all$ScrubletPredicted
liver_all$EitherPredicted = liver_all$PegasusPredicted | liver_all$ScrubletPredicted


# ===== Doublet Association Test =====
pegasusLst = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversPegasusAssocDoublets.RDS"
)

scrubletLst = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletAssocDoublets.RDS"
)

eitherLst = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletORPegasusAssocDoublets.RDS"
)

bothLst= readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletANDPegasusAssocDoublets.RDS"
)


liver_all$ScrubletDoublets = rownames(liver_all@meta.data) %in% c(
  unlist(scrubletLst,use.names = FALSE),
  rownames(liver_all@meta.data)[liver_all$ScrubletPredicted]
) 

liver_all$PegasusDoublets = rownames(liver_all@meta.data) %in% c(
  unlist(pegasusLst,use.names = FALSE),
  rownames(liver_all@meta.data)[liver_all$PegasusPredicted]
) 


liver_all$EitherDoublets = rownames(liver_all@meta.data) %in% c(
  unlist(eitherLst,use.names = FALSE),
  rownames(liver_all@meta.data)[liver_all$EitherPredicted]
) 

liver_all$BothDoublets = rownames(liver_all@meta.data) %in% c(
  unlist(bothLst,use.names = FALSE),
  rownames(liver_all@meta.data)[liver_all$BothPredicted]
) 

# ==== Filter Doublets =====
clust_cols = hue_pal()(5)
names(clust_cols) = as.character(0:4)


liver_tmp = liver_all[,!(liver_all$PegasusDoublets)]

# Normalization
liver_tmp <- NormalizeData(object = liver_tmp)
liver_tmp <- FindVariableFeatures(liver_tmp, selection.method = "vst", nfeatures = 2000)

# PCA
# scale all genes
all.genes <- rownames(liver_tmp)
liver_tmp <- ScaleData(liver_tmp, features = all.genes)
liver_tmp <- RunPCA(liver_tmp, features = VariableFeatures(object = liver_tmp))
rm(all.genes)

# ===== Initial Clustering ======
liver_tmp = FindNeighbors(object=liver_tmp, reduction = "harmony",dims=1:15)
liver_tmp = RunUMAP(object=liver_tmp, reduction = "harmony",dims=1:15)
liver_tmp = FindClusters(liver_tmp,resolution=0.008,random.seed=665,group.singletons = TRUE)


DimPlot(liver_tmp)


saveRDS(
  liver_tmp,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)

