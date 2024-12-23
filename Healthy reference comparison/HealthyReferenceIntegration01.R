rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(harmony)
library(rafalib)
library(tidyverse)
library(ggplot2)
library(clustree)


# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_FullMetadata.RDS"
)
liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2a/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)


# ==== Toronto Dropbox ====
toronto_integrated = readRDS(
  "/data/work/Projects/HCA/AndrewsEtAlLiver/Integrated_with_Subannotations.rds"
)

# filter data: single nuclei
toronto_integrated = toronto_integrated[,toronto_integrated$assay_type == "single_nuc"]
# QC-filters
liver_metadata[,c("percent_mito","nCount_RNA","nFeature_RNA")] %>% summary
toronto_integrated = toronto_integrated[,toronto_integrated$percent.mt <= 20 & toronto_integrated$nCount_RNA >= 100 & toronto_integrated$nFeature_RNA >= 100]

# ==== Joint Object =====
joint_metadata = data.frame(
  Batch = c(liver_metadata$Batch,toronto_integrated$sample),
  Condition = c(
    rep("COVID",nrow(liver_metadata)),
    rep("Control",nrow(toronto_integrated@meta.data))
  )
)

rownames(joint_metadata) = c(
  rownames(liver_metadata),rownames(toronto_integrated@meta.data)
)

common_genes = intersect(
  rownames(toronto_integrated@assays$RNA@counts),
  rownames(liver_counts)
)


joint_so = CreateSeuratObject(
  counts = cbind(liver_counts[common_genes,rownames(liver_metadata)],toronto_integrated@assays$RNA@counts[common_genes,]),
  meta.data = joint_metadata
)

# ===== Processing =====
joint_so = NormalizeData(joint_so)
# standard processing
joint_so = ScaleData(joint_so,features=rownames(joint_so))
joint_so = FindVariableFeatures(joint_so)
joint_so = RunPCA(
  object=joint_so, assay="RNA", 
  features = VariableFeatures(joint_so)
)


# ==== Joint Embedding ====
# Determine percent of variation associated with each PC
pct <- joint_so[["pca"]]@stdev / sum(joint_so[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
# Minimum of the two calculation
pcs <- min(co1, co2)

# Harmony
set.seed(664)
joint_so = RunHarmony(object=joint_so,group.by.vars = "Batch",max.iter.harmony = 30,plot_convergence = T)

pcs=30
joint_so = FindNeighbors(object = joint_so, reduction = "harmony", dims = 1:pcs)
joint_so = RunUMAP(object = joint_so, reduction = "harmony", dims = 1:pcs, return.model = T)

DimPlot(joint_so,group.by = "Condition")

# ===== Clustering (Compartments) =====
# 0.005, 0.007, 0.009, 0.01
joint_so = FindClusters(joint_so,resolution = 0.01)

saveRDS(
  joint_so,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/COVIDControlJointSeuratObject.RDS"
)


