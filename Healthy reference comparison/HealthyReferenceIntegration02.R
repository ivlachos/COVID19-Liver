rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(harmony)
library(rafalib)
library(tidyverse)
library(ggplot2)

# ==== Functions =====
# Helper function to pick number of components
PickComponentNum = function(tmpso){
  # Determine percent of variation associated with each PC
  pct <- tmpso[["pca"]]@stdev / sum(tmpso[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  return(pcs)
}

# ===== Joint Embedding ====
joint_so = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/COVIDControlJointSeuratObject.RDS"
)




for(ck in levels(joint_so$seurat_clusters)){
  cat(
    paste0(
      "\n---------------------------\n",
      "Joint Cluster: \t\t",ck,"\n",
      "\n---------------------------\n"
    )
  )
  
  
  # ==== Compartment Clustering ====
  cat("Normalization...\n")
  tmpso = joint_so[,joint_so$seurat_clusters == ck]
  tmpso = NormalizeData(tmpso)
  # standard processing
  tmpso = ScaleData(tmpso,features=rownames(tmpso))
  tmpso = FindVariableFeatures(tmpso)
  
  cat("PCA/Harmony...\n")
  tmpso = RunPCA(
    object=tmpso, assay="RNA", 
    features = VariableFeatures(tmpso)
  )
  # Harmony
  set.seed(664)
  tmpso = RunHarmony(object=tmpso,group.by.vars = "Batch",max.iter.harmony = 30,plot_convergence = T)
  # pick number of PCs
  pcs = PickComponentNum(tmpso)
  
  cat("Joint embedding...\n")
  tmpso = FindNeighbors(object = tmpso, reduction = "harmony", dims = 1:pcs)
  tmpso = RunUMAP(object = tmpso, reduction = "harmony", dims = 1:pcs, return.model = T)
  tmpso = FindClusters(tmpso, resolution = 0.5)
  
  
  DimPlot(
    tmpso,
    label = TRUE
  )
  
  cat("Saving object...\n")
  saveRDS(
    tmpso,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
      "COVIDControlJointSeuratObject_Cluster",ck,".RDS"
    )
  )
  
  rm(tmpso)
  cat(paste0("\n\nDONE!\n\n"))
}



