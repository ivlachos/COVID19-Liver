rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(genefilter)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(Matrix)
library(parallel)

# ==== Liver Data ====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1.RDS"
)

umapLst = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl3a/data/SelectLiverNoDoubletsUMAPembeddingsLst.RDS"
)

compartment_name = as.character(0:4)
names(compartment_name) = as.character(c(0,2,1,3,4))


# Joint Embedding
joint_so = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/COVIDControlJointSeuratObject.RDS"
  
)

umapComp = do.call(rbind,umapLst[as.character(0:4)])

# ==== Joint Embedding Results =====
TorontoPredNN = lapply(
  as.character(0:4),
  function(ck){
    readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
        "COVIDControlJointSeuratObject_Cluster",ck,"_ControlPred.RDS"
      )
    )
  }
)

TorontoPredNN = do.call(rbind,TorontoPredNN)
TorontoPredNN$Compartment = sapply(strsplit(TorontoPredNN$SubCluster,"_"),head,1) 

TorontoPredSNN = lapply(
  as.character(0:4),
  function(ck){
    readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
        "COVIDControlJointSeuratObject_Cluster",ck,"_snn_ControlPred.RDS"
      )
    )
  }
)

TorontoPredSNN = do.call(rbind,TorontoPredSNN)
TorontoPredSNN$Compartment = sapply(strsplit(TorontoPredSNN$SubCluster,"_"),head,1) 

TorontoPred = rbind(TorontoPredNN,TorontoPredSNN)


missing_controls = setdiff(
  rownames(joint_so@meta.data)[joint_so@meta.data$Condition == "Control"],
  rownames(TorontoPred)
)


# ==== Nearest Neighbors ====
# columns: COVID livers, rownames: Toronto
nn_mat = joint_so@graphs$RNA_nn[ missing_controls, rownames(liver_metadata) ]

control_nn = which(nn_mat == 1, arr.ind = TRUE)
control_nn = control_nn[order(control_nn[,1]),]

# get subclusters and compartment umap coordinates
control_nnCluster = data.frame(
  Control = rownames(nn_mat)[ control_nn[,1] ],
  SubCluster = liver_metadata[colnames(nn_mat)[ control_nn[,2] ],"SubCluster"],
  umapLst[["All"]][colnames(nn_mat)[ control_nn[,2] ],]
)
rownames(control_nnCluster) = NULL

# assign subcluster labels
tmpLst = split(x=control_nnCluster$SubCluster,f=control_nnCluster$Control)
control_predLabels = mclapply(
  tmpLst,
  function(x){
    tab = table(x)
    names(which.max(tab))
  },
  mc.cores=8
)
rm(tmpLst)

# assign global UMAP coordinates
tmpLst = split(x=control_nnCluster,f=control_nnCluster$Control)
control_gUMAP = mclapply(
  seq_along(tmpLst),
  function(ic){
    sck = control_predLabels[[ic]]
    tmp00 = tmpLst[[ic]]
    tmp00 = tmp00[tmp00$SubCluster == sck,]
    colMeans(tmp00[,c("UMAP_1","UMAP_2")])
  },
  mc.cores = 8
)
rm(tmpLst)

control_gUMAP = do.call(rbind,control_gUMAP)
colnames(control_gUMAP) = paste0("g",colnames(control_gUMAP))

# assign compartment UMAP coordinates
control_nnCluster2 = data.frame(
  Control = rownames(nn_mat)[ control_nn[,1] ],
  SubCluster = liver_metadata[colnames(nn_mat)[ control_nn[,2] ],"SubCluster"],
  umapComp[colnames(nn_mat)[ control_nn[,2] ],]
)
rownames(control_nnCluster2) = NULL

tmpLst = split(x=control_nnCluster2,f=control_nnCluster2$Control)
control_cUMAP = mclapply(
  seq_along(tmpLst),
  function(ic){
    sck = control_predLabels[[ic]]
    tmp00 = tmpLst[[ic]]
    tmp00 = tmp00[tmp00$SubCluster == sck,]
    colMeans(tmp00[,c("UMAP_1","UMAP_2")])
  },
  mc.cores = 8
)
rm(tmpLst)

control_cUMAP = do.call(rbind,control_cUMAP)
colnames(control_cUMAP) = paste0("c",colnames(control_cUMAP))



control_pred = data.frame(
  SubCluster = unlist(control_predLabels,use.names = T),
  control_gUMAP,
  control_cUMAP
)



# ==== Shared NN ====
# missing nn
tmp00 = rowSums(nn_mat == 1)
nn_miss = names(which(tmp00 == 0))
snn_mat = joint_so@graphs$RNA_snn[ nn_miss , rownames(liver_metadata) ]
rm(tmp00)

tmp00 = rowSums(snn_mat > 0)
snn_mat = snn_mat[tmp00 > 0,]


# S. Nearest Neighbors
snn_match = apply(snn_mat,1,which.max)

# get subclusters and compartment umap coordinates
control_snnCluster = data.frame(
  Control = names(snn_match),
  SubCluster = liver_metadata[rownames(liver_metadata)[snn_match],"SubCluster"],  
  umapLst[["All"]][rownames(liver_metadata)[snn_match],]
)
rownames(control_snnCluster) = NULL



# assign subcluster labels
tmpLst = split(x=control_snnCluster$SubCluster,f=control_snnCluster$Control)
control_snn_predLabels = mclapply(
  tmpLst,
  function(x){
    tab = table(x)
    names(which.max(tab))
  },
  mc.cores=8
)
rm(tmpLst)


# assign global UMAP coordinates
tmpLst = split(x=control_snnCluster,f=control_snnCluster$Control)
control_snn_gUMAP = mclapply(
  seq_along(tmpLst),
  function(ic){
    sck = control_snn_predLabels[[ic]]
    tmp00 = tmpLst[[ic]]
    tmp00 = tmp00[tmp00$SubCluster == sck,]
    colMeans(tmp00[,c("UMAP_1","UMAP_2")])
  },
  mc.cores = 8
)
rm(tmpLst)

control_snn_gUMAP = do.call(rbind,control_snn_gUMAP)
colnames(control_snn_gUMAP) = paste0("g",colnames(control_snn_gUMAP))

# assign compartment UMAP coordinates
control_snnCluster2 = data.frame(
  Control = names(snn_match),
  SubCluster = liver_metadata[rownames(liver_metadata)[snn_match],"SubCluster"], 
  umapComp[rownames(liver_metadata)[snn_match],]
)
rownames(control_snnCluster2) = NULL
tmpLst = split(x=control_snnCluster2,f=control_snnCluster2$Control)

control_snn_cUMAP = mclapply(
  seq_along(tmpLst),
  function(ic){
    sck = control_snn_predLabels[[ic]]
    tmp00 = tmpLst[[ic]]
    tmp00 = tmp00[tmp00$SubCluster == sck,]
    colMeans(tmp00[,c("UMAP_1","UMAP_2")])
  },
  mc.cores = 8
)
rm(tmpLst)

control_snn_cUMAP = do.call(rbind,control_snn_cUMAP)
colnames(control_snn_cUMAP) = paste0("c",colnames(control_snn_cUMAP))

control_snn_pred = data.frame(
  SubCluster = unlist(control_snn_predLabels,use.names = T),
  control_snn_gUMAP,
  control_snn_cUMAP
)

# ==== Save Results ====
saveRDS(
  control_pred,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "COVIDControlJointSeuratObject_MissingNN_ControlPred.RDS"
  )
)

saveRDS(
  control_snn_pred,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "COVIDControlJointSeuratObject_MissingSNN_ControlPred.RDS"
  )
)






