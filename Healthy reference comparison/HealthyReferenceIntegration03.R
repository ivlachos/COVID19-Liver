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




for(ck in as.character(0:4)){
  cat(
    paste0(
      "\n---------------------------\n",
      "Joint Cluster: \t\t",ck,"\n",
      "COVID Compartment: \t",compartment_name[ck],
      "\n---------------------------\n"
    )
  )
  
  
  # sub set COVID Liver meta data and compartment UMAP
  subset_metadata = liver_metadata[liver_metadata$seurat_clusters == compartment_name[ck] ,]
  subset_umap = umapLst[[ compartment_name[ck] ]]
  
  # read subset of joint embedding
  tmpso = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
      "COVIDControlJointSeuratObject_Cluster",ck,".RDS"
    )
  )
  
  subset_metadata = subset_metadata[rownames(subset_metadata) %in%  rownames(tmpso@meta.data),]
  subset_umap = subset_umap[rownames(subset_metadata),]
  
  
  # columns: COVID livers, rownames: Toronto
  nn_mat = tmpso@graphs$RNA_nn[ rownames(tmpso@meta.data)[tmpso@meta.data$Condition == "Control"] , rownames(subset_metadata) ]
  
  
  # missing nn
  tmp00 = rowSums(nn_mat == 1)
  nn_miss = names(which(tmp00 == 0))
  cat(
    paste0(
      "Missing ",length(nn_miss),"\n\n"
    )
  )
  snn_mat = tmpso@graphs$RNA_snn[ nn_miss , rownames(subset_metadata) ]
  rm(tmp00)
  
  
  tmp00 = rowSums(snn_mat > 0)
  snn_mat = snn_mat[tmp00 > 0,]
  
  # S. Nearest Neighbors
  snn_match = apply(snn_mat,1,which.max)
  
  # get subclusters and compartment umap coordinates
  control_snnCluster = data.frame(
    Control = names(snn_match),
    SubCluster = subset_metadata[rownames(subset_metadata)[snn_match],"SubCluster"],  
    subset_umap[rownames(subset_metadata)[snn_match],]
  )
  rownames(control_snnCluster) = NULL
  
  
  
  # assign subcluster labels
  tmpLst = split(x=control_snnCluster$SubCluster,f=control_snnCluster$Control)
  control_predLabels = mclapply(
    tmpLst,
    function(x){
      tab = table(x)
      names(which.max(tab))
    },
    mc.cores=8
  )
  rm(tmpLst)
  
  # assign compartment UMAP coordinates
  tmpLst = split(x=control_snnCluster,f=control_snnCluster$Control)
  
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
  
  
  # assign global UMAP coordinates
  control_snnCluster2 = data.frame(
    Control = names(snn_match),
    SubCluster = subset_metadata[rownames(subset_metadata)[snn_match],"SubCluster"],  
    umapLst[["All"]][rownames(subset_metadata)[snn_match],]
  )
  rownames(control_snnCluster2) = NULL
  
  
  tmpLst = split(x=control_snnCluster2,f=control_snnCluster2$Control)
  
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
  
  
  control_cUMAP = do.call(rbind,control_cUMAP)
  colnames(control_cUMAP) = paste0("c",colnames(control_cUMAP))
  
  
  control_gUMAP = do.call(rbind,control_gUMAP)
  colnames(control_gUMAP) = paste0("g",colnames(control_gUMAP))
  
  
  control_pred = data.frame(
    SubCluster = unlist(control_predLabels,use.names = T),
    control_gUMAP,
    control_cUMAP
  )
  
  saveRDS(
    control_pred,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
      "COVIDControlJointSeuratObject_Cluster",ck,"_snn_ControlPred.RDS"
    )
  )
  
  cat("\n\nDONE!\n\n")
  
}
