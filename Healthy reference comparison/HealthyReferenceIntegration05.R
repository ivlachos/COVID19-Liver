rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyverse)


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
TorontoPredNN$Projection = "NN"

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
TorontoPredSNN$Projection = "SNN"

TorontoPredMissNN = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "COVIDControlJointSeuratObject_MissingNN_ControlPred.RDS"
  )
)

TorontoPredMissNN$Compartment = sapply(strsplit(TorontoPredMissNN$SubCluster,"_"),head,1) 
TorontoPredMissNN$Projection = "MissNN"


TorontoPredMissSNN = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "COVIDControlJointSeuratObject_MissingSNN_ControlPred.RDS"
  )
)
TorontoPredMissSNN$Compartment = sapply(strsplit(TorontoPredMissSNN$SubCluster,"_"),head,1) 
TorontoPredMissSNN$Projection = "MissSNN"



TorontoPred = rbind(TorontoPredNN,TorontoPredMissNN,TorontoPredSNN,TorontoPredMissSNN)

rm(TorontoPredNN,TorontoPredMissNN,TorontoPredSNN,TorontoPredMissSNN)

TorontoPred$sample = sapply(
  strsplit(
    x=rownames(TorontoPred),split="_"
  ),
  function(x){
    paste(head(x,2),collapse="_")
  }
)

# save results
saveRDS(
  TorontoPred,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "COVIDControlJointSeuratObject_ControlPredAll.RDS"
  )
)
