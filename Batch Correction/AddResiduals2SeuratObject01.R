rm(list=ls())
options(stringsAsFactors = FALSE)


library(Seurat)

# ==== Seurat Object =====
liver_all = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)

# ==== Residuals =====
DevianceResiduals = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoubletsGLMGamPoiDevianceResiduals.RDS"
)

# add assay
liver_all[["DevianceResiduals"]] = CreateAssayObject(data=DevianceResiduals)

# ===== SubClusters ======
clusterLst = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLivers_NoDoublets_ClusterLst.RDS"
)

SubClustLabels = lapply(
  seq_along(clusterLst),
  function(ic){
    x = clusterLst[[ic]]
    tmpNames = names(x)
    res = paste0(names(clusterLst)[ic],"_",as.character(x))
    names(res) = tmpNames
    return(res)
  }
)

SubClustLabels = unlist(SubClustLabels,use.names = TRUE)


liver_all$SubCluster = SubClustLabels[rownames(liver_all@meta.data)] 


saveRDS(
  liver_all,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)