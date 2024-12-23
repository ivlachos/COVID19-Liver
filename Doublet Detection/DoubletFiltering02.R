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
# previous results
liver_prev_metadata = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiver/data/SelectLivermetadata.RDS"
)

# Hepatocyte-Like SubClusters
heplike_subclusters = c("1_3","1_9","2_4","2_6","3_5","3_7")

liver_all = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObj.RDS"
)
# liver_all@assays$PearsonResiduals@data[is.na(liver_all@assays$PearsonResiduals@data)] = 0

liver_all[["PearsonResiduals"]] = NULL
liver_all[["RISC"]] = NULL


# quick checks
DimPlot(
  liver_all
)

DimPlot(
  liver_all[,liver_all$SubCluster %in% heplike_subclusters],
  group.by = "SubCluster"
) + facet_wrap(vars(SubCluster)) + NoLegend()



# ==== Doublet Prediction (Scrublet) =====
liver_scrub = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAllScrubletScores.RDS"
)

setkey(liver_scrub,barcodekey)

# add Scrublet Score to liver meta data
liver_all$ScrubletScore = liver_scrub[colnames(liver_all)]$doublet_scores
liver_all$ScrubletPredicted = liver_scrub[colnames(liver_all)]$predicted_doublets


setkey(liver_scrub,barcodekey)


liver_prev_metadata$ScrubletScore = liver_scrub[rownames(liver_prev_metadata)]$doublet_scores
liver_prev_metadata$ScrubletPredicted = liver_scrub[rownames(liver_prev_metadata)]$predicted_doublets


rm(liver_scrub)

umapLst = readRDS("/data/shiny/shiny-server/apps/SelectLiver/data/SelectLiverUMAPembeddingsLst.RDS")

# ==== Doublet Prediction (Pegasus) =====
liver_pegasus = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAllPegasusScores.RDS"
)

setkey(liver_pegasus,barcodekey)



# add Pegasus Score to liver meta data
liver_all$PegasusScore = liver_pegasus[colnames(liver_all)]$doublet_score
liver_all$PegasusPredicted = liver_pegasus[colnames(liver_all)]$pred_dbl


setkey(liver_pegasus,barcodekey)

liver_prev_metadata$PegasusScore = liver_pegasus[rownames(liver_prev_metadata)]$doublet_scores
liver_prev_metadata$PegasusPredicted = liver_pegasus[rownames(liver_prev_metadata)]$predicted_doublets


rm(liver_pegasus)

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


plot_grid(
  DimPlot(liver_all,group.by = "PegasusDoublets"),
  DimPlot(liver_all,group.by = "ScrubletDoublets"),
  DimPlot(liver_all,group.by = "EitherDoublets"),
  DimPlot(liver_all,group.by = "BothDoublets")
)

liver_prev_metadata$PegasusDoublets = liver_all@meta.data[rownames(liver_prev_metadata),"PegasusDoublets"]
liver_prev_metadata$ScrubletDoublets = liver_all@meta.data[rownames(liver_prev_metadata),"ScrubletDoublets"]
liver_prev_metadata$EitherDoublets = liver_all@meta.data[rownames(liver_prev_metadata),"EitherDoublets"]
liver_prev_metadata$BothDoublets = liver_all@meta.data[rownames(liver_prev_metadata),"BothDoublets"]

# ==== Doublets Doublets Barplots II =====

# Barplot: Scrublet
doubletScrubletTab = 100*sapply(
  split(
    liver_all$ScrubletDoublets,
    liver_all$SubCluster
  ),
  mean
)


PredictScrubletTabSubClustAll = data.frame(
  SubCluster = names(doubletScrubletTab),
  Prop = as.numeric(c(doubletScrubletTab))
)
PredictScrubletTabSubClustAll$Compartment = sapply(strsplit(PredictScrubletTabSubClustAll$SubCluster,"_"),head,1)
PredictScrubletTabSubClustAll$Hepatocyte.Like = PredictScrubletTabSubClustAll$SubCluster %in% heplike_subclusters

bpScrubletSubclust = ggplot(
  PredictScrubletTabSubClustAll,
  aes(x=SubCluster,y=Prop,fill=Compartment)
) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust=0.95
    )
  ) + 
  xlab("") + ylab("%")

bpScrubletSubclust
bpScrubletSubclust + facet_wrap(vars(Hepatocyte.Like))



# Barplot: Pegasus
doubletPegasusTab = 100*sapply(
  split(
    liver_all$PegasusDoublets,
    liver_all$SubCluster
  ),
  mean
)


PredictPegasusTabSubClustAll = data.frame(
  SubCluster = names(doubletPegasusTab),
  Prop = as.numeric(c(doubletPegasusTab))
)
PredictPegasusTabSubClustAll$Compartment = sapply(strsplit(PredictPegasusTabSubClustAll$SubCluster,"_"),head,1)
PredictPegasusTabSubClustAll$Hepatocyte.Like = PredictPegasusTabSubClustAll$SubCluster %in% heplike_subclusters

bpPegasusSubclust = ggplot(
  PredictPegasusTabSubClustAll,
  aes(x=SubCluster,y=Prop,fill=Compartment)
) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust=0.95
    )
  ) + 
  xlab("") + ylab("%")

bpPegasusSubclust
bpPegasusSubclust + facet_wrap(vars(Hepatocyte.Like))



# Barplot: Scrublet OR Pegasus
doubletORTab = 100*sapply(
  split(
    liver_all$EitherDoublets,
    liver_all$SubCluster
  ),
  mean
)


PredictORTabSubClustAll = data.frame(
  SubCluster = names(doubletORTab),
  Prop = as.numeric(c(doubletORTab))
)
PredictORTabSubClustAll$Compartment = sapply(strsplit(PredictORTabSubClustAll$SubCluster,"_"),head,1)
PredictORTabSubClustAll$Hepatocyte.Like = PredictORTabSubClustAll$SubCluster %in% heplike_subclusters

bpORSubclust = ggplot(
  PredictORTabSubClustAll,
  aes(x=SubCluster,y=Prop,fill=Compartment)
) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust=0.95
    )
  ) + 
  xlab("") + ylab("%")

bpORSubclust
bpORSubclust + facet_wrap(vars(Hepatocyte.Like))


# Barplot: Scrublet AND Pegasus
doubletANDTab = 100*sapply(
  split(
    liver_all$BothDoublets,
    liver_all$SubCluster
  ),
  mean
)


PredictANDTabSubClustAll = data.frame(
  SubCluster = names(doubletANDTab),
  Prop = as.numeric(c(doubletANDTab))
)
PredictANDTabSubClustAll$Compartment = sapply(strsplit(PredictANDTabSubClustAll$SubCluster,"_"),head,1)
PredictANDTabSubClustAll$Hepatocyte.Like = PredictANDTabSubClustAll$SubCluster %in% heplike_subclusters

bpANDSubclust = ggplot(
  PredictANDTabSubClustAll,
  aes(x=SubCluster,y=Prop,fill=Compartment)
) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust=0.95
    )
  ) + 
  xlab("") + ylab("%")

bpANDSubclust
bpANDSubclust + facet_wrap(vars(Hepatocyte.Like))



PredictScrubletTabSubClustAll$Class = "Scrublet"
PredictPegasusTabSubClustAll$Class = "Pegasus"
PredictANDTabSubClustAll$Class = "Scrublet AND Pegasus"
PredictORTabSubClustAll$Class = "Scrublet OR Pegasus"

DoubletsDoubletTab = rbind(
  PredictScrubletTabSubClustAll,PredictPegasusTabSubClustAll,
  PredictORTabSubClustAll,PredictANDTabSubClustAll
)

ggplot(
  DoubletsDoubletTab,
  aes(x=SubCluster,y=Prop,fill=Compartment)
) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust=0.95
    )
  ) + 
  xlab("") + ylab("%") +
  facet_wrap(vars(Class))



# ===== UMAP Plots (Doublets) ====
# Pegasus
umapDF = data.frame(
  umapLst[["All"]],
  liver_prev_metadata[rownames(umapLst[["All"]]),]
)

pumapPegasusAll = ggplot(
  umapDF,
  aes(x=UMAP_1,y=UMAP_2)
) + 
  geom_point(size=0.15, color = "lightgray") +
  geom_point(
    data = umapDF %>% dplyr::filter(PegasusDoublets),
    aes(x=UMAP_1,y=UMAP_2),
    color="red",size=0.25,
    inherit.aes = FALSE
  ) +
  theme_cowplot() +
  ggtitle("Pegasus Doublets")

rm(umapDF)



umapPegasusLst = lapply(
  as.character(0:4),
  function(ck){
    umapDF = data.frame(
      umapLst[[ck]],
      liver_prev_metadata[rownames(umapLst[[ck]]),]
    )
    
    ggplot(
      umapDF,
      aes(x=UMAP_1,y=UMAP_2)
    ) + 
      geom_point(size=0.15, color = "lightgray") +
      geom_point(
        data = umapDF %>% dplyr::filter(PegasusDoublets),
        aes(x=UMAP_1,y=UMAP_2),
        color="red",size=0.25,
        inherit.aes = FALSE
      ) +
      theme_cowplot() +
      ggtitle(paste0("Compartment ",ck))
    
    
  }
)

# Scrublet
umapDF = data.frame(
  umapLst[["All"]],
  liver_prev_metadata[rownames(umapLst[["All"]]),]
)

pumapScrubletAll = ggplot(
  umapDF,
  aes(x=UMAP_1,y=UMAP_2)
) + 
  geom_point(size=0.15, color = "lightgray") +
  geom_point(
    data = umapDF %>% dplyr::filter(ScrubletDoublets),
    aes(x=UMAP_1,y=UMAP_2),
    color="red",size=0.25,
    inherit.aes = FALSE
  ) +
  theme_cowplot() +
  ggtitle("Scrublet Doublets")

rm(umapDF)



umapScrubletLst = lapply(
  as.character(0:4),
  function(ck){
    umapDF = data.frame(
      umapLst[[ck]],
      liver_prev_metadata[rownames(umapLst[[ck]]),]
    )
    
    ggplot(
      umapDF,
      aes(x=UMAP_1,y=UMAP_2)
    ) + 
      geom_point(size=0.15, color = "lightgray") +
      geom_point(
        data = umapDF %>% dplyr::filter(ScrubletDoublets),
        aes(x=UMAP_1,y=UMAP_2),
        color="red",size=0.25,
        inherit.aes = FALSE
      ) +
      theme_cowplot() +
      ggtitle(paste0("Compartment ",ck))
    
    
  }
)


plot_grid(plotlist = umapScrubletLst)

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

liver_tmp = FindNeighbors(object=liver_tmp, reduction = "harmony",dims=1:15)
liver_tmp = RunUMAP(object=liver_tmp, reduction = "harmony",dims=1:15)
liver_tmp = FindClusters(liver_tmp,resolution=0.008,random.seed=665,group.singletons = TRUE)


DimPlot(liver_tmp)


saveRDS(
  liver_tmp,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)



# ===== SubClusters ======
SubClustResolution = c(0.33,0.3,0.5,0.4,0.3)
names(SubClustResolution) = as.character(0:4)
clusterLst = list()
umapLstAlt = list()
for(ic in seq_along(SubClustResolution)){
  CurrentCluster = names(SubClustResolution)[ic]
  # Subset data
  tmpso = liver_all[,liver_all$seurat_clusters == CurrentCluster]
  # Normalization
  tmpso <- NormalizeData(object = tmpso)
  tmpso <- FindVariableFeatures(tmpso, selection.method = "vst", nfeatures = 2000)
  # PCA
  # scale all genes
  all.genes <- rownames(tmpso)
  tmpso <- ScaleData(tmpso, features = all.genes)
  # PCA
  tmpso <- RunPCA(tmpso, features = VariableFeatures(object = tmpso))
  # harmony
  tmpso = RunHarmony(tmpso,"Batch", plot_convergence = TRUE,assay.use="RNA")
  tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:15)
  tmpso = RunUMAP(object=tmpso, reduction = "harmony",dims=1:15)
  tmpso = FindClusters(tmpso,resolution=SubClustResolution[ic])
  clusterLst[[CurrentCluster]] = tmpso$seurat_clusters
  umapLstAlt[[CurrentCluster]] = tmpso@reductions$umap@cell.embeddings
  rm(tmpso)
}

saveRDS(
  clusterLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLivers_NoDoublets_ClusterLst.RDS"
)
saveRDS(
  umapLstAlt,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLivers_NoDoublets_umapEmbeddingsLst.RDS"
)

