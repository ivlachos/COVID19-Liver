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

# ==== Functions =====
ClusterHelperFun = function(
  tissue_so, niter = 10,npc=10, 
  cell_fraction = 0.1, resolution_grid = seq(0.01,0.05,0.01),
  RecomputeSNN=TRUE,
  mt_exclude = FALSE
){
  
  if(cell_fraction < 1){
    # subsample cells
    set.seed(17)
    rndCells = sample(x=Cells(tissue_so),size=ceiling( length(Cells(tissue_so))*cell_fraction ))
    tmpso = tissue_so[,rndCells]
  }else{
    tmpso = tissue_so
  }
  
  
  if(RecomputeSNN){
    
    # # Recompute SNN graph using Harmony embedding
    # tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:npc)
    # 
    tmpso <- FindVariableFeatures(tmpso, selection.method = "vst", nfeatures = 2100)
    
    # PCA
    # scale all genes
    all.genes <- rownames(tmpso)
    tmpso <- ScaleData(tmpso, features = all.genes)
    
    # PCA
    if(mt_exclude){
      hvg =  setdiff(VariableFeatures(object = tmpso),mt_genes)
    }else{
      hvg = VariableFeatures(object = tmpso)
    }
    
    tmpso <- RunPCA(tmpso, features = hvg[1:2000])
    
    # harmony
    tmpso = RunHarmony(tmpso,"Batch", plot_convergence = TRUE,assay.use="RNA")
    tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:npc)
  }
  
  
  # Get distance based on Harmony 
  subset_dist = dist(x=Embeddings(object=tmpso[["harmony"]])[,1:npc])
  
  # Test different settings
  silhouetteLst = list()
  dunnLst = list()
  clustLst = list()
  for(rs in resolution_grid){
    silhoutte_coef = c()
    dunn_index = c()
    clust_n = c()
    for(ic in 1:niter){
      # clustering
      tmpso = FindClusters(tmpso,resolution=rs,random.seed=ic)
      clusters=tmpso$seurat_clusters
      if(length(levels(clusters)) > 1){
        # number of clusters
        clust_n = c(clust_n,length(unique(clusters)))
        
        # Silhoutte score
        sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = subset_dist)
        silhoutte_coef = c(silhoutte_coef,mean(sil[,"sil_width"]))
        
        # Dunn Index
        dunn_index = c(dunn_index,clValid::dunn(distance = subset_dist, clusters = as.numeric(x = as.factor(x = clusters))))
        
      }
      
    }
    silhouetteLst[[as.character(rs)]] = silhoutte_coef
    dunnLst[[as.character(rs)]] = dunn_index
    clustLst[[as.character(rs)]] = clust_n
  }
  
  
  
  # get median silhouette score
  median_silhoutte = sapply(silhouetteLst,median)
  
  # get median dunn index
  median_dunn = sapply(dunnLst,median)
  
  # rafalib::mypar(3,1)
  # boxplot(silhouetteLst)
  # boxplot(dunnLst)
  # boxplot(clustLst)
  
  # pick optimal parameter based on borda vote on Dunn index and Silhoutte score
  index_rank = rank(median_silhoutte,ties.method = "min") + rank(median_dunn,ties.method = "min")
  final_resolution = which.max(index_rank) %>% names %>% as.numeric
  
  
  resLst = list(
    silhouetteLst = silhouetteLst,
    dunnLst = dunnLst,
    clustSizeLst = clustLst,
    final_resolution = final_resolution
  )
  
  return(resLst)
}




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


# ===== Doublet Prediction (COVID Paper) =====
# final calls
bo_doublets = readLines(
  "/data/work/Projects/BrScRNAseq/data/covid_autopsy_liver_doublets.txt"
)

liver_all$BoDoublets = colnames(liver_all) %in% bo_doublets



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


# ==== Predicted Doublets Barplots =====

# Barplot: Scrublet
doubletScrubletTab = 100*sapply(
  split(
    liver_all$ScrubletPredicted,
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
    liver_all$PegasusPredicted,
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
    liver_all$PegasusPredicted | liver_all$ScrubletPredicted,
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
    liver_all$PegasusPredicted & liver_all$ScrubletPredicted,
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

PredictedDoubletTab = rbind(
  PredictScrubletTabSubClustAll,PredictPegasusTabSubClustAll,
  PredictORTabSubClustAll,PredictANDTabSubClustAll
)

ggplot(
  PredictedDoubletTab,
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


# ===== Doublet Association Test =====
liver_all$BothPredicted = liver_all$PegasusPredicted & liver_all$ScrubletPredicted
liver_all$EitherPredicted = liver_all$PegasusPredicted | liver_all$ScrubletPredicted

DimPlot(liver_all,group.by = "ScrubletPredicted")
DimPlot(liver_all,group.by = "PegasusPredicted")
DimPlot(liver_all,group.by = "BothPredicted")
DimPlot(liver_all,group.by = "EitherPredicted")

pegasusLst = list()
scrubletLst = list()
eitherLst = list()
bothLst = list()

for(ck in as.character(0:4)){
  # subset
  tmpso = liver_all[,liver_all$seurat_clusters == ck]
  
  # Normalization
  tmpso <- NormalizeData(object = tmpso)
  tmpso <- FindVariableFeatures(tmpso, selection.method = "vst", nfeatures = 2000)
  
  # PCA
  # scale all genes
  all.genes <- rownames(tmpso)
  tmpso <- ScaleData(tmpso, features = all.genes)
  tmpso <- RunPCA(tmpso, features = VariableFeatures(object = tmpso))
  rm(all.genes)
  
  # Harmony: after adjustment
  tmpso = RunHarmony(tmpso,"Batch", plot_convergence = FALSE,assay.use="RNA")
  tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:30)
  tmpso = RunUMAP(object=tmpso, reduction = "harmony",dims=1:30)
  
  print(plot_grid(
    DimPlot(tmpso,group.by = "ScrubletPredicted"),
    DimPlot(tmpso,group.by = "PegasusPredicted"),
    DimPlot(tmpso,group.by = "BothPredicted"),
    DimPlot(tmpso,group.by = "EitherPredicted")
  ))
  
  
  tmpso =  FindClusters(tmpso,resolution=3,group.singletons = TRUE)
  
  # Pegasus
  pegasus_doublet_pct = 100*sapply(
    split(tmpso$PegasusPredicted,tmpso$seurat_clusters),
    mean
  )
  
  pegasus_fisher_pvals = mclapply(
    levels(tmpso$seurat_clusters),
    function(sck){
      pegTab = table(data.frame(
        PegasusPredicted = tmpso$PegasusPredicted,
        CurrentCluster = tmpso$seurat_clusters == sck
      ))
      fisher.test(pegTab)$p.value
    },
    mc.cores=6
  )
  names(pegasus_fisher_pvals) = levels(tmpso$seurat_clusters)
  pegasus_fisher_pvals = unlist(pegasus_fisher_pvals)
  
  
  # Scrublet
  scrublet_doublet_pct = 100*sapply(
    split(tmpso$ScrubletPredicted,tmpso$seurat_clusters),
    mean
  )
  
  scrublet_fisher_pvals = mclapply(
    levels(tmpso$seurat_clusters),
    function(sck){
      pegTab = table(data.frame(
        ScrubletPredicted = tmpso$ScrubletPredicted,
        CurrentCluster = tmpso$seurat_clusters == sck
      ))
      fisher.test(pegTab)$p.value
    },
    mc.cores=6
  )
  names(scrublet_fisher_pvals) = levels(tmpso$seurat_clusters)
  scrublet_fisher_pvals = unlist(scrublet_fisher_pvals)
  
  
  # Either
  # Scrublet or Pegasus
  either_doublet_pct = 100*sapply(
    split(tmpso$EitherPredicted,tmpso$seurat_clusters),
    mean
  )
  
  either_fisher_pvals = mclapply(
    levels(tmpso$seurat_clusters),
    function(sck){
      pegTab = table(data.frame(
        EitherPredicted = tmpso$EitherPredicted,
        CurrentCluster = tmpso$seurat_clusters == sck
      ))
      fisher.test(pegTab)$p.value
    },
    mc.cores=6
  )
  names(either_fisher_pvals) = levels(tmpso$seurat_clusters)
  either_fisher_pvals = unlist(either_fisher_pvals)
  
  
  # Both
  # Scrublet and Pegasus
  both_doublet_pct = 100*sapply(
    split(tmpso$BothPredicted,tmpso$seurat_clusters),
    mean
  )
  
  both_fisher_pvals = mclapply(
    levels(tmpso$seurat_clusters),
    function(sck){
      pegTab = table(data.frame(
        BothPredicted = tmpso$BothPredicted,
        CurrentCluster = tmpso$seurat_clusters == sck
      ))
      fisher.test(pegTab)$p.value
    },
    mc.cores=6
  )
  names(both_fisher_pvals) = levels(tmpso$seurat_clusters)
  both_fisher_pvals = unlist(both_fisher_pvals)
  
  # Get results
  pegasus_test = data.frame(
    Percent = pegasus_doublet_pct,
    Pval = pegasus_fisher_pvals
  )
  pegasus_test$FDR = p.adjust(p=pegasus_test$Pval,method = "fdr")
  
  scrublet_test = data.frame(
    Percent = scrublet_doublet_pct,
    Pval = scrublet_fisher_pvals
  )
  scrublet_test$FDR = p.adjust(p=scrublet_test$Pval,method = "fdr")
  
  
  either_test = data.frame(
    Percent = either_doublet_pct,
    Pval = either_fisher_pvals
  )
  either_test$FDR = p.adjust(p=either_test$Pval,method = "fdr")
  
  
  
  both_test = data.frame(
    Percent = both_doublet_pct,
    Pval = both_fisher_pvals
  )
  both_test$FDR = p.adjust(p=both_test$Pval,method = "fdr")
  
  # assign doublets
  pegasus_picks = pegasus_test$Percent > 60 & pegasus_test$FDR < 0.05
  scrublet_picks = scrublet_test$Percent > 60 & scrublet_test$FDR < 0.05
  either_picks = either_test$Percent > 60 & either_test$FDR < 0.05
  both_picks = both_test$Percent > 60 & both_test$FDR < 0.05
  
  if(sum(pegasus_picks) > 0){
    pegasus_assoc_doublets = rownames(tmpso@meta.data)[tmpso@meta.data$seurat_clusters %in% rownames(pegasus_test)[pegasus_picks]]
  }else{
    pegasus_assoc_doublets = NA
  }
  
  if(sum(scrublet_picks) > 0){
    scrublet_assoc_doublets = rownames(tmpso@meta.data)[tmpso@meta.data$seurat_clusters %in% rownames(scrublet_test)[scrublet_picks]]
  }else{
    scrublet_assoc_doublets = NA
  }
  
  if(sum(either_picks) > 0){
    either_assoc_doublets = rownames(tmpso@meta.data)[tmpso@meta.data$seurat_clusters %in% rownames(either_test)[either_picks]]
  }else{
    either_assoc_doublets = NA
  }
  
  if(sum(both_picks) > 0){
    both_assoc_doublets = rownames(tmpso@meta.data)[tmpso@meta.data$seurat_clusters %in% rownames(both_test)[both_picks]]
  }else{
    both_assoc_doublets = NA
  }
  
  pegasusLst[[ck]] = pegasus_assoc_doublets
  scrubletLst[[ck]] = scrublet_assoc_doublets
  eitherLst[[ck]] = either_assoc_doublets
  bothLst[[ck]] = both_assoc_doublets
  
}

rm(tmpso)


sapply(pegasusLst,length)
sapply(scrubletLst,length)
sapply(eitherLst,length)
sapply(bothLst,length)


saveRDS(
  pegasusLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversPegasusAssocDoublets.RDS"
)

saveRDS(
  scrubletLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletAssocDoublets.RDS"
)

saveRDS(
  eitherLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletORPegasusAssocDoublets.RDS"
)

saveRDS(
  bothLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversScrubletANDPegasusAssocDoublets.RDS"
)





