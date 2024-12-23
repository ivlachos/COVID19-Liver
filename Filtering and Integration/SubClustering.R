rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(rafalib)
library(cluster)
library(clValid)


# ==== Functions =====
ClusterHelperFun = function(tissue_so, niter = 10, cell_fraction = 0.1, resolution_grid = seq(0.01,0.05,0.01),RecomputeSNN=TRUE){
  # subsample cells
  set.seed(17)
  rndCells = sample(x=Cells(tissue_so),size=ceiling( length(Cells(tissue_so))*cell_fraction ))
  tmpso = tissue_so[,rndCells]
  
  if(RecomputeSNN){
    
    # # Recompute SNN graph using Harmony embedding
    # tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:30)
    # 
    tmpso <- FindVariableFeatures(tmpso, selection.method = "vst", nfeatures = 2100)
    
    # PCA
    # scale all genes
    all.genes <- rownames(tmpso)
    tmpso <- ScaleData(tmpso, features = all.genes)
    
    # PCA
    hvg =  setdiff(VariableFeatures(object = tmpso),mt_genes)
    tmpso <- RunPCA(tmpso, features = hvg[1:2000])
    
    # harmony
    tmpso = RunHarmony(tmpso,"Batch", plot_convergence = TRUE,assay.use="RNA")
    tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:30)
  }
  
  
  # Get distance based on Harmony 
  subset_dist = dist(x=Embeddings(object=tmpso[["harmony"]])[,1:30])
  
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
    # SubClusters = tmpso$seurat_clusters,
    final_resolution = final_resolution
  )
  
  return(resLst)
}

# ==== Liver Data ====
liver_all = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAll.RDS"
)

mt_genes = grep("^MT-",rownames(liver_all),value=T)

# ==== Processing ====
# Filter high MT% Cells
liver_all = liver_all[,liver_all$percent_mito <= 20]

# Normalization
liver_all <- NormalizeData(object = liver_all)
liver_all <- FindVariableFeatures(liver_all, selection.method = "vst", nfeatures = 2100)

# PCA
# scale all genes
all.genes <- rownames(liver_all)
liver_all <- ScaleData(liver_all, features = all.genes)

# PCA
hvg =  setdiff(VariableFeatures(object = liver_all),mt_genes)
liver_all <- RunPCA(liver_all, features = hvg)#VariableFeatures(object = liver_all))

# harmony
liver_all = RunHarmony(liver_all,"Batch", plot_convergence = TRUE,assay.use="RNA")
liver_all = FindNeighbors(object=liver_all, reduction = "harmony",dims=1:30)
liver_all = RunUMAP(object=liver_all, reduction = "harmony",dims=1:30)

DimPlot(
  liver_all,
  group.by = "Source"
)

# Initial Clustering
InitialClust = ClusterHelperFun(
  tissue_so=liver_all, 
  niter = 10, 
  cell_fraction = 0.1, 
  RecomputeSNN = TRUE,
  resolution_grid = seq(0.001,0.01,0.001)
)

MyColor = rep("black",length(InitialClust$silhouetteLst))
MyColor[names(InitialClust$silhouetteLst) == InitialClust$final_resolution] = "red"
mypar(3,1)
boxplot(InitialClust$silhouetteLst,main="Silhoutte Score",border=MyColor)
boxplot(InitialClust$dunnLst,main="Dunn Index",border=MyColor)
boxplot(InitialClust$clustSizeLst,main="Number of Clusters",border=MyColor)


InitialClust$clustSizeLst[[as.character(InitialClust$final_resolution)]]

SelectedResolution = 0.01

liver_all = FindClusters(liver_all,resolution=SelectedResolution,random.seed=665)

DimPlot(
  liver_all,
  group.by = "seurat_clusters",
  label=T,
  repel=T
)

saveRDS(
  liver_all,
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/Clusters/CumulusCellbenderLiverAllMtThr.RDS"
)


# ==== SubClustering ===
CurrentCluster = "4"
# Subset data
tmpso = liver_all[,liver_all$seurat_clusters == CurrentCluster]

# Normalization
tmpso <- NormalizeData(object = tmpso)
tmpso <- FindVariableFeatures(tmpso, selection.method = "vst", nfeatures = 2100)

# PCA
# scale all genes
all.genes <- rownames(tmpso)
tmpso <- ScaleData(tmpso, features = all.genes)

# PCA
hvg =  setdiff(VariableFeatures(object = tmpso),mt_genes)
tmpso <- RunPCA(tmpso, features = hvg[1:2000])

# harmony
tmpso = RunHarmony(tmpso,"Batch", plot_convergence = TRUE,assay.use="RNA")
tmpso = FindNeighbors(object=tmpso, reduction = "harmony",dims=1:30)
tmpso = RunUMAP(object=tmpso, reduction = "harmony",dims=1:30)

# Select clustering parameters
tmpLst = ClusterHelperFun(
  tissue_so=tmpso, 
  niter = 10, 
  cell_fraction = 1, 
  RecomputeSNN = TRUE,
  resolution_grid = seq(0.01,0.1,0.01)
)


MyColor = rep("black",length(tmpLst$silhouetteLst))
MyColor[names(tmpLst$silhouetteLst) == tmpLst$final_resolution] = "red"
mypar(3,1)
boxplot(tmpLst$silhouetteLst,main="Silhoutte Score",border=MyColor)
boxplot(tmpLst$dunnLst,main="Dunn Index",border=MyColor)
boxplot(tmpLst$clustSizeLst,main="Number of Clusters",border=MyColor)


tmpLst$clustSizeLst[[as.character(tmpLst$final_resolution)]]

# 0: 0.3
# 1: 0.01
# 2: 0.05
# 3: 0.03
# 4: 0.08
SelectedResolution = 0.08

tmpso = FindClusters(tmpso,resolution=SelectedResolution)


DimPlot(
  tmpso,
  group.by = "seurat_clusters"
)

saveRDS(
  tmpso,
  paste0(
    "/data/work/Projects/BrScRNAseq/data/BroadTerra/Clusters/CumulusCellbenderLiverAll",
    "Cluster",CurrentCluster,".RDS"
  )
)

rm(tmpso)



