rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(patchwork)
library(uwot)
library(rafalib)


# ===== Liver Data ====
liver_tmp = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)
liver_tmp[["DevianceResiduals"]] = NULL

umapLst = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVIDPlus/data/SelectLiverNoDoubletsUMAPembeddingsLst.RDS"
)


# ===== Compartment 0 ====
CurrentCluster="0"
# Subset data
tmpso = liver_tmp[,liver_tmp$seurat_clusters == CurrentCluster]
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

tmpso = RunUMAP(
  tmpso,
  reduction = "harmony", 
  dims = 1:15,
  umap.method = "uwot",
  return.model = T
)

p01 = DimPlot(tmpso)

p00 = ggplot(
  data.frame(
    umapLst[["0"]]
  ),
  aes(x=UMAP_1,y=UMAP_2)
) +
  geom_point(size=0.25)

p00 + p01



saveRDS(
  tmpso,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets_Compartment0UWOT.RDS"
)



# ===== Compartment 1 ====
CurrentCluster="1"
# Subset data
tmpso = liver_tmp[,liver_tmp$seurat_clusters == CurrentCluster]
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

tmpso = RunUMAP(
  tmpso,
  reduction = "harmony", 
  dims = 1:15,
  umap.method = "uwot",
  return.model = T
)

p01 = DimPlot(tmpso)

p00 = ggplot(
  data.frame(
    umapLst[["1"]]
  ),
  aes(x=UMAP_1,y=UMAP_2)
) +
  geom_point(size=0.25)

p00 + p01

saveRDS(
  tmpso,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets_Compartment1UWOT.RDS"
)

# ===== Compartment 2 ====
CurrentCluster="2"
# Subset data
tmpso = liver_tmp[,liver_tmp$seurat_clusters == CurrentCluster]
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

tmpso = RunUMAP(
  tmpso,
  reduction = "harmony", 
  dims = 1:15,
  umap.method = "uwot",
  return.model = T
)

p01 = DimPlot(tmpso)

p00 = ggplot(
  data.frame(
    umapLst[[CurrentCluster]]
  ),
  aes(x=UMAP_1,y=UMAP_2)
) +
  geom_point(size=0.25)

p00 + p01

saveRDS(
  tmpso,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets_Compartment2UWOT.RDS"
)



# ===== Compartment 3 ====
CurrentCluster="3"
# Subset data
tmpso = liver_tmp[,liver_tmp$seurat_clusters == CurrentCluster]
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

tmpso = RunUMAP(
  tmpso,
  reduction = "harmony", 
  dims = 1:15,
  umap.method = "uwot",
  return.model = T
)

p01 = DimPlot(tmpso)

p00 = ggplot(
  data.frame(
    umapLst[[CurrentCluster]]
  ),
  aes(x=UMAP_1,y=UMAP_2)
) +
  geom_point(size=0.25)

p00 + p01

saveRDS(
  tmpso,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets_Compartment3UWOT.RDS"
)




# ===== Compartment 4 ====
CurrentCluster="4"
# Subset data
tmpso = liver_tmp[,liver_tmp$seurat_clusters == CurrentCluster]
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

tmpso = RunUMAP(
  tmpso,
  reduction = "harmony", 
  dims = 1:15,
  umap.method = "uwot",
  return.model = T
)

p01 = DimPlot(tmpso)

p00 = ggplot(
  data.frame(
    umapLst[[CurrentCluster]]
  ),
  aes(x=UMAP_1,y=UMAP_2)
) +
  geom_point(size=0.25)

p00 + p01

saveRDS(
  tmpso,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets_Compartment4UWOT.RDS"
)


