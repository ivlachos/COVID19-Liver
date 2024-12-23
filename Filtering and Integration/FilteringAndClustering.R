rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(cluster)
library(clValid)


# ====> ncLab Livers <=====
# ===== Broad Cumulus Results =====
liver1_so = readRDS("/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver1/Liver1.GRCh38premrna-rna.Seurat.RDS")
liver2_so = readRDS("/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver2_and_2a/Liver2_and_2a.GRCh38premrna-rna.Seurat.RDS")
liver3_so = readRDS("/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver3/Liver3.GRCh38premrna-rna.Seurat.RDS")

# match assignment in our data
liver2a_so = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver2/Liver2.GRCh38premrna-rna.Seurat.RDS"
)
liver2b_so = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver2a/Liver2a.GRCh38premrna-rna.Seurat.RDS"
)
# add batch information
liver2a_so$orig.ident = "Liver2"
liver2a_so$Batch = "ncLabLiver2a"
liver2b_so$orig.ident = "Liver2"
liver2b_so$Batch = "ncLabLiver2b"



# add batch information
liver1_so$orig.ident =  sapply(strsplit(
  Cells(liver1_so),
  "-"
),head,1)
liver1_so$Batch = paste0("ncLab",liver1_so$orig.ident)

liver2_so$orig.ident = sapply(strsplit(
  Cells(liver2_so),
  "-"
),head,1)
liver2_so$orig.ident = dplyr::recode(
  liver2_so$orig.ident,
  "Liver" = "Liver2a",
  "Liver2" = "Liver2b"
)
liver2_so$Batch = paste0("ncLab",liver2_so$orig.ident)


liver3_so$orig.ident =  sapply(strsplit(
  Cells(liver3_so),
  "-"
),head,1)
liver3_so$Batch = paste0("ncLab",liver3_so$orig.ident)

# ====> Broad Livers <=====
# read in Seurat objects
rds_files = system(
  command = 'find /data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus -type f -name "*.RDS"',
  intern = TRUE
)

# 3pNG_BIDMC2_liver_nuclei 	12-P485759-S020-R01
# 3pNG_BIDMC1_liver_nuclei 	12-P230638-S003-R01
# 3pNG_BIDMC3_liver_nuclei 	12-P617758-S003-R01
# 3pNG_BIDMC4_liver_nuclei 	12-P852049-S003-R01
# 3pNG_BIDMC5_liver_nuclei 	12-P890292-S003-R01

soLst = lapply(
  rds_files,
  readRDS
)

liverBroadCumulus = merge(
  x=soLst[[1]],
  y=soLst[-1]
)
rm(soLst)


BroadIDs = sapply(strsplit(
  Cells(liverBroadCumulus),
  "-"
),function(x){paste(x[1:4],collapse="-")})


liverBroadCumulus$Batch = dplyr::recode(
  BroadIDs,
  "12-P485759-S020-R01"="BroadLiver2",
  "12-P230638-S003-R01"="BroadLiver1",
  "12-P617758-S003-R01"="BroadLiver3",
  "12-P852049-S003-R01"="BroadLiver4",
  "12-P890292-S003-R01"="BroadLiver5"
)

liverBroadCumulus$orig.ident = dplyr::recode(
  BroadIDs,
  "12-P485759-S020-R01"="Liver2",
  "12-P230638-S003-R01"="Liver1",
  "12-P617758-S003-R01"="Liver3",
  "12-P852049-S003-R01"="Liver4",
  "12-P890292-S003-R01"="Liver5"
)



# ====> Ben Izar <====
liver_02_cov = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_02_cov_Ben/liver_02_cov_Ben.GRCh38premrna-rna.Seurat.RDS"
)
liver_12_cov = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_12_cov_Ben/liver_12_cov_Ben.GRCh38premrna-rna.Seurat.RDS"
)
liver_01_cov = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_01_cov_Ben/liver_01_cov_Ben.GRCh38premrna-rna.Seurat.RDS"
)
liver_08_cov = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_08_cov_Ben/liver_08_cov_Ben.GRCh38premrna-rna.Seurat.RDS"
)

# add batch information
liver_02_cov$orig.ident = "liver_02_cov"
liver_02_cov$Batch = "Izarliver_02_cov"

liver_12_cov$orig.ident = "liver_12_cov"
liver_12_cov$Batch = "Izarliver_12_cov"

liver_01_cov$orig.ident = "liver_01_cov"
liver_01_cov$Batch = "Izarliver_01_cov"

liver_08_cov$orig.ident = "liver_08_cov"
liver_08_cov$Batch = "Izarliver_08_cov"

# ==== Single Seurat Object ====
liver_all <- merge(
  x=liver1_so,
  y=list(
    liver2a_so,liver2b_so,
    liver3_so,
    liverBroadCumulus,
    liver_02_cov,
    liver_12_cov,
    liver_01_cov,
    liver_08_cov
  )
)

liver_all$Source = "NA"
liver_all$Source[grep("^02",liver_all$Batch,value=F)] = "MGH"
liver_all$Source[grep("^04",liver_all$Batch,value=F)] = "BWH"
liver_all$Source[grep("^Izar",liver_all$Batch,value=F)] = "Izar"
liver_all$Source[grep("^ncLab",liver_all$Batch,value=F)] = "ncLabBI"
liver_all$Source[grep("^Broad",liver_all$Batch,value=F)] = "BroadBI"



# MT Content
VlnPlot(
  liver_all,
  features = "percent_mito",
  group.by="Source",
  pt.size=0.045
)

# ==== Normalization ====
liver_all <- NormalizeData(object = liver_all)

# ==== Highly Variable Genes ====
liver_all <- FindVariableFeatures(liver_all, selection.method = "vst", nfeatures = 2000)


# ===== PCA ====
# scale all genes
all.genes <- rownames(liver_all)
liver_all <- ScaleData(liver_all, features = all.genes)
# PCA
liver_all <- RunPCA(liver_all, features = VariableFeatures(object = liver_all))

ElbowPlot(liver_all,ndims=50)



# ==== Dimensionality Reduction (UMAP) ====
liver_all <- FindNeighbors(liver_all, dims = 1:30)
liver_all <- RunUMAP(liver_all, dims = 1:30)

beforeUmap00 = DimPlot(
  liver_all,
  group.by="Source",
  label=F,
  repel=F
)

beforeUmap01 = DimPlot(
  liver_all,
  group.by="Batch",
  label=T,
  repel=T
)

beforeUmap02 = DimPlot(
  liver_all,
  group.by="orig.ident",
  label=T,
  repel=T
)

# ==== Harmony ====
liver_all = RunHarmony(liver_all,"Batch", plot_convergence = TRUE,assay.use="RNA")
liver_all = FindNeighbors(object=liver_all, reduction = "harmony",dims=1:30)
liver_all = RunUMAP(object=liver_all, reduction = "harmony",dims=1:30)


afterUmap00 = DimPlot(
  liver_all,
  group.by="Source",
  label=F,
  repel=F
)

afterUmap01 = DimPlot(
  liver_all,
  group.by="Batch",
  label=T,
  repel=T
)

afterUmap02 = DimPlot(
  liver_all,
  group.by="orig.ident",
  label=T,
  repel=T
)

DimPlot(
  liver_all,
  group.by="Batch",
  split.by="Batch",
  label=F,
  repel=F,
  ncol=6
) + NoLegend()

# ==== Clustering ====
# initial clustering
liver_all <- FindClusters(liver_all,resolution = 0.01,random.seed = 665)

DimPlot(
  liver_all,
  group.by="seurat_clusters",
  label=T,
  repel=T
)

liver_all[["percent.mt"]] <- PercentageFeatureSet(liver_all, pattern = "^MT-")

FeaturePlot(
  liver_all,
  features = "percent.mt"
)

VlnPlot(
  liver_all,
  features = "percent.mt",
  group.by="Batch",
  pt.size=0.045
)

VlnPlot(
  liver_all,
  features = "percent.mt",
  group.by="seurat_clusters",
  pt.size=0.045
)

VlnPlot(
  liver_all,
  features = "percent.mt",
  group.by="Source",
  pt.size=0.045
)

mtLst = split(liver_all$percent.mt,liver_all$Source)

sapply(mtLst,range)


# save Seurat Object with all cells
saveRDS(
  liver_all,
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAll.RDS"
)
