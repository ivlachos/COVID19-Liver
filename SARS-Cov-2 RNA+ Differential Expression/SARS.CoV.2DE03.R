rm(list=ls())
options(stringsAsFactors = FALSE)


library(Matrix)
library(rafalib)
library(parallel)
library(limma)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(ggrepel)

library(googledrive)
library(RColorBrewer)
library(DT)
Sys.setenv("RSTUDIO_PANDOC" = "/data/resources/tools/pandoc/")

# ==== Pathways ====
tmpGenes = toupper(c(
  'cdk1',
  'pcna',
  'e2f1',
  'pold3',
  'mki67',
  'orc2',
  'rrm1',
  'cdt1',
  'cdc6',
  'mcm5',
  'rrm2',
  'orc1',
  'ticrr'
))


gsLst = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathways.RDS"
)

gsLst = lapply(gsLst,unique)

# Sanchez-Taltavull Signature
gsLst[["KC Proliferation"]] = tmpGenes

gsDF = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathwaysDF.RDS"
)

library(readxl)

# ==== Cluster Manual Annotation ====
subcluster_annot = read_xlsx(
  path="/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/subcluster_annotation.xlsx",
  sheet=1
) %>% as.data.frame()

# cluster names used in manuscript
short_names = subcluster_annot$Label
names(short_names) = subcluster_annot$`Clustering v.1.1`

# ==== Liver Data ====
# COVID+ Boot Test Results
tmp00 = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_BootCovidTest.RDS"
)
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1.RDS"
)
liver_metadata$Boot.SARS.CoV.2 = tmp00[rownames(liver_metadata),"Boot.SARS.CoV.2"]
rm(tmp00)

liver_metadata$Batch2 = gsub("ncLab","",liver_metadata$Batch)
liver_metadata$Batch2 = gsub("Broad","",liver_metadata$Batch2)

liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)

covid_umi_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsCOVIDcounts.RDS"
)

liver_metadata$SARS.CoV.2.UMI = covid_umi_counts[,rownames(liver_metadata)]
rm(covid_umi_counts)

liver_metadata$COVID.plus = liver_metadata$Boot.SARS.CoV.2 == "RNA.Plus" & liver_metadata$SARS.CoV.2.UMI >= 2

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2a/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)



# ==== Nuclei Selection ====
# at least 10 COVID+ nuclei
rnaplus_sum = tapply(X=liver_metadata$COVID.plus, INDEX= liver_metadata$SubCluster, FUN =sum)
selected_subck = names(which(rnaplus_sum >= 10))


MrkUpLst = list()
MrkDownLst = list()
PathUpLst = list()
PathDownLst = list()

for(sck in  selected_subck){
  # ====== Sars-Cov-2 RNA+ =====
  subcluster_metadata = liver_metadata[liver_metadata$SubCluster == sck,]
  rnaplus_nuclei = rownames(subcluster_metadata)[subcluster_metadata$COVID.plus]
  # keep only COVID+ and matched COVID-
  subcluster_metadata = subcluster_metadata[rnaplus_nuclei,]
  subcluster_counts = liver_counts[,rownames(subcluster_metadata)]
  subcluster_counts = subcluster_counts[-1*grep("sars",rownames(subcluster_counts),ignore.case = T,value=F),]
  
  
  
  # ===== Expression Summaries =====
  # - percent expressed
  RNAMinusPctExprs = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_PctEsxprs.RDS"
    )
  )
  RNAPlusPctExprs = 100*Matrix::rowMeans(subcluster_counts > 0)
  
  CountsPctExprs = data.frame(
    RNA.plus = RNAPlusPctExprs,
    RNA.minus = RNAMinusPctExprs[names(RNAPlusPctExprs)]
  )
  
  # ===== Pathway Summaries =====
  PathwayPctExprs = mclapply(
    gsLst,
    function(x){
      colMeans(CountsPctExprs[rownames(CountsPctExprs) %in% x,,drop=F])
    },
    mc.cores = 12
  )
  PathwayPctExprs = do.call(rbind,PathwayPctExprs)
  
  
  
  # ===== Markers Up =====
  MarkerUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_MarkersUp.RDS"
    )
  )
  
  MarkerUpDF = MarkerUpDF[MarkerUpDF$AvgDiff  > 0 & MarkerUpDF$FDR < 0.05,]
  
  if(nrow(MarkerUpDF) > 0){
    COVIDPlusMarkersUp = data.frame(
      Gene = MarkerUpDF$Gene,
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = MarkerUpDF$AvgDiff,
      Pval = MarkerUpDF$Pval,
      AdjPval = MarkerUpDF$FDR,
      Avg.RNA.Plus = MarkerUpDF$RNA.plus,
      Avg.RNA.Minus = MarkerUpDF$RNA.minus,
      PctExprs = CountsPctExprs[MarkerUpDF$Gene,]
    )
    COVIDPlusMarkersUp = COVIDPlusMarkersUp[order(abs(COVIDPlusMarkersUp$LogFC),decreasing = T),]
    COVIDPlusMarkersUp = COVIDPlusMarkersUp[apply(COVIDPlusMarkersUp[,c("PctExprs.RNA.plus","PctExprs.RNA.minus")],1,function(x){any(x >= 5)}),]
    rownames(COVIDPlusMarkersUp) = NULL
    
    MrkUpLst[[sck]] = COVIDPlusMarkersUp
  }
  
  # ===== Markers Down =====
  MarkerDownDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_MarkersDown.RDS"
    )
  )
  
  MarkerDownDF = MarkerDownDF[MarkerDownDF$AvgDiff  < 0 & MarkerDownDF$FDR < 0.05,,drop=F]
  
  if(nrow(MarkerDownDF) > 0 ){
    COVIDPlusMarkersDown = data.frame(
      Gene = MarkerDownDF$Gene,
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = MarkerDownDF$AvgDiff,
      Pval = MarkerDownDF$Pval,
      AdjPval = MarkerDownDF$FDR,
      Avg.RNA.Plus = MarkerDownDF$RNA.plus,
      Avg.RNA.Minus = MarkerDownDF$RNA.minus,
      PctExprs = CountsPctExprs[MarkerDownDF$Gene,]
    )
    COVIDPlusMarkersDown = COVIDPlusMarkersDown[order(abs(COVIDPlusMarkersDown$LogFC),decreasing = T),]
    COVIDPlusMarkersDown = COVIDPlusMarkersDown[apply(COVIDPlusMarkersDown[,c("PctExprs.RNA.plus","PctExprs.RNA.minus")],1,function(x){any(x >= 5)}),]
    rownames(COVIDPlusMarkersDown) = NULL
    
    MrkDownLst[[sck]] = COVIDPlusMarkersDown
  }
  
  # ===== Pathways Up =====
  PathwayUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_PathwaysUp.RDS"
    )
  )
  
  PathwayUpDF = PathwayUpDF[PathwayUpDF$AvgDiff  > 0 & PathwayUpDF$FDR < 0.05,]
  
  if(nrow(PathwayUpDF) > 0){
    COVIDPlusPathwaysUp = data.frame(
      Pathway = PathwayUpDF$Pathway,
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = PathwayUpDF$AvgDiff,
      Pval = PathwayUpDF$Pval,
      AdjPval = PathwayUpDF$FDR,
      Avg.RNA.Plus = PathwayUpDF$RNA.plus,
      Avg.RNA.Minus = PathwayUpDF$RNA.minus,
      PctExprs.RNA.Plus = PathwayPctExprs[PathwayUpDF$Pathway,"RNA.plus"],
      PctExprs.RNA.Minus = PathwayPctExprs[PathwayUpDF$Pathway,"RNA.minus"]
    )
    COVIDPlusPathwaysUp = COVIDPlusPathwaysUp[order(abs(COVIDPlusPathwaysUp$LogFC),decreasing = T),]
    rownames(COVIDPlusPathwaysUp) = NULL
    
    PathUpLst[[sck]] = COVIDPlusPathwaysUp
    
  }
  
  
  # ===== Pathways Down =====
  PathwayDownDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_PathwaysDown.RDS"
    )
  )
  
  PathwayDownDF = PathwayDownDF[PathwayDownDF$AvgDiff  < 0 & PathwayDownDF$FDR < 0.05,,drop=F]
  
  if(nrow(PathwayDownDF) > 0 ){
    COVIDPlusPathwaysDown = data.frame(
      Pathway = PathwayDownDF$Pathway,
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = PathwayDownDF$AvgDiff,
      Pval = PathwayDownDF$Pval,
      AdjPval = PathwayDownDF$FDR,
      Avg.RNA.Plus = PathwayDownDF$RNA.plus,
      Avg.RNA.Minus = PathwayDownDF$RNA.minus,
      PctExprs.RNA.Plus = PathwayPctExprs[PathwayDownDF$Pathway,"RNA.plus"],
      PctExprs.RNA.Minus = PathwayPctExprs[PathwayDownDF$Pathway,"RNA.minus"]
    )
    COVIDPlusPathwaysDown = COVIDPlusPathwaysDown[order(abs(COVIDPlusPathwaysDown$LogFC),decreasing = T),]
    rownames(COVIDPlusPathwaysDown) = NULL
    
    PathDownLst[[sck]] = COVIDPlusPathwaysDown
  }
}


# ===== Aggregate Results ====
dribblePaths = c(
  "0" = "https://drive.google.com/drive/u/0/folders/1At8Uj1MPF7AbVBXMlRleUjpgpZOl1ffO",
  "1" = "https://drive.google.com/drive/u/0/folders/1rcBk5K-Wew2_eCRvU4fI_es8O7CnFKXc"
)

clustLst = split(
  x=selected_subck,
  f=sapply(strsplit(selected_subck,"_"),head,1)
)

for(ck in names(clustLst)){
  
  MarkersUpDF = do.call(
    rbind,
    MrkUpLst[clustLst[[ck]]]
  )
  rownames(MarkersUpDF) = NULL
  
  MarkersDownDF = do.call(
    rbind,
    MrkDownLst[clustLst[[ck]]]
  )
  rownames(MarkersDownDF) = NULL
  
  PathwaysUpDF = do.call(
    rbind,
    PathUpLst[clustLst[[ck]]]
  )
  rownames(PathwaysUpDF) = NULL
  
  PathwaysDownDF = do.call(
    rbind,
    PathDownLst[clustLst[[ck]]]
  )
  rownames(PathwaysDownDF) = NULL
  
  
  
  # Upload results
  write_tsv(
    x=MarkersUpDF,
    file=paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_MarkersUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_MarkersUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  write_tsv(
    x=MarkersDownDF,
    file=paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_MarkersDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_MarkersDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  write_tsv(
    x=PathwaysUpDF,
    file=paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_PathwaysUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_PathwaysUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  write_tsv(
    x=PathwaysDownDF,
    file=paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_PathwaysDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "Compartment",ck,"_PathwaysDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}


# get group coefficients
geneCoeffLst = list()
pathwayCoeffLst = list()
for(sck in selected_subck){
  
  # get scaled coefficients
  MarkerUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_MarkersUp.RDS"
    )
  )
  
  subclust_gene_coeff = MarkerUpDF[,c("RNA.plus", "RNA.minus")]
  coeff_avg = mean(unlist(subclust_gene_coeff))
  coeff_sd = sd(unlist(subclust_gene_coeff))
  subclust_gene_coeff = (subclust_gene_coeff-coeff_avg)/coeff_sd
  rownames(subclust_gene_coeff) = MarkerUpDF$Gene
  
  geneCoeffLst[[sck]] = subclust_gene_coeff
  
  
  PathwayUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/",
      "SubCluster",sck,"_PathwaysUp.RDS"
    )
  )
  
  subclust_path_coeff = PathwayUpDF[,c("RNA.plus", "RNA.minus")]
  coeff_avg = mean(unlist(subclust_path_coeff))
  coeff_sd = sd(unlist(subclust_path_coeff))
  subclust_path_coeff = (subclust_path_coeff-coeff_avg)/coeff_sd
  rownames(subclust_path_coeff) = PathwayUpDF$Pathway
  
  pathwayCoeffLst[[sck]] = subclust_path_coeff
}







