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

library(scater)
library(sva)
library(iasva)

library(RColorBrewer)
library(DT)
Sys.setenv("RSTUDIO_PANDOC" = "/data/resources/tools/pandoc/")

library(googledrive)


# ==== Functions =====
GetEmpiricalZRanks = function(ExprsMat,num_cores=12){
  zrnk = mclapply(
    1:ncol(ExprsMat),
    function(ic){
      # get rank
      ri = rank(ExprsMat[,ic],ties.method = "min")
      # scale and center
      c(scale(ri))
    },
    mc.cores = num_cores
  )
  
  zrnk = do.call(cbind,zrnk)
  
  dimnames(zrnk) = dimnames(ExprsMat)
  
  return(zrnk)
}

GetPathwayScore = function(zrnkMat,gs){
  # make sure all genes are present
  gs = gs[gs %in% rownames(zrnkMat)]
  ng = length(gs)
  
  if(ng > 1){
    pathway_score = matrixStats::colMeans2(zrnkMat[gs,])*sqrt(ng)
    names(pathway_score) = colnames(zrnkMat)
  }else if(ng == 1){
    pathway_score = zrnkMat[gs,]
  }else{
    pathway_score = rep(NA,ncol(zrnkMat))
    names(pathway_score) = colnames(zrnkMat)
  }
  return(pathway_score)
}

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


# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1.RDS"
)
liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2a/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)

liver_zrnk = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
    "AdjustedEmpiricalZRanks.RDS"
  )
)

liver_zrnk = liver_zrnk[,rownames(liver_metadata)]
# ===== Toronto Livers ====
toronto_counts = readRDS(
  "/data/work/Projects/BrScRNAseq/data/Toronto/Healthy_Liver_Toronto/toronto_umi_counts.RDS"
)

toronto_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/data/Toronto/Healthy_Liver_Toronto/toronto_full_metadata.RDS"
)


toronto_metadata$nUMI = toronto_metadata$nCount_RNA
toronto_metadata$nGene = toronto_metadata$nFeature_RNA

toronto_metadata$lognUMI = log(toronto_metadata$nUMI)
toronto_metadata$lognGene = log(toronto_metadata$nGene)

toronto_predictions = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "Toronto_SecretJointSeuratObject_ControlPred.RDS"
  )
)

toronto_metadata = cbind(toronto_predictions[rownames(toronto_predictions), ],toronto_predictions)
toronto_counts = toronto_counts[,rownames(toronto_metadata)]

toronto_metadata$batch = sapply(strsplit(toronto_metadata$sample,"_"),tail,1)
toronto_metadata$donor = sapply(strsplit(toronto_metadata$sample,"_"),head,1)

toronto_zrnk = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Toronto/ExpressionRanksEmpirical/",
    "TorontoAdjustedEmpiricalZRanks.RDS"
  )
)

toronto_zrnk = toronto_zrnk[,rownames(toronto_metadata)]

# at least 5 cells per sample, at least 3 samples
tmptab = colSums(apply(table(toronto_metadata[,c("SubCluster","sample")]),1,function(x){x >= 5}))
AllSubCusters = names(which(tmptab >= 3))
rm(tmptab)

current_obj = c(ls(),"current_obj")


for(sck in AllSubCusters){
  cat(
    paste0(
      "\n---------------------\n",
      "Cluster:\t",sck,"\n",
      "---------------------\n"
    )
  )
  
  cat(
    paste0(
      "\n[",sck,"] Pseudobulking samples ...\n"
    )
  )
  # ====PseudoBulk Approach ====
  # subset by cluster
  toronto_pseudobulk = scater::sumCountsAcrossCells(
    x=toronto_counts[,rownames(toronto_metadata)[toronto_metadata$SubCluster == sck]],
    ids = toronto_metadata$sample[toronto_metadata$SubCluster == sck],
    average=FALSE
  )
  liver_pseudobulk = scater::sumCountsAcrossCells(
    x=liver_counts[,rownames(liver_metadata)[liver_metadata$SubCluster == sck]],
    ids = liver_metadata$Batch[liver_metadata$SubCluster == sck],
    average=FALSE
  )
  
  common_genes = intersect(rownames(toronto_pseudobulk),rownames(liver_pseudobulk))
  
  pseudobulk_counts = cbind(
    toronto_pseudobulk[common_genes,],
    liver_pseudobulk[common_genes,]
  )
  
  pseudobulk_metadata = data.frame(
    Condition = ifelse(
      colnames(pseudobulk_counts) %in% colnames(liver_pseudobulk),
      "COVID","Control"
    ),
    Batch = liver_metadata$Source[match(colnames(pseudobulk_counts),liver_metadata$Batch)]
  )
  
  tmpind = match(colnames(pseudobulk_counts),toronto_metadata$sample)
  pseudobulk_metadata$Batch[ !is.na(tmpind) ] = toronto_metadata$batch[tmpind[!is.na(tmpind)]]
  rm(tmpind)
  rownames(pseudobulk_metadata) = colnames(pseudobulk_counts)
  
  
  cat(
    paste0(
      "\n[",sck,"] Fitting pseudobulk linear model ...\n"
    )
  )
  # ==== Pseudobulk Linear Model ====
  pseudobulk_dge = DGEList(
    counts = pseudobulk_counts
  )
  pseudobulk_dge = calcNormFactors(pseudobulk_dge)
  
  
  pseudobulk_logcpm = edgeR::cpm(pseudobulk_dge,log=T)
  
  pseudobulk_se <- SummarizedExperiment(assays = pseudobulk_logcpm)
  pseudobulk_model = model.matrix( ~ Condition,pseudobulk_metadata)
  pseudobulk_nsv <- num.sv(pseudobulk_logcpm,pseudobulk_model,seed=42)
  pseudobulk_iasvobj = iasva(pseudobulk_se, pseudobulk_model,verbose = FALSE, 
                             permute = FALSE, num.sv = pseudobulk_nsv)
  pseudobulk_sv = pseudobulk_iasvobj$sv
  rownames(pseudobulk_sv) = colnames(pseudobulk_logcpm)
  pseudobulk_sv_design = model.matrix(
    ~ Condition + pseudobulk_sv,pseudobulk_metadata
  )
  
  pseudobulk_voom = voom(pseudobulk_dge, pseudobulk_sv_design, plot = T)
  pseudobulk_lfit = lmFit(object=pseudobulk_voom,design=pseudobulk_sv_design)
  pseudobulk_lfit = eBayes(pseudobulk_lfit)
  pseudobulk_toptab = topTable(fit = pseudobulk_lfit,coef =  "ConditionCOVID", number=Inf)
  
  
  cat(
    paste0(
      "\n[",sck,"] Pseudobulking expression ranks ...\n"
    )
  )
  # ==== Ranks ====
  toronto_zrnk_pseudobulk = scater::sumCountsAcrossCells(
    x=toronto_zrnk[,rownames(toronto_metadata)[toronto_metadata$SubCluster == sck]],
    ids = toronto_metadata$sample[toronto_metadata$SubCluster == sck],
    average=TRUE
  )
  liver_zrnk_pseudobulk = scater::sumCountsAcrossCells(
    x=liver_zrnk[,rownames(liver_metadata)[liver_metadata$SubCluster == sck]],
    ids = liver_metadata$Batch[liver_metadata$SubCluster == sck],
    average=TRUE
  )
  
  
  pseudobulk_zrnk = cbind(
    toronto_zrnk_pseudobulk[common_genes,],
    liver_zrnk_pseudobulk[common_genes,]
  )
  
  cat(
    paste0(
      "\n[",sck,"] Estimating pathway scores ...\n"
    )
  )
  # ==== Pathway Scores ====
  gsLst = lapply(
    gsLst,
    function(x){x[x %in% rownames(pseudobulk_zrnk)]}
  )
  gsLst = gsLst[sapply(gsLst,length) > 0]
  
  go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
  go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])
  
  comp_pathwayScores = mclapply(
    gsLst,
    function(x){
      GetPathwayScore(zrnkMat = pseudobulk_zrnk, gs = x)
    },
    mc.cores=8
  )
  comp_pathwayScores = do.call(rbind,comp_pathwayScores)
  
  comp_pathwayScores = comp_pathwayScores[,rownames(pseudobulk_metadata)]
  
  
  cat(
    paste0(
      "\n[",sck,"] Estimating pathway linear model ...\n"
    )
  )
  pathwayScore_se <- SummarizedExperiment(assays = comp_pathwayScores)
  pathwayScore_model = model.matrix( ~ Condition,pseudobulk_metadata)
  pathwayScore_nsv <- num.sv(comp_pathwayScores,pathwayScore_model,seed=42)
  pathwayScore_iasvobj = iasva(pathwayScore_se, pathwayScore_model,verbose = FALSE, 
                               permute = FALSE, num.sv = pathwayScore_nsv)
  pathwayScore_sv = pathwayScore_iasvobj$sv
  rownames(pathwayScore_sv) = colnames(comp_pathwayScores)
  pathwayScore_sv_design = model.matrix(
    ~ Condition + pathwayScore_sv,pseudobulk_metadata
  )
  
  pathway_lfit = lmFit( object = comp_pathwayScores, design = pathwayScore_sv_design)
  pathway_lfit = eBayes(pathway_lfit)
  pathway_toptab = topTable(fit=pathway_lfit,coef="ConditionCOVID",number=Inf)
  
  
  cat(
    paste0(
      "\n[",sck,"] Estimating expression summaries ...\n"
    )
  )
  # ==== Expression Summaries ====
  # - percent expressed
  CountsPctExprs = data.frame(
    Control = 100*Matrix::rowMeans(toronto_counts[common_genes,rownames(toronto_metadata)[toronto_metadata$SubCluster == sck]] > 0),
    COVID = 100*Matrix::rowMeans(liver_counts[common_genes,rownames(liver_metadata)[liver_metadata$SubCluster == sck]] > 0)
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
  
  
  cat(
    paste0(
      "\n[",sck,"] Formatting markers ...\n"
    )
  )
  # ==== Markers ====
  up_mrk_ind = pseudobulk_toptab$logFC > 0 & pseudobulk_toptab$adj.P.Val < 0.05
  if(any(up_mrk_ind)){
    COVIDMarkersUp = data.frame(
      Gene = rownames(pseudobulk_toptab)[ up_mrk_ind ],
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = pseudobulk_toptab$logFC[up_mrk_ind],
      Pval = pseudobulk_toptab$P.Value[up_mrk_ind],
      AdjPval = pseudobulk_toptab$adj.P.Val[up_mrk_ind],
      PctExprs = CountsPctExprs[rownames(pseudobulk_toptab)[ up_mrk_ind ],,drop=F]
    )
    rownames(COVIDMarkersUp) = NULL
    
    saveRDS(
      COVIDMarkersUp,
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/TorontoDE/",
        "Cluster_",sck,"_COVIDMarkersUp.RDS"
      )
    )
  }
  
  
  
  down_mrk_ind = pseudobulk_toptab$logFC < 0 & pseudobulk_toptab$adj.P.Val < 0.05
  if(any(down_mrk_ind)){
    
    COVIDMarkersDown = data.frame(
      Gene = rownames(pseudobulk_toptab)[ down_mrk_ind ],
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = pseudobulk_toptab$logFC[down_mrk_ind],
      Pval = pseudobulk_toptab$P.Value[down_mrk_ind],
      AdjPval = pseudobulk_toptab$adj.P.Val[down_mrk_ind],
      PctExprs = CountsPctExprs[rownames(pseudobulk_toptab)[ down_mrk_ind ],,drop=F]
    )
    rownames(COVIDMarkersDown) = NULL
    
    saveRDS(
      COVIDMarkersDown,
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/TorontoDE/",
        "Cluster_",sck,"_COVIDMarkersDown.RDS"
      )
    )
    
  }
  
  
  cat(
    paste0(
      "\n[",sck,"] Formatting pathways ...\n"
    )
  )
  
  # ==== Pathways ====
  up_path_ind = pathway_toptab$logFC > 0 & pathway_toptab$adj.P.Val < 0.05
  if(any(up_path_ind)){
    
    COVIDPathwaysUp = data.frame(
      Pathway = rownames(pathway_toptab)[ up_path_ind ],
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = pathway_toptab$logFC[up_path_ind],
      Pval = pathway_toptab$P.Value[up_path_ind],
      AdjPval = pathway_toptab$adj.P.Val[up_path_ind],
      PctExprs = PathwayPctExprs[rownames(pathway_toptab)[ up_path_ind ],,drop=F]
    )
    rownames(COVIDPathwaysUp) = NULL
    
    saveRDS(
      COVIDPathwaysUp,
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/TorontoDE/",
        "Cluster_",sck,"_COVIDPathwaysUp.RDS"
      )
    )
    
  }
  
  
  
  down_path_ind = pathway_toptab$logFC < 0 & pathway_toptab$adj.P.Val < 0.05
  if(any(down_path_ind)){
    
    COVIDPathwaysDown = data.frame(
      Pathway = rownames(pathway_toptab)[ down_path_ind ],
      Cluster = sck,
      ClusterName = short_names[sck],
      LogFC = pathway_toptab$logFC[down_path_ind],
      Pval = pathway_toptab$P.Value[down_path_ind],
      AdjPval = pathway_toptab$adj.P.Val[down_path_ind],
      PctExprs = PathwayPctExprs[rownames(pathway_toptab)[ down_path_ind ],,drop=F]
    )
    rownames(COVIDPathwaysDown) = NULL
    
    saveRDS(
      COVIDPathwaysDown,
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/TorontoDE/",
        "Cluster_",sck,"_COVIDPathwaysDown.RDS"
      )
    )
    
  }
  
  cat(
    paste0(
      "\n[",sck,"] Estimating coefficients for each group ...\n"
    )
  )
  # coefficients for each group
  pseudobulk_coeff = cbind(
    Control = pseudobulk_lfit$coefficients[,"(Intercept)"],
    COVID = rowSums(pseudobulk_lfit$coefficients[,c("(Intercept)","ConditionCOVID")])
  )
  pseudobulk_coeff = (pseudobulk_coeff - pseudobulk_toptab[rownames(pseudobulk_coeff),"AveExpr"])/sqrt(pseudobulk_lfit$s2.post)
  
  saveRDS(
    pseudobulk_coeff,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Cluster_",sck,"_PseudobulkGeneGroupCoeff.RDS"
    )
  )
  
  pathway_coeff = cbind(
    Control = pathway_lfit$coefficients[,"(Intercept)"],
    COVID = rowSums(pathway_lfit$coefficients[,c("(Intercept)","ConditionCOVID")])
  )
  pathway_coeff = (pathway_coeff - pathway_toptab[rownames(pathway_coeff),"AveExpr"])/sqrt(pathway_lfit$s2.post)
  
  saveRDS(
    pathway_coeff,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Cluster_",sck,"_PseudobulkPathwayGroupCoeff.RDS"
    )
  )
  
  
  cat(
    paste0(
      "\n[",sck,"] Done >:( \n"
    )
  )
  
  cat("\n Clearing objects ... \n")
  rm(list = setdiff(objects(),current_obj) )
  
  
}
  

# ==== Aggregate Results by Compartment =====
# split clusters by compartment
clustLst = split(
  x=AllSubCusters,
  f=sapply(strsplit(AllSubCusters,"_"),head,1)
)


# ==== Upload Results ====
dribblePaths = c(
  "0" = "https://drive.google.com/drive/u/0/folders/1fc-k7H9n2oArJHU89ngAiTqCvSqIE0gA",
  "1" = "https://drive.google.com/drive/u/0/folders/1nV-F1p86qf7ytVywwBLP2lfUUbepWr0_",
  "2" = "https://drive.google.com/drive/u/0/folders/1WxbRyaje2TN2it34-Zqn_u2kNh9b654T",
  "3" = "https://drive.google.com/drive/u/0/folders/1Cjb92oEGGQmtSlM2sT_2266RvUXB_T6Z",
  "4" = "https://drive.google.com/drive/u/0/folders/1UJ-_o_PvvhpUHT_nOjItUSdoWT-asAaQ"
)


for(ck in names(clustLst)){
  
  # Markers
  COVIDMarkersUp = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDMarkersUp.RDS"
        )
      ))
    }
  )
  COVIDMarkersUp = COVIDMarkersUp[sapply(COVIDMarkersUp,class) == "data.frame"]
  COVIDMarkersUp = do.call(rbind,COVIDMarkersUp)
  
  write_tsv(
    x = COVIDMarkersUp,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDMarkersUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDMarkersUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  COVIDMarkersDown = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDMarkersDown.RDS"
        )
      ))
    }
  )
  COVIDMarkersDown = COVIDMarkersDown[sapply(COVIDMarkersDown,class) == "data.frame"]
  COVIDMarkersDown = do.call(rbind,COVIDMarkersDown)
  
  write_tsv(
    x = COVIDMarkersDown,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDMarkersDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDMarkersDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  # Pathways
  COVIDPathwaysUp = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDPathwaysUp.RDS"
        )
      ))
    }
  )
  COVIDPathwaysUp = COVIDPathwaysUp[sapply(COVIDPathwaysUp,class) == "data.frame"]
  COVIDPathwaysUp = do.call(rbind,COVIDPathwaysUp)
  
  write_tsv(
    x = COVIDPathwaysUp,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDPathwaysUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDPathwaysUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  COVIDPathwaysDown = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDPathwaysDown.RDS"
        )
      ))
    }
  )
  COVIDPathwaysDown = COVIDPathwaysDown[sapply(COVIDPathwaysDown,class) == "data.frame"]
  COVIDPathwaysDown = do.call(rbind,COVIDPathwaysDown)
  
  write_tsv(
    x = COVIDPathwaysDown,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDPathwaysDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/",
      "Compartment_",ck,"_COVIDPathwaysDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}



  

for( ck in names(clustLst) ){
  # ==== Interactive Table ( Marker Set Up ) ====
  scaled_coeff = lapply(
    clustLst[[ck]],
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkGeneGroupCoeff.RDS"
        )
      )
    }
  )
  scaled_coeff = do.call(rbind,scaled_coeff)
  
  paletteLength <- 50
  myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
  my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})
  rm(scaled_coeff)
  
  # ==== Interactive Table (Markers Up) ====
  # Markers
  COVIDMarkersUp = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDMarkersUp.RDS"
        )
      ))
    }
  )
  names(COVIDMarkersUp) = clustLst[[ck]]
  COVIDMarkersUp = COVIDMarkersUp[sapply(COVIDMarkersUp,class) == "data.frame"]
  
  
  MarkerUpSummaryDF = lapply(
    names(COVIDMarkersUp),
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkGeneGroupCoeff.RDS"
        )
      )
      res = data.frame(
        COVIDMarkersUp[[sck]],
        scaled_coeff[COVIDMarkersUp[[sck]]$Gene,,drop=F]
      )
      rownames(res) = NULL
      return(res)
    }
  )
  
  MarkerUpSummaryDF = do.call(rbind,MarkerUpSummaryDF)
  
  font.size <- "7pt"
  MarkerUpDT = datatable(
    MarkerUpSummaryDF,
    extensions = 'FixedColumns',
    # width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 1),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(4,7:ncol(MarkerUpSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(c("Control","COVID"), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerUpDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/DT/",
      "Compartment_",ck,"_COVIDMarkersUp.html"
    )
  )
  
  
  
  # ==== Interactive Table (Markers Down) ====
  # Markers
  COVIDMarkersDown = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDMarkersDown.RDS"
        )
      ))
    }
  )
  names(COVIDMarkersDown) = clustLst[[ck]]
  COVIDMarkersDown = COVIDMarkersDown[sapply(COVIDMarkersDown,class) == "data.frame"]
  
  
  MarkerDownSummaryDF = lapply(
    names(COVIDMarkersDown),
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkGeneGroupCoeff.RDS"
        )
      )
      res = data.frame(
        COVIDMarkersDown[[sck]],
        scaled_coeff[COVIDMarkersDown[[sck]]$Gene,,drop=F]
      )
      rownames(res) = NULL
      return(res)
    }
  )
  
  MarkerDownSummaryDF = do.call(rbind,MarkerDownSummaryDF)
  
  font.size <- "7pt"
  MarkerDownDT = datatable(
    MarkerDownSummaryDF,
    extensions = 'FixedColumns',
    # width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 1),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(4,7:ncol(MarkerDownSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(c("Control","COVID"), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerDownDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/DT/",
      "Compartment_",ck,"_COVIDMarkersDown.html"
    )
  )
  
  
  # ==== Interactive Table ( Pathway Set Up ) ====
  scaled_coeff = lapply(
    clustLst[[ck]],
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkPathwayGroupCoeff.RDS"
        )
      )
    }
  )
  scaled_coeff = do.call(rbind,scaled_coeff)
  
  paletteLength <- 50
  myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
  my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})
  rm(scaled_coeff)
  
  
  # ==== Interactive Table (Pathways Up) ====
  # Pathways
  COVIDPathwaysUp = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDPathwaysUp.RDS"
        )
      ))
    }
  )
  names(COVIDPathwaysUp) = clustLst[[ck]]
  COVIDPathwaysUp = COVIDPathwaysUp[sapply(COVIDPathwaysUp,class) == "data.frame"]
  
  
  PathwayUpSummaryDF = lapply(
    names(COVIDPathwaysUp),
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkPathwayGroupCoeff.RDS"
        )
      )
      res = data.frame(
        COVIDPathwaysUp[[sck]],
        scaled_coeff[COVIDPathwaysUp[[sck]]$Pathway,,drop=F]
      )
      rownames(res) = NULL
      return(res)
    }
  )
  
  PathwayUpSummaryDF = do.call(rbind,PathwayUpSummaryDF)
  
  font.size <- "7pt"
  PathwayUpDT = datatable(
    PathwayUpSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 1),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(4,7:ncol(PathwayUpSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(c("Control","COVID"), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayUpDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/DT/",
      "Compartment_",ck,"_COVIDPathwaysUp.html"
    )
  )
  
  
  
  
  # ==== Interactive Table (Pathways Down) ====
  # Pathways
  COVIDPathwaysDown = lapply(
    clustLst[[ck]],
    function(sck){
      try(readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_COVIDPathwaysDown.RDS"
        )
      ))
    }
  )
  names(COVIDPathwaysDown) = clustLst[[ck]]
  COVIDPathwaysDown = COVIDPathwaysDown[sapply(COVIDPathwaysDown,class) == "data.frame"]
  
  
  PathwayDownSummaryDF = lapply(
    names(COVIDPathwaysDown),
    function(sck){
      scaled_coeff = readRDS(
        paste0(
          "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
          "SelectLiversNoDoublets/TorontoDE/",
          "Cluster_",sck,"_PseudobulkPathwayGroupCoeff.RDS"
        )
      )
      res = data.frame(
        COVIDPathwaysDown[[sck]],
        scaled_coeff[COVIDPathwaysDown[[sck]]$Pathway,,drop=F]
      )
      rownames(res) = NULL
      return(res)
    }
  )
  
  PathwayDownSummaryDF = do.call(rbind,PathwayDownSummaryDF)
  
  font.size <- "7pt"
  PathwayDownDT = datatable(
    PathwayDownSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 1),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(4,7:ncol(PathwayDownSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(c("Control","COVID"), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayDownDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/TorontoDE/DT/",
      "Compartment_",ck,"_COVIDPathwaysDown.html"
    )
  )
}
