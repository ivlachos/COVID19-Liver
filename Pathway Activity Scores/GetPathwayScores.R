rm(list=ls())
options(stringsAsFactors = FALSE)


library(Matrix)
library(parallel)
library(limma)
library(tidyverse)



# ==== Functions =====
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

mmec = c(
  "PBK",
  "SMC2",
  "RACGAP1",
  "UBE2C",
  "CKS2",
  "CCNB2",
  "TOP2A",
  "LMNB1",
  "CENPA",
  "PRC1",
  "NUSAP1",
  "SMC4",
  "RRM2",
  "SPC25",
  "CENPF",
  "CDK1",
  "CDC20",
  "BIRC5",
  "PIMREG"
)

gsLst = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathways.RDS"
)

gsLst = lapply(gsLst,unique)

# Sanchez-Taltavull Signature
gsLst[["KC Proliferation"]] = tmpGenes
# Niethammer Signature
gsLst[["Proliferative EC"]] = mmec

gsDF = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathwaysDF.RDS"
)




# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1_DonorIDs.RDS"
)
liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)



# ===== Adjusted Ranks ====
comp_empAdjZRnks = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
    "AdjustedEmpiricalZRanks.RDS"
  )
)

comp_empAdjZRnks = comp_empAdjZRnks[,rownames(liver_metadata)]

# ===== Pathway Scores ======
gsLst = lapply(
  gsLst,
  function(x){x[x %in% rownames(comp_empAdjZRnks)]}
)
gsLst = gsLst[sapply(gsLst,length) > 0]


comp_pathwayScores = mclapply(
  gsLst,
  function(x){
    GetPathwayScore(zrnkMat = comp_empAdjZRnks, gs = x)
  },
  mc.cores=8
)
comp_pathwayScores = do.call(rbind,comp_pathwayScores)

rm(comp_empAdjZRnks)

# adjust pathway scores for batch
comp_pathwayScoresAdj = removeBatchEffect(
  x=comp_pathwayScores,
  batch = liver_metadata[colnames(comp_pathwayScores),"Batch"],
  covariates = liver_metadata[colnames(comp_pathwayScores),c("nGene","lognUMI")]
)

rm(comp_pathwayScores)

saveRDS(
  comp_pathwayScoresAdj,
  "/data/shiny/shiny-server/apps/SelectLiverNoDblPathways/data/AdjustedAllPathwayScoresAdjEmpZRanksUpdate.RDS"
)




