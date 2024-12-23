#!/usr/bin/env Rscript

rm(list=ls())
options(stringsAsFactors = FALSE)


suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rafalib))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    opt_str = c("--SUBCLUSTER"),
    action="store",
    type="integer",
    default=27,
    help="subcluster index",
    metavar = "integer"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


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


cat("\n\nLoading pathways ... \n\n")
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

rm(tmpGenes)

cat("\n\nLoading liver data ... \n\n")
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


empAdjZRnks = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
    "AdjustedEmpiricalZRanks.RDS"
  )
)

empAdjZRnks = empAdjZRnks[,rownames(liver_metadata)]

cat("\n\nNuclei selection ... \n\n")
# ==== Nuclei Selection ====
# at least 10 COVID+ nuclei
rnaplus_sum = tapply(X=liver_metadata$COVID.plus, INDEX= liver_metadata$SubCluster, FUN =sum)
selected_subck = names(which(rnaplus_sum >= 10))

# function to match nuclei
MatchRNAPlus = function(sck,ShowPlot=T,Verbose = T){
  if(Verbose){
    cat(
      paste0(
        "\n-----------------\n",
        "Cluster: ",sck,
        "\n-----------------\n\n"
      )
    )
    
    cat("\nSubsetting data ... \n\n")
  }
  
  # ==== Subset SubCluster ====
  subcluster_metadata = liver_metadata[liver_metadata$SubCluster == sck,]
  # 1) remove ambient
  table(subcluster_metadata$Boot.SARS.CoV.2)
  subcluster_metadata = subcluster_metadata[subcluster_metadata$Boot.SARS.CoV.2 != "RNA.Ambient",]
  table(subcluster_metadata$Boot.SARS.CoV.2)
  
  
  # 2) At least 10 COVID+ nuclei (checked above already)
  sum(subcluster_metadata$COVID.plus) >= 10
  
  if(Verbose){
    cat("\nSelect donors ... \n\n")
  }
  
  # 3) Donors with at least 2 COVID+ nuclei
  batch_tab = tapply(
    X=subcluster_metadata$COVID.plus,
    INDEX=subcluster_metadata$Batch,
    FUN=sum
  )
  keep_donors = names(batch_tab > 2)
  subcluster_metadata = subcluster_metadata[subcluster_metadata$Batch %in% keep_donors,]
  
  # recheck for 2)
  if(sum(subcluster_metadata$COVID.plus) >= 10){
    if(Verbose){
      cat("\nMatch complexity distribution ... \n\n")
    }
    
    # 4) Match complexity
    # Target
    rnaplus_complexity = log10(subcluster_metadata$nGene[subcluster_metadata$COVID.plus])
    # Distribution to match
    rnaminus_complexity = log10(subcluster_metadata$nGene[subcluster_metadata$Boot.SARS.CoV.2 == "RNA.Minus"])
    names(rnaminus_complexity) = rownames(subcluster_metadata)[subcluster_metadata$Boot.SARS.CoV.2 == "RNA.Minus"]
    
    # get breaks
    target_breaks = cut(rnaplus_complexity,5)
    # split target
    targetLst = split(rnaplus_complexity,target_breaks)
    
    # get bins from breaks
    bins = gsub("(","",levels(target_breaks),fixed=T)
    bins = gsub("]","",bins,fixed=T)
    bins = do.call(rbind,strsplit(bins,","))
    colnames(bins) = c("Lo","Hi")
    bins = as.data.frame(bins)
    
    # split data to match
    splitLst = lapply(
      1:nrow(bins),
      function(ic){
        names(rnaminus_complexity)[rnaminus_complexity > bins$Lo[ic] & rnaminus_complexity <= bins$Hi[ic]]
      }
    )
    # get match
    matchLst = lapply(
      seq_along(splitLst),
      function(ic){
        if(length(splitLst[[ic]]) > 0){
          sample(x=splitLst[[ic]],size=length(targetLst[[ic]]))
        }else{
          NA
        }
      }
    )
    
    
    if(ShowPlot){
      # quick plot
      mypar(3,1,mar=c(7,2,2,1))
      barplot(sapply(targetLst,length),las=2,main=sck)
      barplot(sapply(splitLst,length),las=2,main="Obs")
      barplot(sapply(matchLst,length),las=2,main="Matched")
      mypar(1,1)
    }
    
    
    matchLst = matchLst[!sapply(matchLst,function(x){any(is.na(x))})]
    if(Verbose){
      cat("DONE!\n\n")
    }
    return( unlist(matchLst))
  }else{
    return(NA)
  }
}

# ==== Setup ====
sck = selected_subck[as.numeric(opt$SUBCLUSTER)]
niter = 1E3

cat(
  paste0(
    "\n-----------------\n",
    "Cluster: ",sck,
    "\n-----------------\n\n"
  )
)

exprsLst = list()
pathLst = list()
pb = txtProgressBar(style=3,max=niter)

cat("\n\n Sars-Cov-2 RNA- \n\n")
for(ic in 1:niter){
  
  
  # selected_subck
  matched_rnaminus = MatchRNAPlus(sck=sck,ShowPlot=F,Verbose=F)
  
  
  # ==== Subset SubCluster ====
  matched_metadata = liver_metadata[matched_rnaminus,]
  
  matched_counts = liver_counts[,rownames(matched_metadata)]
  matched_counts = matched_counts[-1*grep("sars",rownames(matched_counts),ignore.case = T,value=F),]
  
  
  
  cat(paste0("\n[",ic,"] Estimating pathway scores ... \n\n"))
  # ===== Adjusted Ranks ====
  matched_empAdjZRnks = empAdjZRnks[,rownames(matched_metadata)]
  
  
  # ===== Pathway Scores =====-
  gsLst = lapply(
    gsLst,
    function(x){x[x %in% rownames(matched_empAdjZRnks)]}
  )
  gsLst = gsLst[sapply(gsLst,length) > 0]
  
  go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
  go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])
  
  matched_pathwayScores = mclapply(
    gsLst,
    function(x){
      GetPathwayScore(zrnkMat = matched_empAdjZRnks, gs = x)
    },
    mc.cores=8
  )
  matched_pathwayScores = do.call(rbind,matched_pathwayScores)
  
  rm(matched_empAdjZRnks)
  
  # adjust pathway scores for batch
  matched_pathwayScoresAdj = removeBatchEffect(
    x=matched_pathwayScores,
    batch = liver_metadata[colnames(matched_pathwayScores),"Batch"],
    covariates = liver_metadata[colnames(matched_pathwayScores),c("nGene","lognUMI")]
  )
  
  rm(matched_pathwayScores)
  
  
  cat(paste0("\n[",ic,"] Estimating average gene expression ... \n\n"))
  # ===== Gene Linear Model  =====
  # TMM Normalization
  matched_dge <- DGEList(
    counts=matched_counts[,rownames(matched_metadata)],
    group=matched_metadata$ReSubCluster
  )
  matched_dge <- calcNormFactors(object=matched_dge,method="TMM")
  matched_lcpm <- edgeR::cpm(matched_dge, log = TRUE, prior.count = 3)
  
  matched_batchAdj_lcpm = removeBatchEffect(
    x=matched_lcpm,
    batch = liver_metadata[colnames(matched_lcpm),"Batch"]
  )
  
  
  
  # store results
  exprsLst[[ic]] = rowMeans(matched_batchAdj_lcpm)
  pathLst[[ic]] = rowMeans(matched_pathwayScoresAdj)
  setTxtProgressBar(pb, ic)
  cat(paste0("\n DONE [",ic,"/",niter,"] \n\n"))
}



cat("\n\n Sars-Cov-2 RNA+ \n\n")
# ====== Sars-Cov-2 RNA+ =====
subcluster_metadata = liver_metadata[liver_metadata$SubCluster == sck,]
rnaplus_nuclei = rownames(subcluster_metadata)[subcluster_metadata$COVID.plus]
# keep only COVID+ and matched COVID-
subcluster_metadata = subcluster_metadata[rnaplus_nuclei,]

subcluster_counts = liver_counts[,rownames(subcluster_metadata)]
subcluster_counts = subcluster_counts[-1*grep("sars",rownames(subcluster_counts),ignore.case = T,value=F),]


cat(paste0("\n Estimating pathway scores ... \n\n"))
# ===== Adjusted Ranks ====
subcluster_empAdjZRnks = empAdjZRnks[,rownames(subcluster_metadata)]


# ===== Pathway Scores =====-
gsLst = lapply(
  gsLst,
  function(x){x[x %in% rownames(subcluster_empAdjZRnks)]}
)
gsLst = gsLst[sapply(gsLst,length) > 0]

go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])

subcluster_pathwayScores = mclapply(
  gsLst,
  function(x){
    GetPathwayScore(zrnkMat = subcluster_empAdjZRnks, gs = x)
  },
  mc.cores=8
)
subcluster_pathwayScores = do.call(rbind,subcluster_pathwayScores)

rm(subcluster_empAdjZRnks)

# adjust pathway scores for batch
subcluster_pathwayScoresAdj = removeBatchEffect(
  x=subcluster_pathwayScores,
  batch = liver_metadata[colnames(subcluster_pathwayScores),"Batch"],
  covariates = liver_metadata[colnames(subcluster_pathwayScores),c("nGene","lognUMI")]
)

rm(subcluster_pathwayScores)



cat(paste0("\n Estimating average gene expression ... \n\n"))
# ===== Gene Linear Model  =====
# TMM Normalization
subcluster_dge <- DGEList(
  counts=subcluster_counts[,rownames(subcluster_metadata)],
  group=subcluster_metadata$ReSubCluster
)
subcluster_dge <- calcNormFactors(object=subcluster_dge,method="TMM")
subcluster_lcpm <- edgeR::cpm(subcluster_dge, log = TRUE, prior.count = 3)

subcluster_batchAdj_lcpm = removeBatchEffect(
  x=subcluster_lcpm,
  batch = liver_metadata[colnames(subcluster_lcpm),"Batch"]
)

subcluster_exprs = rowMeans(subcluster_batchAdj_lcpm)
subcluster_path = rowMeans(subcluster_pathwayScoresAdj)

exprsMat = do.call(cbind,exprsLst)
pathMat = do.call(cbind,pathLst)


# ==== Aggregate Results ====
# estimate difference
exprdiff_boot = subcluster_exprs - rowMeans(exprsMat[names(subcluster_exprs),])
pathdiff_boot = subcluster_path - rowMeans(pathMat[names(subcluster_path),])
# pvalues (greater)
pexprg_boot = rowMeans(exprsMat[names(subcluster_exprs),] > subcluster_exprs)
ppathg_boot = rowMeans(pathMat[names(subcluster_path),] > subcluster_path)
# pvalues (less)
pexprl_boot = rowMeans(exprsMat[names(subcluster_exprs),] < subcluster_exprs)
ppathl_boot = rowMeans(pathMat[names(subcluster_path),] < subcluster_path)

MarkerUpDF = data.frame(
  Gene = names(subcluster_exprs),
  RNA.plus = subcluster_exprs,
  RNA.minus = rowMeans(exprsMat[names(subcluster_exprs),]),
  AvgDiff = exprdiff_boot,
  Pval = pexprg_boot,
  FDR = p.adjust(p=pexprg_boot,method="fdr")
)
rownames(MarkerUpDF) = NULL
head(MarkerUpDF[order(MarkerUpDF$AvgDiff,decreasing=T),],10)

MarkerDownDF = data.frame(
  Gene = names(subcluster_exprs),
  RNA.plus = subcluster_exprs,
  RNA.minus = rowMeans(exprsMat[names(subcluster_exprs),]),
  AvgDiff = exprdiff_boot,
  Pval = pexprl_boot,
  FDR = p.adjust(p=pexprl_boot,method="fdr")
)
rownames(MarkerDownDF) = NULL
head(MarkerDownDF[order(MarkerDownDF$AvgDiff,decreasing=F),],10)


PathwayUpDF = data.frame(
  Pathway = names(subcluster_path),
  RNA.plus = subcluster_path,
  RNA.minus = rowMeans(pathMat[names(subcluster_path),]),
  AvgDiff = pathdiff_boot,
  Pval = ppathg_boot,
  FDR = p.adjust(p=ppathg_boot,method="fdr")
)
PathwayUpDF = PathwayUpDF[go_pathways,]
rownames(PathwayUpDF) = NULL
head(PathwayUpDF[order(PathwayUpDF$AvgDiff,decreasing=T),],10)

PathwayDownDF = data.frame(
  Pathway = names(subcluster_path),
  RNA.plus = subcluster_path,
  RNA.minus = rowMeans(pathMat[names(subcluster_path),]),
  AvgDiff = pathdiff_boot,
  Pval = ppathl_boot,
  FDR = p.adjust(p=ppathl_boot,method="fdr")
)
PathwayDownDF = PathwayDownDF[go_pathways,]
rownames(PathwayDownDF) = NULL
head(PathwayDownDF[order(PathwayDownDF$AvgDiff,decreasing=F),],10)

sum(PathwayUpDF$FDR < 0.01)
mean(PathwayUpDF$FDR < 0.01)

sum(PathwayDownDF$FDR < 0.01)
mean(PathwayDownDF$FDR < 0.01)

sum(MarkerUpDF$FDR < 0.01)
mean(MarkerUpDF$FDR < 0.01)

sum(MarkerDownDF$FDR < 0.01)
mean(MarkerDownDF$FDR < 0.01)



cat("\n\n Saving Results \n\n")
# ===== Save Results ====
if(!dir.exists(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/"
  )
)){
  
  dir.create(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
      "SubCluster",sck,"/"
    )
  )
}

saveRDS(
  MarkerUpDF,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/",
    "SubCluster",sck,"_MarkersUp.RDS"
  )
)

saveRDS(
  MarkerDownDF,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/",
    "SubCluster",sck,"_MarkersDown.RDS"
  )
)

saveRDS(
  PathwayUpDF,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/",
    "SubCluster",sck,"_PathwaysUp.RDS"
  )
)

saveRDS(
  PathwayDownDF,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/",
    "SubCluster",sck,"_PathwaysDown.RDS"
  )
)

cat("\n\nDONE! >:(\n\n")

rm(list=ls())
gc()