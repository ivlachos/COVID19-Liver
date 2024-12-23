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

# 
# cat("\n\nLoading pathways ... \n\n")
# # ==== Pathways ====
# tmpGenes = toupper(c(
#   'cdk1',
#   'pcna',
#   'e2f1',
#   'pold3',
#   'mki67',
#   'orc2',
#   'rrm1',
#   'cdt1',
#   'cdc6',
#   'mcm5',
#   'rrm2',
#   'orc1',
#   'ticrr'
# ))
# 
# 
# gsLst = readRDS(
#   "/data/work/Collabs/Covid/data/CuratedPathways.RDS"
# )
# 
# gsLst = lapply(gsLst,unique)
# 
# # Sanchez-Taltavull Signature
# gsLst[["KC Proliferation"]] = tmpGenes
# 
# gsDF = readRDS(
#   "/data/work/Collabs/Covid/data/CuratedPathwaysDF.RDS"
# )
# 
# rm(tmpGenes)

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


# empAdjZRnks = readRDS(
#   paste0(
#     "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
#     "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
#     "AdjustedEmpiricalZRanks.RDS"
#   )
# )
# 
# empAdjZRnks = empAdjZRnks[,rownames(liver_metadata)]

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


cat(paste0("\n\n [",sck,"] Estimating percent expressed ...\n\n"))
TmpPctExprs = mclapply(
  1:niter,
  function(ic){
    rnaplus_matches = MatchRNAPlus(sck=sck,ShowPlot=F,Verbose = F)
    rowMeans(liver_counts[,rnaplus_matches] > 0)*100
  },
  mc.cores = 4
)
TmpPctExprs = do.call(cbind,TmpPctExprs)

RNAMinusPctExprs = rowMeans(TmpPctExprs)


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
  RNAMinusPctExprs,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/COVIDPlusDE/",
    "SubCluster",sck,"/",
    "SubCluster",sck,"_PctEsxprs.RDS"
  )
)

cat("\n\nDONE! >:(\n\n")

rm(list=ls())
gc()