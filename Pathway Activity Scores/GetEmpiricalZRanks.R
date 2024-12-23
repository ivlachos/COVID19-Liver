#!/usr/bin/env Rscript


rm(list=ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    opt_str = c("--CHUNK"),
    action="store",
    type="integer",
    default=27,
    help="Chunk index",
    metavar = "integer"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# ==== Functions =====
chunk2 = function (x, n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}

# original function to compute scaled and centered ranks
GetZRanks = function(ExprsMat){
  zrnk = mclapply(
    1:ncol(ExprsMat),
    function(ic){
      c(scale(rank(ExprsMat[,ic],ties.method = "min")))
    },
    mc.cores = 12
  )
  
  zrnk = do.call(cbind,zrnk)
  
  dimnames(zrnk) = dimnames(ExprsMat)
  
  return(zrnk)
}



# helper functions to compute expected value and variance of ranks
ExpectedMeanRank = function(n0,n1){
  ((n1+1)*(n1+2*n0))/(2*(n1+n0))
}


# CORRECTED VERSION
ExpectedVarianceRank = function(n0,n1){
  n = n0 + n1
  
  A = n1/n
  B = (n1**2 -1)/12
  C = (n1*n0)/(n**2)
  D = ((n1+2*n0-1)**2)/4
  
  A*B + C*D
}




# # Version Below is WRONG!!
# ExpectedVarianceRank = function(n0,n1){
#   n = n0 + n1
#   A = (n1*(n1**2 -1))/12
#   B = (n1*n0)/(n**2)
#   C = ((n1+2*n0-1)**2)/2
#   
#   A + B*C
# }


# scale and center ranks based on theoretical mean and variance
GetExpectedZRanks = function(ExprsMat,num_cores=12){
  zrnk = mclapply(
    1:ncol(ExprsMat),
    function(ic){
      # get rank
      ri = rank(ExprsMat[,ic],ties.method = "min")
      # get mean and variance
      n0 = sum(ExprsMat[,ic] == 0)
      n1 = nrow(ExprsMat) - n0
      
      mean_r = ExpectedMeanRank(n0,n1)
      var_r = ExpectedVarianceRank(n0,n1)
      
      c((ri - mean_r)/sqrt(var_r))
    },
    mc.cores = num_cores
  )
  
  zrnk = do.call(cbind,zrnk)
  
  dimnames(zrnk) = dimnames(ExprsMat)
  
  return(zrnk)
}

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


# function to estimate pathway scores
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



cat("\nReading data ...\n")
# ===== Liver Data =====
liver_logcounts <- readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2/data/SelectLiverNoDoubletsSeuratLogCounts.RDS"
)

cat("\nSplitting data ...\n")
ici = as.numeric( opt$CHUNK )
# ==== Split Genes ====
AllSamples = colnames(liver_logcounts)
# shuffle genes 
set.seed(664)
AllSamples = sample(AllSamples,length(AllSamples))

SampleChunk = chunk2(AllSamples,100)

cat(paste0("\n\nEstimating ZRanks for ",ici,"/100\n\n"))
# ==== Expression Ranks ======
zrnkMat = GetEmpiricalZRanks(
  ExprsMat=liver_logcounts[,SampleChunk[[ici]]],
  num_cores=4
)

cat(paste0("\n\nSaving ZRanks for ",ici,"/100\n\n"))
saveRDS(
  zrnkMat,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
    "EmpiricalZRanks_Part",ici,".RDS"
  )
)

cat("\n\nDONE! >:(\n\n")