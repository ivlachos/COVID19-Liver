
rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(limma)
library(parallel)
library(stackoverflow)


# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_FullMetadata.RDS"
)

liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)


# ===== Scaled and Centered Ranks (Empirical) =====
empZRnks = lapply(
  1:100,
  function(ic){
    readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
        "EmpiricalZRanks_Part",ic,".RDS"
      )
    )
    
  }
)


empZRnks = do.call(cbind,empZRnks)
empZRnks = empZRnks[,rownames(liver_metadata)]



# ===== Adjusted Ranks ====
HelperFun =  function(ki){
  
  # check that at least 5 observations in the smallest group, otherwise don't split
  zero_ind = empZRnks[ki,] <= 0
  if(min(sum(zero_ind),sum(!zero_ind)) >= 5){
    
    lfitZero = lm(empZRnks[ki,zero_ind] ~ liver_metadata$nGene[zero_ind] + liver_metadata$lognUMI[zero_ind])
    lfitPos = lm(empZRnks[ki,!zero_ind] ~ liver_metadata$nGene[!zero_ind] + liver_metadata$lognUMI[!zero_ind])
    
    
    ZeroXB = as.matrix(lfitZero$model[,-1])%*%lfitZero$coefficients[-1]
    ZeroAdj = empZRnks[ki,zero_ind] - ZeroXB
    PosXB = as.matrix(lfitPos$model[,-1])%*%lfitPos$coefficients[-1]
    PosAdj = empZRnks[ki,!zero_ind] - PosXB
    
    
    tmpAdj = c(ZeroAdj[,,drop=T],PosAdj[,,drop=T])
    tmpAdj = tmpAdj[rownames(liver_metadata)]
    
    # PosAdj2 = scale(tmpAdj[!zero_ind]) + mean(empZRnks[ki,!zero_ind])
    # tmpAdj2 = c(
    #   PosAdj2[,,drop=T],
    #   tmpAdj[zero_ind]
    # )
    # tmpAdj2 = tmpAdj2[rownames(liver_metadata)
    
  }else{
    cat("Not enough samples! \nPooling samples. \n")
    lfitAll = lm(empZRnks[ki,] ~ liver_metadata$nGene + liver_metadata$lognUMI)
    
    
    
    AllXB = as.matrix(lfitAll$model[,-1])%*%lfitAll$coefficients[-1]
    AllAdj = empZRnks[ki,] - AllXB
    
    
    tmpAdj = AllAdj[,,drop=T]
    tmpAdj = tmpAdj[rownames(liver_metadata)]
    # tmpAdj2 = tmpAdj
  }
  
  cat(paste0(ki,"/",nrow(empZRnks),"\n"))
  return(tmpAdj)
}




empAdjZRnks = mclapply(
  1:nrow(empZRnks),
  HelperFun,
  mc.cores = 4
)
names(empAdjZRnks) = rownames(empZRnks)

empAdjZRnks = do.call(rbind,empAdjZRnks)

# save
saveRDS(
  empAdjZRnks,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Pathways/ExpressionRanksEmpirical/",
    "AdjustedEmpiricalZRanks.RDS"
  )
)
