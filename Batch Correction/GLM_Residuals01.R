
rm(list=ls())
options(stringsAsFactors = FALSE)

library(parallel)
library(edgeR)
library(SummarizedExperiment)


# ==== Functions  ====
GetPearsonResidualsGamPoi <- function(fit){
  # fit: a GLM fit from glmGamPoi 
  y = assay(as.list(fit)[["data"]])
  mu = as.list(fit)[["Mu"]]
  phi = as.list(fit)[["overdispersions"]]
  
  
  v <- mu*(1+phi*mu)
  res <- (y-mu) / sqrt(v)
  return(res)
}


GetDevianceResidualsGamPoi <- function(fit){
  # fit: a GLM fit from glmGamPoi 
  y = assay(as.list(fit)[["data"]])
  mu = as.list(fit)[["Mu"]]
  phi = as.list(fit)[["overdispersions"]]
  
  d <- nbinomUnitDeviance(y,mu,phi)
  res <- sign(y-mu) * sqrt(d)
  return(res)
}

GetMidPQuantileResidualsGamPoi <- function(fit){
  # fit: a GLM fit from glmGamPoi 
  y = assay(as.list(fit)[["data"]])
  mu = as.list(fit)[["Mu"]]
  phi = as.list(fit)[["overdispersions"]]
  
  
  res <- zscoreNBinom(y,mu,size=1/phi)
  return(res)
}



# ===== Pearson Residuals ====
pearsonLst = mclapply(
  1:10,
  function(ic){
    glmRes = readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/GLM/SelectLiverNoDoublets_glm_gpFit_",
        ic,".RDS"
      )
    )
    
    return(GetPearsonResidualsGamPoi(glmRes))
    
  },
  mc.cores = 12
)

PearsonResiduals = do.call(rbind,pearsonLst)

PearsonResiduals[is.na(PearsonResiduals)] = 0

saveRDS(
  PearsonResiduals,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoubletsGLMGamPoiPearsonResiduals.RDS"
)

rm(PearsonResiduals,pearsonLst)



# ===== Deviance Residuals ====
devianceLst = mclapply(
  1:10,
  function(ic){
    glmRes = readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/GLM/SelectLiverNoDoublets_glm_gpFit_",
        ic,".RDS"
      )
    )
    
    return(GetDevianceResidualsGamPoi(glmRes))
    
  },
  mc.cores = 12
)


DevianceResiduals = do.call(rbind,devianceLst)

DevianceResiduals[is.na(DevianceResiduals)] = 0

saveRDS(
  DevianceResiduals,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoubletsGLMGamPoiDevianceResiduals.RDS"
)

rm(DevianceResiduals,devianceLst)


# ===== MidP Quantile Residuals ====
MidPQLst = mclapply(
  1:10,
  function(ic){
    glmRes = readRDS(
      paste0(
        "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
        "SelectLiversNoDoublets/GLM/SelectLiverNoDoublets_glm_gpFit_",
        ic,".RDS"
      )
    )
    
    return(GetMidPQuantileResidualsGamPoi(glmRes))
    
  },
  mc.cores=8
)


MidPQResiduals = do.call(rbind,MidPQLst)



saveRDS(
  MidPQResiduals,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoubletsGLMGamPoiMidPQResiduals.RDS"
)

rm(MidPQResiduals,MidPQLst)