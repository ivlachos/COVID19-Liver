rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(edgeR)
library(iasva)
library(sva)
library(SummarizedExperiment)
library(parallel)
library(tidyverse)
library(rafalib)
library(pheatmap)
library(patchwork)
library(googledrive)
library(RColorBrewer)
library(DT)
Sys.setenv("RSTUDIO_PANDOC" = "/data/resources/tools/pandoc/")

# ==== Functions ====
# functions to get residuals from edgeR GLM fits
GetPearsonResiduals <- function(fit){
  # fit: a GLM fit from edgeR [DGEGLM]
  y <- fit$counts
  mu <- fit$fitted.values
  phi <- fit$dispersion
  v <- mu*(1+phi*mu)
  res <- (y-mu) / sqrt(v)
  return(res)
}

GetDevianceResiduals <- function(fit){
  # fit: a GLM fit from edgeR [DGEGLM]
  y <- fit$counts
  mu <- fit$fitted.values
  phi <- fit$dispersion
  
  d <- nbinomUnitDeviance(y,mu,phi)
  res <- sign(y-mu) * sqrt(d)
  return(res)
}

GetMidPQuantileResiduals <- function(fit){
  # fit: a GLM fit from edgeR [DGEGLM]
  y <- fit$counts
  mu <- fit$fitted.values
  phi <- fit$dispersion
  
  res <- zscoreNBinom(y,mu,size=1/phi)
  return(res)
}

BatchDesignMatrix = function(batch = NULL, batch2 = NULL, covariates = NULL, 
                             design = matrix(1, ncol(x), 1)) 
{
  if (is.null(batch) && is.null(batch2) && is.null(covariates)) 
    return(as.matrix(x))
  if (!is.null(batch)) {
    batch <- as.factor(batch)
    contrasts(batch) <- contr.sum(levels(batch))
    batch <- model.matrix(~batch)[, -1, drop = FALSE]
  }
  if (!is.null(batch2)) {
    batch2 <- as.factor(batch2)
    contrasts(batch2) <- contr.sum(levels(batch2))
    batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
  }
  if (!is.null(covariates)) 
    covariates <- as.matrix(covariates)
  X.batch <- cbind(batch, batch2, covariates)
  return(cbind(design, X.batch))
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


# ==== Liver DSP Data ====

# ==== Raw Data ====
# **SegmentProperties.txt - The current annotations that we have available for each segment within each ROI. 
# This will be updated once the BIDMC team provides complete annotations for the slides which will then be used for analysis
SegmentProperties = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_SegmentProperties.txt",
  data.table = FALSE
)

rownames(SegmentProperties) = SegmentProperties$Sample_ID
SegmentProperties$Patient = factor(SegmentProperties$Patient)
SegmentProperties$lobule = factor(SegmentProperties$lobule)


# **TargetProperties.txt - Annotation of the probes & pools from which they are derived - 1:1 with the rows of the BioProbeCountMatrix
TargetProperties = fread("/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_TargetProperties.txt")

# Count Matrices (targets/probes as rows, samples as columns; column 1 is target name):

# *BIDMC_Covid_BioProbeCountMatrix.txt - uncollapsed probe count matrix. 
# This is effectively the deduplicated count matrix after processing FASTQ files. 
# No QC has been performed on this beyond standard FASTQ processing. Pipeline information can be provided.
BioProbeCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_BioProbeCountMatrix.txt",
  data.table = FALSE
)


# *BIDMC_Covid_TargetCountMatrix.txt - collapsed target count matrix. 
# This is effectively the raw data file if you do not want to work at a per probe level 
# (only impacts negative controls and COVID spike-ins which have more than 1 probe per target). 
# Targets are reported as the geometric mean of the probes after QC'ing probes for outliers which 
# may add noise to the target expression in certain AOIs. More information can be provided about outlier calling methods.
TargetCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/M-317 BIDMC_COVID WTA/preliminary_data_20200714/BIDMC_Covid_TargetCountMatrix.txt",
  data.table = FALSE
)

# **NegNorm* - Count matrix normalized against the geometric mean of the negative control probes. 
NegNormCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_NegNorm_TargetCountMatrix.txt",
  data.table = FALSE
)

SampleIDs = grep("DSP",colnames(BioProbeCountMatrix),value=T)


# ==== WTA Negative Probes ====
NegProbeIDs = BioProbeCountMatrix$ProbeDisplayName[BioProbeCountMatrix$TargetName == "Neg Probe"]
TargetsWTA = unique(TargetProperties$TargetName[TargetProperties$Pooling == "1"])

NegProbeCounts = BioProbeCountMatrix[BioProbeCountMatrix$TargetName == "Neg Probe",SampleIDs] 
rownames(NegProbeCounts) = NegProbeIDs

TotalCountsWTA = colSums(BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsWTA,SampleIDs])

# ===== COVID Negative Probes =====
COVIDNegProbeIDs = BioProbeCountMatrix$ProbeDisplayName[BioProbeCountMatrix$TargetName == "SARS-CoV-2 Neg"]
TargetsCOVID = unique(TargetProperties$TargetName[TargetProperties$Pooling == "2"])

COVIDNegProbeCounts = BioProbeCountMatrix[BioProbeCountMatrix$TargetName == "SARS-CoV-2 Neg",SampleIDs] 
rownames(COVIDNegProbeCounts) = COVIDNegProbeIDs

TotalCountsCOVID = colSums(BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsCOVID,SampleIDs])

# ==== Liver annotation ===== 
zone_annot = readRDS(
  "/data/work/Projects/BrScRNAseq/NanoString/NK/Liver/ZoneAnnotation.RDS"
)
rownames(zone_annot) = gsub(".","-",rownames(zone_annot),fixed=T)

liver_annot = readRDS(
  "/data/work/Projects/BrScRNAseq/NanoString/LiverAnnotation_CD45Extended.RDS"
)
rownames(liver_annot) = gsub(".","-",rownames(liver_annot),fixed=T)

table(liver_annot[,c("SlideAnnotation","Patient")])


liver_annot$log_aoi_size = log(liver_annot$aoi_size)

# ===== Select ROIs ====
liver_annot = liver_annot[!grepl("CD45",liver_annot$SlideAnnotation),]
liver_annot$nuclei_counts = SegmentProperties[rownames(liver_annot),"nuclei_counts"]
liver_annot$log_nuclei_counts = log(as.numeric(liver_annot$nuclei_counts))



# ==== WTA Raw Counts ====
TargetsWTA = setdiff(TargetsWTA,"Neg Probe")
RawCountsWTA = BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsWTA,rownames(liver_annot)]
rownames(RawCountsWTA) = BioProbeCountMatrix$TargetName[BioProbeCountMatrix$TargetName %in% TargetsWTA]

# ==== COVID Raw Counts ====
TargetsCOVID = setdiff(TargetsCOVID,"SARS-CoV-2 Neg")
covid_genome = c("ORF1ab","ORF1ab_REV","S")

RawCountsCOVID = TargetCountMatrix[TargetCountMatrix$TargetName %in% TargetsCOVID,rownames(liver_annot)]
rownames(RawCountsCOVID) = TargetCountMatrix$TargetName[TargetCountMatrix$TargetName %in% TargetsCOVID]


# ==== Ranks ====
zrnkMatUnAdj = GetEmpiricalZRanks(
  ExprsMat=rbind(RawCountsWTA,RawCountsCOVID),
  num_cores=8
)


# ==== Pathway Scores =====
gsLst = lapply(gsLst,function(x){x[x %in% rownames(zrnkMatUnAdj)]})
gsLst = gsLst[sapply(gsLst,length) > 0]

go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])


unadj_markerScores = mclapply(
  gsLst,
  function(x){
    GetPathwayScore(zrnkMat = zrnkMatUnAdj, gs = x)
  },
  mc.cores=12
)
unadj_markerScores = do.call(rbind,unadj_markerScores)

# ===== Adjusted Pathway Scores (SVs) ====
unadj_markerScores = unadj_markerScores[,rownames(liver_annot)]
pathway_se = SummarizedExperiment(assays = unadj_markerScores)
pathway_model = model.matrix( ~ SlideAnnotation + Patient + log_aoi_size + log_nuclei_counts,liver_annot)
pathway_nsv <- num.sv(unadj_markerScores,pathway_model,seed=42)

pathway_iasvobj = iasva(pathway_se, pathway_model,verbose = FALSE, 
                        permute = FALSE, num.sv = pathway_nsv)

pathway_sv = pathway_iasvobj$sv
rownames(pathway_sv) = colnames(unadj_markerScores)

pathwayScores_limma_adj <- removeBatchEffect(
  x = unadj_markerScores,
  batch = liver_annot$Patient,
  covariates = cbind(pathway_sv,liver_annot[colnames(unadj_markerScores),c("log_aoi_size","log_nuclei_counts")]),
  design = model.matrix(~ -1 + region,data = liver_annot)
)




# ===== Pathway Linear Model ====
# All Pathways
pathway_design =  model.matrix(~ -1 + region,data = liver_annot)
pathway_lfit = lmFit( object = pathwayScores_limma_adj, design = pathway_design)
pathway_lfit = eBayes(pathway_lfit)


AllRegions = unique(liver_annot$region)
# ==== Region Up Pathways ====
regionPathUpLst = mclapply(
  AllRegions,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "RegionContMat = makeContrasts(\n",
          "region",sck,"vsAll = ",
          "region",sck," - (",
          paste(
            paste0(
              "region",setdiff(AllRegions,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllRegions) -1,
          ",\nlevels = colnames(pathway_lfit$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(pathway_lfit,RegionContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "RegionContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllRegions,sck),
              function(x){
                paste0(
                  "region",sck,"vs",x," = ",
                  "region",sck," - ","region",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(pathway_lfit$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(pathway_lfit,RegionContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_up_ind = apply( 
      subclust_contrast2_toptab[,grep("region",colnames(subclust_contrast2_toptab))] ,1,
      function(x){all(x > 0)}
    )
    
    subclust_up_pathways = rownames(subclust_contrast2_toptab)[subclust_up_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_up_pathways,1:2],
      subclust_contrast2_toptab[subclust_up_pathways,c("F","P.Value","adj.P.Val")]
    )
    
    
    tmptab = tmptab[order(tmptab[,1],decreasing=T),]
    
    tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
    
    return(tmptab)
    
  },
  mc.cores = 4
)
names(regionPathUpLst) = AllRegions

# ==== Region Down Pathways ====
regionPathDownLst = mclapply(
  AllRegions,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "RegionContMat = makeContrasts(\n",
          "region",sck,"vsAll = ",
          "region",sck," - (",
          paste(
            paste0(
              "region",setdiff(AllRegions,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllRegions) -1,
          ",\nlevels = colnames(pathway_lfit$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(pathway_lfit,RegionContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "RegionContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllRegions,sck),
              function(x){
                paste0(
                  "region",sck,"vs",x," = ",
                  "region",sck," - ","region",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(pathway_lfit$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(pathway_lfit,RegionContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_down_ind = apply( 
      subclust_contrast2_toptab[,grep("region",colnames(subclust_contrast2_toptab))] ,1,
      function(x){all(x < 0)}
    )
    
    subclust_down_pathways = rownames(subclust_contrast2_toptab)[subclust_down_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_down_pathways,1:2],
      subclust_contrast2_toptab[subclust_down_pathways,c("F","P.Value","adj.P.Val")]
    )
    
    
    tmptab = tmptab[order(tmptab[,1],decreasing=F),]
    
    tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
    
    return(tmptab)
    
  },
  mc.cores = 4
)
names(regionPathDownLst) = AllRegions

# ===== PCA Plot (Up/Down Pathways) =====
region_pathways = unlist(c(
  lapply(regionPathUpLst,function(x){x$Pathway}),
  lapply(regionPathDownLst,function(x){x$Pathway})
),use.names = FALSE)

region_pathways = unique(region_pathways)

pathways_pcout = prcomp(t(pathwayScores_limma_adj[region_pathways,]),scale. = T, center = T)
pathways_varexp = pathways_pcout$sdev**2
pathways_varexp = 100*pathways_varexp/sum(pathways_varexp)

regionDF = data.frame(
  pathways_pcout$x,
  region = liver_annot[rownames(pathways_pcout$x),"region"]
) %>% group_by(region) %>% summarise(across(everything(), list(median)))

pathways_ppca = ggplot(
  data.frame(
    pathways_pcout$x,
    liver_annot[rownames(pathways_pcout$x),]
  ),
  aes(x=PC1,y=PC2,color=region) 
) +
  geom_point()  +
  xlab(
    paste0(
      "PC1 (",round(pathways_varexp[1],2),"%)"
    )
  ) +
  ylab(
    paste0(
      "PC2 (",round(pathways_varexp[2],2),"%)"
    )
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

pathways_ppca = pathways_ppca + geom_text_repel(data=regionDF,aes(x=PC1_1,y=PC2_1,label=region),inherit.aes = FALSE)


(pathways_ppca + ggtitle("Pathways"))

sapply(regionPathUpLst,function(x){length(x$Pathway)})
sapply(regionPathDownLst,function(x){length(x$Pathway)})


# ==== Pathway Summaries ====
region_pathway_coeff = pathway_lfit$coefficients[,grep("region",colnames(pathway_lfit$coefficients))]
# - rank/quantile of mean
PathwayQuantExprs = lapply(
  AllRegions,
  function(sck){
    CurrentRegion = paste0("region",sck)
    OtherRegions = setdiff(colnames(region_pathway_coeff),CurrentRegion)
    
    tmpQ = mclapply(
      1:nrow(region_pathway_coeff),
      function(ic){
        mean(region_pathway_coeff[ic,OtherRegions] <= region_pathway_coeff[ic,CurrentRegion])
      },
      mc.cores = 8
    )
    
    tmpQ = unlist(tmpQ)
    names(tmpQ) = rownames(region_pathway_coeff)
    return(tmpQ)
  }
)
PathwayQuantExprs = do.call(cbind,PathwayQuantExprs)
colnames(PathwayQuantExprs) = AllRegions


PathwaySummaries = mclapply(
  AllRegions,
  function(x){
    res = cbind(
      region_pathway_coeff[,paste0("region",x)],
      PathwayQuantExprs[rownames(region_pathway_coeff),x]
    )
    colnames(res) = paste0(x,c(".Avg",".Qnt"))
    return(res)
  },
  mc.cores=8
)
PathwaySummaries = do.call(cbind,PathwaySummaries)



# ===== Pathway Helper =====
HelperPathway = function(mrkLst,ic=1){
  res = data.frame(
    Pathway = mrkLst[[ic]]$Pathway,
    Region = names(mrkLst)[ic],
    OverallDiff = mrkLst[[ic]][,2],
    Pval = mrkLst[[ic]][,"P.Value"],
    AdjPval = mrkLst[[ic]][,"adj.P.Val"],
    QuantExprs = PathwayQuantExprs[mrkLst[[ic]]$Pathway,names(mrkLst)[ic]]
  )
  
  res = res[order(res$QuantExprs, decreasing=T),]
  rownames(res) = NULL
  return(res)
}

# ==== Pathways Up ====
RegionPathwaysUp = mclapply(
  seq_along(regionPathUpLst),
  function(ic){
    HelperPathway(mrkLst=regionPathUpLst,ic=ic)
  },
  mc.cores = 4
)
names(RegionPathwaysUp) = names(regionPathUpLst)


# ==== Pathways Down ====
RegionPathwaysDown = mclapply(
  seq_along(regionPathDownLst),
  function(ic){
    HelperPathway(mrkLst=regionPathDownLst,ic=ic)
  },
  mc.cores = 4
)
names(RegionPathwaysDown) = names(regionPathDownLst)

RegionPathwaysDown = lapply(RegionPathwaysDown,function(x){x[order(x$QuantExprs, decreasing=FALSE),]})

# ==== Raw Probe Counts ====
probe_raw_dge = DGEList(
  counts = rbind(RawCountsWTA[,rownames(liver_annot)],RawCountsCOVID[,rownames(liver_annot)])
)
probe_raw_dge = calcNormFactors(probe_raw_dge)

# ==== Raw Probe SVs ====
probe_logcpm = edgeR::cpm(probe_raw_dge,log=T)
probe_se <- SummarizedExperiment(assays = probe_logcpm)
probe_model = model.matrix( ~ SlideAnnotation + Patient,liver_annot)
probe_nsv <- num.sv(probe_logcpm,probe_model,seed=42)

probe_iasvobj = iasva(probe_se, probe_model,verbose = FALSE, 
                    permute = FALSE, num.sv = probe_nsv)

probe_sv = probe_iasvobj$sv
rownames(probe_sv) = colnames(probe_logcpm)

# ==== Limma: Remove Batch Effects =====
probe_limma_adj <- removeBatchEffect(
  x = probe_logcpm,
  batch = liver_annot$Patient,
  covariates = probe_sv,
  design = model.matrix(~ -1 + region,data = liver_annot)
)


# ===== Gene Linear Model =====
region_design = model.matrix(~ -1 + region,data = liver_annot)
gene_lfit = lmFit( object = probe_limma_adj, design = region_design)
gene_lfit = eBayes(gene_lfit, trend = TRUE, robust = TRUE)



# ==== Region Up Markers ====
regionMrkUpLst = mclapply(
  AllRegions,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "RegionContMat = makeContrasts(\n",
          "region",sck,"vsAll = ",
          "region",sck," - (",
          paste(
            paste0(
              "region",setdiff(AllRegions,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllRegions) -1,
          ",\nlevels = colnames(gene_lfit$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(gene_lfit,RegionContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "RegionContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllRegions,sck),
              function(x){
                paste0(
                  "region",sck,"vs",x," = ",
                  "region",sck," - ","region",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(gene_lfit$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(gene_lfit,RegionContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_up_ind = apply( 
      subclust_contrast2_toptab[,grep("region",colnames(subclust_contrast2_toptab))] ,1,
      function(x){all(x > 0)}
    )
    
    subclust_up_markers = rownames(subclust_contrast2_toptab)[subclust_up_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_up_markers,1:2],
      subclust_contrast2_toptab[subclust_up_markers,c("F","P.Value","adj.P.Val")]
    )
    
    
    tmptab = tmptab[order(tmptab[,1],decreasing=T),]
    
    tmptab = tibble::rownames_to_column(tmptab,var="Gene")
    
    return(tmptab)
    
  },
  mc.cores = 4
)
names(regionMrkUpLst) = AllRegions

# ==== Region Down Markers ====
regionMrkDownLst = mclapply(
  AllRegions,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "RegionContMat = makeContrasts(\n",
          "region",sck,"vsAll = ",
          "region",sck," - (",
          paste(
            paste0(
              "region",setdiff(AllRegions,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllRegions) -1,
          ",\nlevels = colnames(gene_lfit$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(gene_lfit,RegionContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "RegionContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllRegions,sck),
              function(x){
                paste0(
                  "region",sck,"vs",x," = ",
                  "region",sck," - ","region",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(gene_lfit$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(gene_lfit,RegionContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_down_ind = apply( 
      subclust_contrast2_toptab[,grep("region",colnames(subclust_contrast2_toptab))] ,1,
      function(x){all(x < 0)}
    )
    
    subclust_down_markers = rownames(subclust_contrast2_toptab)[subclust_down_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_down_markers,1:2],
      subclust_contrast2_toptab[subclust_down_markers,c("F","P.Value","adj.P.Val")]
    )
    
    
    tmptab = tmptab[order(tmptab[,1],decreasing=F),]
    
    tmptab = tibble::rownames_to_column(tmptab,var="Gene")
    
    return(tmptab)
    
  },
  mc.cores = 4
)
names(regionMrkDownLst) = AllRegions


sapply(regionMrkUpLst,function(x){length(x$Gene)})
sapply(regionMrkDownLst,function(x){length(x$Gene)})
sapply(regionPathUpLst,function(x){length(x$Pathway)})
sapply(regionPathDownLst,function(x){length(x$Pathway)})


# ==== Expression Summaries ====
region_coeff = gene_lfit$coefficients[,grep("region",colnames(gene_lfit$coefficients))]
# - rank/quantile of mean
CoeffQuantExprs = lapply(
  AllRegions,
  function(sck){
    CurrentRegion = paste0("region",sck)
    OtherRegions = setdiff(colnames(region_coeff),CurrentRegion)
    
    tmpQ = mclapply(
      1:nrow(region_coeff),
      function(ic){
        mean(region_coeff[ic,OtherRegions] <= region_coeff[ic,CurrentRegion])
      },
      mc.cores = 8
    )
    
    tmpQ = unlist(tmpQ)
    names(tmpQ) = rownames(region_coeff)
    return(tmpQ)
  }
)
CoeffQuantExprs = do.call(cbind,CoeffQuantExprs)
colnames(CoeffQuantExprs) = AllRegions


RegionSummaries = mclapply(
  AllRegions,
  function(x){
    res = cbind(
      region_coeff[,paste0("region",x)],
      CoeffQuantExprs[rownames(region_coeff),x]
    )
    colnames(res) = paste0(x,c(".Avg",".Qnt"))
    return(res)
  },
  mc.cores=8
)
RegionSummaries = do.call(cbind,RegionSummaries)

# ===== Marker Helper =====
HelperMarker = function(mrkLst,ic=1){
  res = data.frame(
    Gene = mrkLst[[ic]]$Gene,
    Region = names(mrkLst)[ic],
    OverallDiff = mrkLst[[ic]][,2],
    Pval = mrkLst[[ic]][,"P.Value"],
    AdjPval = mrkLst[[ic]][,"adj.P.Val"],
    QuantExprs = CoeffQuantExprs[mrkLst[[ic]]$Gene,names(mrkLst)[ic]]
  )
  
  res = res[order(res$QuantExprs, decreasing=T),]
  rownames(res) = NULL
  return(res)
}

# ==== Markers Up ====
RegionMarkersUp = mclapply(
  seq_along(regionMrkUpLst),
  function(ic){
    HelperMarker(mrkLst=regionMrkUpLst,ic=ic)
  },
  mc.cores = 4
)
names(RegionMarkersUp) = names(regionMrkUpLst)


# ==== Markers Down ====
RegionMarkersDown = mclapply(
  seq_along(regionMrkDownLst),
  function(ic){
    HelperMarker(mrkLst=regionMrkDownLst,ic=ic)
  },
  mc.cores = 4
)
names(RegionMarkersDown) = names(regionMrkDownLst)

RegionMarkersDown = lapply(RegionMarkersDown,function(x){x[order(x$QuantExprs, decreasing=FALSE),]})


# ===== Upload Markers ====
dribblePaths = c(
  "Portal" = "https://drive.google.com/drive/u/0/folders/1jwMz_pR9z8J4a5SeTDVpP5zsV8UWCcaN",
  "Zone1" = "https://drive.google.com/drive/u/0/folders/19YtBwqqbRARzujjNINf0hf4o3aztyne_",
  "Zone2" = "https://drive.google.com/drive/u/0/folders/1oFb4oEXI-I0DxPQKzjS89yWjxDDWXvkz",
  "Zone3" = "https://drive.google.com/drive/u/0/folders/1kRSdLDAhHBBruON5bpucd2PlL8nrZurz"
)


# Pathways Up
for(ck in AllRegions){
  tmpDF = RegionPathwaysUp[[ck]]
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"PathwaysSVsUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"PathwaysSVsUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  
  write_tsv(
    x = tmpDF[tmpDF$Pathway %in% go_pathways,],
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"CuratedPathwaysSVsUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"CuratedPathwaysSVsUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

# Pathways Down
for(ck in AllRegions){
  tmpDF = RegionPathwaysDown[[ck]]
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"PathwaysSVsDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"PathwaysSVsDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  write_tsv(
    x = tmpDF[tmpDF$Pathway %in% go_pathways,],
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"CuratedPathwaysSVsDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"CuratedPathwaysSVsDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}


# Markers Up
for(ck in AllRegions){
  tmpDF = RegionMarkersUp[[ck]]
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"MarkersSVsUp.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"MarkersSVsUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

# Markers Down
for(ck in AllRegions){
  tmpDF = RegionMarkersDown[[ck]]
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"MarkersSVsDown.tsv"
    )
  )
  
  drive_upload(
    media = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/",
      ck,"MarkersSVsDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

# ===== Interactive Tables =====
scaled_coeff = t(scale(t(
  region_pathway_coeff
)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})


# Pathways Up
for(ck in AllRegions){
  tmpDF = RegionPathwaysUp[[ck]]
  rownames(tmpDF) = NULL
  
  PathwayUpSummaryDF = data.frame(
    tmpDF,
    scaled_coeff[tmpDF$Pathway,]
  )
  
  rownames(PathwayUpSummaryDF) = NULL
  colnames(PathwayUpSummaryDF) = gsub("region","",colnames(PathwayUpSummaryDF))
  colnames(PathwayUpSummaryDF) = gsub("Quant","Qnt",colnames(PathwayUpSummaryDF) )
  
  font.size <- "7.5pt"
  PathwayUpDT = datatable(
    PathwayUpSummaryDF,
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(PathwayUpSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(PathwayUpSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayUpDT,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"PathwaysSVsUp.html"
    )
    
  )
  
  # curated pathways
  font.size <- "7.5pt"
  PathwayUpDT2 = datatable(
    PathwayUpSummaryDF[PathwayUpSummaryDF$Pathway %in% go_pathways,],
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(PathwayUpSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(PathwayUpSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayUpDT2,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"CuratedPathwaysSVsUp.html"
    )
    
  )
}

# Pathways Down
for(ck in AllRegions){
  tmpDF = RegionPathwaysDown[[ck]]
  rownames(tmpDF) = NULL
  
  PathwayDownSummaryDF = data.frame(
    tmpDF,
    scaled_coeff[tmpDF$Pathway,]
  )
  
  rownames(PathwayDownSummaryDF) = NULL
  colnames(PathwayDownSummaryDF) = gsub("region","",colnames(PathwayDownSummaryDF))
  colnames(PathwayDownSummaryDF) = gsub("Quant","Qnt",colnames(PathwayDownSummaryDF) )
  
  font.size <- "7.5pt"
  PathwayDownDT = datatable(
    PathwayDownSummaryDF,
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(PathwayDownSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(PathwayDownSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayDownDT,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"PathwaysSVsDown.html"
    )
    
  )
  
  # curated pathways
  
  font.size <- "7.5pt"
  PathwayDownDT2 = datatable(
    PathwayDownSummaryDF[PathwayDownSummaryDF$Pathway %in% go_pathways,],
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(PathwayDownSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(PathwayDownSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayDownDT2,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"CuratedPathwaysSVsDown.html"
    )
    
  )
}


scaled_coeff = t(scale(t(region_coeff)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})


# Markers Up
for(ck in AllRegions){
  tmpDF = RegionMarkersUp[[ck]]
  rownames(tmpDF) = NULL
  
  MarkerUpSummaryDF = data.frame(
    tmpDF,
    scaled_coeff[tmpDF$Gene,]
  )
  
  rownames(MarkerUpSummaryDF) = NULL
  colnames(MarkerUpSummaryDF) = gsub("region","",colnames(MarkerUpSummaryDF))
  colnames(MarkerUpSummaryDF) = gsub("Quant","Qnt",colnames(MarkerUpSummaryDF) )
  
  font.size <- "7.5pt"
  MarkerUpDT = datatable(
    MarkerUpSummaryDF,
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(MarkerUpSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(MarkerUpSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerUpDT,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"MarkersSVsUp.html"
    )
    
  )
}

# Markers Down
for(ck in AllRegions){
  tmpDF = RegionMarkersDown[[ck]]
  rownames(tmpDF) = NULL
  
  MarkerDownSummaryDF = data.frame(
    tmpDF,
    scaled_coeff[tmpDF$Gene,]
  )
  
  rownames(MarkerDownSummaryDF) = NULL
  colnames(MarkerDownSummaryDF) = gsub("region","",colnames(MarkerDownSummaryDF))
  colnames(MarkerDownSummaryDF) = gsub("Quant","Qnt",colnames(MarkerDownSummaryDF) )
  
  font.size <- "7.5pt"
  MarkerDownDT = datatable(
    MarkerDownSummaryDF,
    extensions = 'FixedColumns',
    #width=1800,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 2),
      initComplete = htmlwidgets::JS(
        "function(settings, json) {",
        paste0("$(this.api().table().container()).css({'font-size': '", font.size, "'});"),
        "}")
    ),
    rownames=FALSE
  ) %>%
    formatRound(columns = c(3,6:ncol(MarkerDownSummaryDF)),digits=2) %>%
    formatRound(columns=c("Pval","AdjPval"),digits=4) %>% 
    formatStyle(which(colnames(MarkerDownSummaryDF) %in% AllRegions), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerDownDT,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/LiverDSP/MarkersPathways/DT/",
      ck,"MarkersSVsDown.html"
    )
    
  )
}



# ==== Save Results ====
saveRDS(
  probe_limma_adj,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbelimmaAdj.RDS"
)

saveRDS(
  pathwayScores_limma_adj,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbePathwaySVsScores.RDS"
)

saveRDS(
  pathway_lfit,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbePathwaySVslimmaLFit.RDS"
)


saveRDS(
  regionPathUpLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbePathwaysSVsUpLst.RDS"
)

saveRDS(
  regionPathDownLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbePathwaysSVsDownLst.RDS"
)

saveRDS(
  regionMrkUpLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbeMarkersSVsUpLst.RDS"
)

saveRDS(
  regionMrkDownLst,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbeMarkersSVsDownLst.RDS"
)

saveRDS(
  gene_lfit,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/RegionProbelimmaSVsLFit.RDS"
)

