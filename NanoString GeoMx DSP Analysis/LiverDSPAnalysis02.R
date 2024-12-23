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
liver_annot = readRDS(
  "/data/work/Projects/BrScRNAseq/NanoString/LiverAnnotation_CD45Extended.RDS"
)
rownames(liver_annot) = gsub(".","-",rownames(liver_annot),fixed=T)

table(liver_annot[,c("SlideAnnotation","Patient")])


liver_annot$log_aoi_size = log(liver_annot$aoi_size)

# ===== Select ROIs ====
liver_annot = liver_annot[grepl("Zone",liver_annot$SlideAnnotation),]
liver_annot$log_nuclei_counts = log(as.numeric(liver_annot$nuclei_counts))
liver_annot$nuclei_counts = SegmentProperties[rownames(liver_annot),"nuclei_counts"]

# ==== WTA Raw Counts ====
TargetsWTA = setdiff(TargetsWTA,"Neg Probe")
RawCountsWTA = BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsWTA,rownames(liver_annot)]
rownames(RawCountsWTA) = BioProbeCountMatrix$TargetName[BioProbeCountMatrix$TargetName %in% TargetsWTA]

# ==== COVID Raw Counts ====
TargetsCOVID = setdiff(TargetsCOVID,"SARS-CoV-2 Neg")
covid_genome = c("ORF1ab","ORF1ab_REV","S")

RawCountsCOVID = TargetCountMatrix[TargetCountMatrix$TargetName %in% TargetsCOVID,rownames(liver_annot)]
rownames(RawCountsCOVID) = TargetCountMatrix$TargetName[TargetCountMatrix$TargetName %in% TargetsCOVID]

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

# ==== Ranks ====
zrnkMat = GetEmpiricalZRanks(
  ExprsMat=rbind(RawCountsWTA,RawCountsCOVID),
  num_cores=8
)


# ==== Pathway Scores =====
gsLst = lapply(gsLst,function(x){x[x %in% rownames(zrnkMat)]})
gsLst = gsLst[sapply(gsLst,length) > 0]

go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])

zone_pathways = gsDF$Pathway[grep("Zone",gsDF$Category)]

obs_markerScores = mclapply(
  gsLst,
  function(x){
    GetPathwayScore(zrnkMat = zrnkMat, gs = x)
  },
  mc.cores=12
)
obs_markerScores = do.call(rbind,obs_markerScores)

obs_markerScores= obs_markerScores[,rownames(liver_annot)]

pathway_se = SummarizedExperiment(assays = obs_markerScores)
pathway_model = model.matrix( ~ SlideAnnotation + Patient + log_aoi_size + log_nuclei_counts,liver_annot)
pathway_nsv <- num.sv(obs_markerScores,pathway_model,seed=42)

pathway_iasvobj = iasva(pathway_se, pathway_model,verbose = FALSE, 
                        permute = FALSE, num.sv = pathway_nsv)

pathway_sv = pathway_iasvobj$sv
rownames(pathway_sv) = colnames(obs_markerScores)

pathwayScores_limma_adj <- removeBatchEffect(
  x = obs_markerScores,
  batch = liver_annot$Patient,
  covariates = cbind(pathway_sv,liver_annot[colnames(obs_markerScores),c("log_aoi_size","log_nuclei_counts")]),
  design = model.matrix(~ -1 + region,data = liver_annot)
)


# variance explained
pslimmaadj_varexp = scater::getVarianceExplained(
  x=pathwayScores_limma_adj,
  liver_annot[colnames(zrnkMat),c("aoi_size","log_aoi_size","Patient","nuclei_counts","log_nuclei_counts","SlideAnnotation")]
)

mypar(1,1)
boxplot(pslimmaadj_varexp,main="Pathway Score",ylim=c(0,100))


# ==== Zonation Gradient Setup  ====
zone_annot = readRDS("/data/work/Projects/BrScRNAseq/NanoString/NK/Liver/ZoneAnnotation.RDS")
zone_annot = zone_annot[rownames(liver_annot),]

# normalize distance to Zone 1 ROIs
z1d_max = tapply(zone_annot$Z1Dist,zone_annot$SpatialGroup,max)
zone_annot$Z1Dist2 = zone_annot$Z1Dist/z1d_max[zone_annot$SpatialGroup]


# ====> Zonation Gradient [Genes] <====
# ==== Zonation Gradient Distance (limma)  ====
zone_gene_design =  model.matrix(~ Z1Dist2,data = zone_annot)
zone_gene_lfit = lmFit( object = probe_limma_adj, design = zone_gene_design)
zone_gene_lfit = eBayes(zone_gene_lfit, trend = TRUE, robust = TRUE)
zone_gene_ttab = topTable(zone_gene_lfit,coef="Z1Dist2",number=Inf,adjust.method = "BH")



# ====> Zonation Gradient [Pathways] <====
# ==== Zonation Gradient Distance (limma)  ====
zone_pathway_design = model.matrix(~ Z1Dist2,data = zone_annot)
zone_pathway_lfit = lmFit( object = pathwayScores_limma_adj, design = zone_pathway_design)
zone_pathway_lfit = eBayes(zone_pathway_lfit)#, trend = TRUE, robust = TRUE)
zone_pathway_ttab = topTable(zone_pathway_lfit,coef="Z1Dist2",number=Inf,adjust.method = "BH")



zonecat_design = model.matrix(~ -1 + region,data = zone_annot)
colnames(zonecat_design) = gsub("region","",colnames(zonecat_design))
zonecat_gene_lfit = lmFit( object = probe_limma_adj, design = zonecat_design)
zonecat_pathway_lfit = lmFit( object = pathwayScores_limma_adj, design = zonecat_design)


# ====> Save Results <=====
# for shinyapp
# PCA embeddings
gVar = genefilter::rowVars(probe_limma_adj)
gVar = sort(gVar,decreasing=T)
zone_pcout = prcomp(x=t(probe_limma_adj[names(gVar)[1:1000],]),scale. = T, center = T)

ggplot(
  data.frame(
    zone_pcout$x,
    zone_annot[rownames(zone_pcout$x),]
  ),
  aes(x=PC1,y=PC2,color=region)
) +
  geom_point()

saveRDS(
  zone_pcout,
  "/data/shiny/shiny-server/apps/ZonationDSP/data/ZonesProbeLimmaAdjPCAtop1000HVGs.RDS"
)


# saveRDS(
#   zone_annot,
#   "/data/shiny/shiny-server/apps/ZonationDSP/data/ZoneAnnotation.RDS"
# )

saveRDS(
  probe_limma_adj,
  "/data/shiny/shiny-server/apps/ZonationDSP/data/ZoneProbeLimmaAdj.RDS"
)

saveRDS(
  pathwayScores_limma_adj,
  "/data/shiny/shiny-server/apps/ZonationDSP/data/ZoneProbePathwaySVsScores.RDS"
)

saveRDS(
  zone_pathway_lfit,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/ZonationProbelimmaSVsPathwayLFit.RDS"
)

saveRDS(
  zone_gene_lfit,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/LiverDSP/ZonationProbelimmaSVsLFit.RDS"
)

# ====> Upload Results <====
zonation_dribble = "https://drive.google.com/drive/u/0/folders/1SsOVY1iZZ-lvBRFtWm8UM-Jpl2BlzQ1V"

# ==== Genes Up ====
GenesUp = zone_gene_ttab %>% 
  dplyr::filter(
    logFC > 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Gene") %>%
  dplyr::select(-c(AveExpr,t,B))
GenesUp = cbind(GenesUp,zonecat_gene_lfit$coefficients[GenesUp$Gene,])
GenesUp = GenesUp[order(abs(GenesUp$logFC),decreasing = T),]
rownames(GenesUp) = NULL


write_tsv(
  x = GenesUp,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingGenes.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingGenes.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)

# ==== Genes Down ====
GenesDown = zone_gene_ttab %>% 
  dplyr::filter(
    logFC < 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Gene") %>%
  dplyr::select(-c(AveExpr,t,B))

GenesDown = GenesDown[order(abs(GenesDown$logFC),decreasing = T),]
GenesDown = cbind(GenesDown,zonecat_gene_lfit$coefficients[GenesDown$Gene,])
GenesDown = GenesDown[order(abs(GenesDown$logFC),decreasing = T),]
rownames(GenesDown) = NULL

write_tsv(
  x = GenesDown,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingGenes.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingGenes.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)



# ==== Pathways Up ====
PathwaysUp = zone_pathway_ttab %>% 
  dplyr::filter(
    logFC > 0,adj.P.Val < 0.1
  ) %>%
  tibble::rownames_to_column(var="Pathway") %>%
  dplyr::select(-c(AveExpr,t,B))

PathwaysUp = cbind(PathwaysUp,zonecat_pathway_lfit$coefficients[PathwaysUp$Pathway,])
PathwaysUp = PathwaysUp[order(abs(PathwaysUp$logFC),decreasing = T),]
rownames(PathwaysUp) = NULL


write_tsv(
  x = PathwaysUp,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingPathways.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingPathways.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)

write_tsv(
  x = PathwaysUp[PathwaysUp$Pathway %in% go_pathways,],
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingCuratedPathways.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsIncreasingCuratedPathways.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)




# ==== Pathways Down ====
PathwaysDown = zone_pathway_ttab %>% 
  dplyr::filter(
    logFC < 0,adj.P.Val < 0.1
  ) %>%
  tibble::rownames_to_column(var="Pathway") %>%
  dplyr::select(-c(AveExpr,t,B))

PathwaysDown = cbind(PathwaysDown,zonecat_pathway_lfit$coefficients[PathwaysDown$Pathway,])
PathwaysDown = PathwaysDown[order(abs(PathwaysDown$logFC),decreasing = T),]
rownames(PathwaysDown) = NULL


write_tsv(
  x = PathwaysDown,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingPathways.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingPathways.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)

write_tsv(
  x = PathwaysDown[PathwaysDown$Pathway %in% go_pathways,],
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingCuratedPathways.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/",
    "ZonationSVsDecreasingCuratedPathways.tsv"
  ),
  path=as_dribble(zonation_dribble),
  type = "spreadsheet",
  overwrite = TRUE
)





# ===== Interactive Tables =====
scaled_coeff = t(scale(t(
  zonecat_gene_lfit$coefficients
)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})


# Genes Up
GenesUpSummaryDF = zone_gene_ttab %>% 
  dplyr::filter(
    logFC > 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Gene") %>%
  dplyr::select(-c(AveExpr,t,B))
GenesUpSummaryDF = cbind(GenesUpSummaryDF,scaled_coeff[GenesUpSummaryDF$Gene,])
GenesUpSummaryDF = GenesUpSummaryDF[order(abs(GenesUpSummaryDF$logFC),decreasing = T),]
rownames(GenesUpSummaryDF) = NULL

font.size <- "7.5pt"
GenesUpDT = datatable(
  GenesUpSummaryDF,
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(GenesUpSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  GenesUpDT,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsIncreasingGenes.html"
  )
)

# Genes Down
GenesDownSummaryDF = zone_gene_ttab %>% 
  dplyr::filter(
    logFC < 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Gene") %>%
  dplyr::select(-c(AveExpr,t,B))
GenesDownSummaryDF = cbind(GenesDownSummaryDF,scaled_coeff[GenesDownSummaryDF$Gene,])
GenesDownSummaryDF = GenesDownSummaryDF[order(abs(GenesDownSummaryDF$logFC),decreasing = T),]
rownames(GenesDownSummaryDF) = NULL

font.size <- "7.5pt"
GenesDownDT = datatable(
  GenesDownSummaryDF,
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(GenesDownSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  GenesDownDT,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsDecreasingGenes.html"
  )
)


scaled_coeff = t(scale(t(
  zonecat_pathway_lfit$coefficients
)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})




# Pathways Up
PathwaysUpSummaryDF = zone_pathway_ttab %>% 
  dplyr::filter(
    logFC > 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Pathway") %>%
  dplyr::select(-c(AveExpr,t,B))
PathwaysUpSummaryDF = cbind(PathwaysUpSummaryDF,scaled_coeff[PathwaysUpSummaryDF$Pathway,])
PathwaysUpSummaryDF = PathwaysUpSummaryDF[order(abs(PathwaysUpSummaryDF$logFC),decreasing = T),]
rownames(PathwaysUpSummaryDF) = NULL

font.size <- "7.5pt"
PathwaysUpDT = datatable(
  PathwaysUpSummaryDF,
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(PathwaysUpSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  PathwaysUpDT,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsIncreasingPathways.html"
  )
)

font.size <- "7.5pt"
PathwaysUpDT2 = datatable(
  PathwaysUpSummaryDF[PathwaysUpSummaryDF$Pathway %in% go_pathways,],
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(PathwaysUpSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  PathwaysUpDT2,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsIncreasingCuratedPathways.html"
  )
)

# Pathways Down
PathwaysDownSummaryDF = zone_pathway_ttab %>% 
  dplyr::filter(
    logFC < 0,adj.P.Val < 0.05
  ) %>%
  tibble::rownames_to_column(var="Pathway") %>%
  dplyr::select(-c(AveExpr,t,B))
PathwaysDownSummaryDF = cbind(PathwaysDownSummaryDF,scaled_coeff[PathwaysDownSummaryDF$Pathway,])
PathwaysDownSummaryDF = PathwaysDownSummaryDF[order(abs(PathwaysDownSummaryDF$logFC),decreasing = T),]
rownames(PathwaysDownSummaryDF) = NULL

font.size <- "7.5pt"
PathwaysDownDT = datatable(
  PathwaysDownSummaryDF,
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(PathwaysDownSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  PathwaysDownDT,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsDecreasingPathways.html"
  )
)

font.size <- "7.5pt"
PathwaysDownDT2 = datatable(
  PathwaysDownSummaryDF[PathwaysDownSummaryDF$Pathway %in% go_pathways,],
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
  formatRound(columns = c("logFC","Zone1","Zone2","Zone3"),digits=2) %>%
  formatRound(columns=c("P.Value","adj.P.Val"),digits=4) %>% 
  formatStyle(which(colnames(PathwaysDownSummaryDF) %in% c("Zone1","Zone2","Zone3")), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))

DT::saveWidget(
  PathwaysDownDT2,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/LiverDSP/Zonation/DT/",
    "ZonationSVsDecreasingCuratedPathways.html"
  )
)
