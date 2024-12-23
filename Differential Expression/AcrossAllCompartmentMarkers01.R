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


gsLst = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathways.RDS"
)

gsLst = lapply(gsLst,unique)

# Sanchez-Taltavull Signature
gsLst[["KC Proliferation"]] = tmpGenes

gsDF = readRDS(
  "/data/work/Collabs/Covid/data/CuratedPathwaysDF.RDS"
)



# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1.RDS"
)
liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2a/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)

AllSubClusters = sort(unique(liver_metadata$SubCluster))

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

go_pathways = gsDF$Pathway[ grep("GO",gsDF$Category)]
go_pathways = unique(go_pathways[go_pathways %in% names(gsLst)])

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
  "/data/shiny/shiny-server/apps/SelectLiverNoDblPathways/data/AdjustedAllPathwayScoresAdjEmpZRanks.RDS"
)



# ==== Pathway Linear Model =====
compartment_design2 <- model.matrix(~ -1 + SubCluster, liver_metadata )
comp_path_lfit2 = lmFit(
  object = comp_pathwayScoresAdj,#[go_pathways,],
  design = compartment_design2
)
comp_path_lfit2 = eBayes(comp_path_lfit2)
rm(compartment_design2)

# ===== Gene Linear Model  =====
# filter genes
subpct = mclapply(
  AllSubClusters,
  function(sck){
    Matrix::rowMeans(liver_counts[,liver_metadata$SubCluster == sck] > 0)
  },
  mc.cores = 12
)
subpct = do.call(cbind,subpct)
colnames(subpct) = AllSubClusters
# keep genes where at least one subcluster has 5% cells expressing the gene
keep_ind = apply(subpct,1,function(x){any(x >= 0.05)})
rm(subpct)

# TMM Normalization
compartment_dge <- DGEList(
  counts=liver_counts[keep_ind,rownames(liver_metadata)], 
  group=liver_metadata$SubCluster
)
compartment_dge <- calcNormFactors(object=compartment_dge,method="TMM")

compartment_design2 <- model.matrix(~ -1 + SubCluster + Batch, liver_metadata )
colnames(compartment_design2) = make.names(colnames(compartment_design2) )

y <- new("EList")
y$E <- edgeR::cpm(compartment_dge, log = TRUE, prior.count = 3)
compartment_fit2 <- lmFit(y, design = compartment_design2)
compartment_fit2 <- eBayes(compartment_fit2, trend = TRUE, robust = TRUE)

rm(compartment_design2,y,compartment_dge)



# ==== Expression Summaries ====
# - percent expressed
CountsPctExprs = mclapply(
  AllSubClusters,
  function(sck){
    CurrentCells = rownames(liver_metadata)[liver_metadata$SubCluster == sck]
    rowMeans(liver_counts[,CurrentCells] > 0)*100
  },
  mc.cores = 8
)
CountsPctExprs = do.call(cbind,CountsPctExprs)
colnames(CountsPctExprs) = AllSubClusters

subclust_coeff = compartment_fit2$coefficients[,grep("SubCluster",colnames(compartment_fit2$coefficients))]
# - rank/quantile of mean
CoeffQuantExprs = lapply(
  AllSubClusters,
  function(sck){
    CurrentSubCluster = paste0("SubCluster",sck)
    OtherSubClusters = setdiff(colnames(subclust_coeff),CurrentSubCluster)
    
    tmpQ = mclapply(
      1:nrow(subclust_coeff),
      function(ic){
        mean(subclust_coeff[ic,OtherSubClusters] <= subclust_coeff[ic,CurrentSubCluster])
      },
      mc.cores = 8
    )
    
    tmpQ = unlist(tmpQ)
    names(tmpQ) = rownames(subclust_coeff)
    return(tmpQ)
  }
)
CoeffQuantExprs = do.call(cbind,CoeffQuantExprs)
colnames(CoeffQuantExprs) = AllSubClusters
rm(subclust_coeff)


subclust_coeff = compartment_fit2$coefficients[,grep("SubCluster",colnames(compartment_fit2$coefficients))]
colnames(subclust_coeff) = gsub("SubCluster","",colnames(subclust_coeff))



SubClusterSummaries = mclapply(
  AllSubClusters,
  function(x){
    res = cbind(
      subclust_coeff[,x],
      CountsPctExprs[rownames(subclust_coeff),x],
      CoeffQuantExprs[rownames(subclust_coeff),x]
    )
    colnames(res) = paste0(x,c(".Avg",".Pct",".Qnt"))
    return(res)
  },
  mc.cores=12
)
SubClusterSummaries = do.call(cbind,SubClusterSummaries)

saveRDS(
  SubClusterSummaries,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1_SubClusterLimmaSummaries.RDS"
)

# - percent expressed by compartment
CountsPctExprsCompartment = mclapply(
  as.character(0:4),
  function(ck){
    CurrentCells = rownames(liver_metadata)[liver_metadata$seurat_clusters == ck]
    rowMeans(liver_counts[,CurrentCells] > 0)*100
  },
  mc.cores = 8
)
CountsPctExprsCompartment = do.call(cbind,CountsPctExprsCompartment)
colnames(CountsPctExprsCompartment) = as.character(0:4)



pathway_coeff = comp_path_lfit2$coefficients[,grep("SubCluster",colnames(comp_path_lfit2$coefficients))]
# ===== Pathway Summaries =====
PathwayQuantExprs = lapply(
  AllSubClusters,
  function(sck){
    CurrentSubCluster = paste0("SubCluster",sck)
    OtherSubClusters = setdiff(colnames(pathway_coeff),CurrentSubCluster)
    
    tmpQ = mclapply(
      1:nrow(pathway_coeff),
      function(ic){
        mean(pathway_coeff[ic,OtherSubClusters] <= pathway_coeff[ic,CurrentSubCluster])
      },
      mc.cores = 8
    )
    
    tmpQ = unlist(tmpQ)
    names(tmpQ) = rownames(pathway_coeff)
    return(tmpQ)
  }
)
PathwayQuantExprs = do.call(cbind,PathwayQuantExprs)
colnames(PathwayQuantExprs) = AllSubClusters
rm(pathway_coeff)


PathwayPctExprs = mclapply(
  gsLst,
  function(x){
    colMeans(CountsPctExprs[x,,drop=F])
  },
  mc.cores = 12
)
PathwayPctExprs = do.call(rbind,PathwayPctExprs)

PathwayPctExprsCompartment = mclapply(
  gsLst,
  function(x){
    colMeans(CountsPctExprsCompartment[x,,drop=F])
  },
  mc.cores = 12
)
PathwayPctExprsCompartment = do.call(rbind,PathwayPctExprsCompartment)


pathway_coeff = comp_path_lfit2$coefficients[,grep("SubCluster",colnames(comp_path_lfit2$coefficients))]
colnames(pathway_coeff) = gsub("SubCluster","",colnames(pathway_coeff))

SubClusterPathwaySummaries = mclapply(
  AllSubClusters,
  function(x){
    res = cbind(
      pathway_coeff[,x],
      PathwayPctExprs[rownames(pathway_coeff),x],
      PathwayQuantExprs[rownames(pathway_coeff),x]
    )
    colnames(res) = paste0(x,c(".Score",".Pct",".Qnt"))
    return(res)
  },
  mc.cores=12
)
SubClusterPathwaySummaries = do.call(cbind,SubClusterPathwaySummaries)

saveRDS(
  SubClusterPathwaySummaries,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1_PathwayLimmaSummaries.RDS"
)

PathwayAllPctExprs = cbind(
  PathwayPctExprsCompartment,
  PathwayPctExprs[rownames(PathwayPctExprsCompartment),]
)

saveRDS(
  PathwayAllPctExprs,
  "/data/shiny/shiny-server/apps/SelectLiverNoDblPathways/data/AdjustedAllPathwayPctExprsAdjEmpZRanks.RDS"
)



# ==== SubCluster Up Markers ====
resubclustMrkUpLst = mclapply(
  AllSubClusters,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "SubClustContMat = makeContrasts(\n",
          "SubCluster",sck,"vsAll = ",
          "SubCluster",sck," - (",
          paste(
            paste0(
              "SubCluster",setdiff(AllSubClusters,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllSubClusters) -1,
          ",\nlevels = colnames(compartment_fit2$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(compartment_fit2,SubClustContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "SubClustContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllSubClusters,sck),
              function(x){
                paste0(
                  "SubCluster",sck,"vs",x," = ",
                  "SubCluster",sck," - ","SubCluster",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(compartment_fit2$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(compartment_fit2,SubClustContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_up_ind = apply( 
      subclust_contrast2_toptab[,grep("SubCluster",colnames(subclust_contrast2_toptab))] ,1,
      function(x){mean(x > 0) >= 0.75}
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
names(resubclustMrkUpLst) = AllSubClusters

# ==== SubCluster Down Markers ====
resubclustMrkDownLst = mclapply(
  AllSubClusters,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "SubClustContMat = makeContrasts(\n",
          "SubCluster",sck,"vsAll = ",
          "SubCluster",sck," - (",
          paste(
            paste0(
              "SubCluster",setdiff(AllSubClusters,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllSubClusters) -1,
          ",\nlevels = colnames(compartment_fit2$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(compartment_fit2,SubClustContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "SubClustContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllSubClusters,sck),
              function(x){
                paste0(
                  "SubCluster",sck,"vs",x," = ",
                  "SubCluster",sck," - ","SubCluster",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(compartment_fit2$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(compartment_fit2,SubClustContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_down_ind = apply( 
      subclust_contrast2_toptab[,grep("SubCluster",colnames(subclust_contrast2_toptab))] ,1,
      function(x){mean(x < 0) >= 0.75}
    )
    
    subclust_down_markers = rownames(subclust_contrast2_toptab)[subclust_down_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_down_markers,1:2],
      subclust_contrast2_toptab[subclust_down_markers,c("F","P.Value","adj.P.Val")]
    )
    
    
    tmptab = tmptab[order(tmptab[,1],decreasing=T),]
    tmptab = tibble::rownames_to_column(tmptab,var="Gene")
    
    return(tmptab)
    
  },
  mc.cores = 4
)
names(resubclustMrkDownLst) = AllSubClusters

# ==== SubCluster Up Pathways ====
resubclustPathUpLst = mclapply(
  AllSubClusters,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "SubClustContMat = makeContrasts(\n",
          "SubCluster",sck,"vsAll = ",
          "SubCluster",sck," - (",
          paste(
            paste0(
              "SubCluster",setdiff(AllSubClusters,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllSubClusters) -1,
          ",\nlevels = colnames(comp_path_lfit2$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(comp_path_lfit2,SubClustContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "SubClustContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllSubClusters,sck),
              function(x){
                paste0(
                  "SubCluster",sck,"vs",x," = ",
                  "SubCluster",sck," - ","SubCluster",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(comp_path_lfit2$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(comp_path_lfit2,SubClustContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_up_ind = apply( 
      subclust_contrast2_toptab[,grep("SubCluster",colnames(subclust_contrast2_toptab))] ,1,
      function(x){mean(x > 0) >= 0.75}
    )
    
    subclust_up_markers = rownames(subclust_contrast2_toptab)[subclust_up_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_up_markers,1:2],
      subclust_contrast2_toptab[subclust_up_markers,c("F","P.Value","adj.P.Val")]
    )
    
    tmptab = tmptab[order(tmptab[,1],decreasing=T),]
    tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
    
    return(tmptab)
  },
  mc.cores = 4
)
names(resubclustPathUpLst) = AllSubClusters


# ==== SubCluster Down Pathways ====
resubclustPathDownLst = mclapply(
  AllSubClusters,
  function(sck){
    # use contrast to compare to all other clusters
    eval(
      parse(
        text = paste0(
          "SubClustContMat = makeContrasts(\n",
          "SubCluster",sck,"vsAll = ",
          "SubCluster",sck," - (",
          paste(
            paste0(
              "SubCluster",setdiff(AllSubClusters,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllSubClusters) -1,
          ",\nlevels = colnames(comp_path_lfit2$coefficients)\n)"
        )
      )
    )
    
    # fit contrasts
    subclust_lfit <- contrasts.fit(comp_path_lfit2,SubClustContMat)
    subclust_lfit <- eBayes(subclust_lfit)
    
    subclust_contrast_toptab <- topTableF(subclust_lfit,number=Inf,adjust.method = "BH")
    
    
    eval(
      parse(
        text = paste0(
          "SubClustContMat2 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllSubClusters,sck),
              function(x){
                paste0(
                  "SubCluster",sck,"vs",x," = ",
                  "SubCluster",sck," - ","SubCluster",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(comp_path_lfit2$coefficients)\n)"
        )
      )
    )
    
    
    subclust_lfit2 <- contrasts.fit(comp_path_lfit2,SubClustContMat2)
    subclust_lfit2 <- eBayes(subclust_lfit2)
    
    subclust_contrast2_toptab <- topTableF(subclust_lfit2,number=Inf,adjust.method = "BH")
    
    
    subclust_down_ind = apply( 
      subclust_contrast2_toptab[,grep("SubCluster",colnames(subclust_contrast2_toptab))] ,1,
      function(x){mean(x < 0) >= 0.75}
    )
    
    subclust_down_markers = rownames(subclust_contrast2_toptab)[subclust_down_ind & subclust_contrast2_toptab$adj.P.Val < 0.1]
    
    tmptab = cbind(
      subclust_contrast_toptab[subclust_down_markers,1:2],
      subclust_contrast2_toptab[subclust_down_markers,c("F","P.Value","adj.P.Val")]
    )
    
    tmptab = tmptab[order(tmptab[,1],decreasing=T),]
    tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
    
    return(tmptab)
  },
  mc.cores = 4
)
names(resubclustPathDownLst) = AllSubClusters



cbind(
  sapply(resubclustMrkUpLst,nrow),
  sapply(resubclustMrkDownLst,nrow),
  sapply(resubclustPathUpLst,nrow),
  sapply(resubclustPathDownLst,nrow)
)


# ===== Marker Helper =====
HelperMarker = function(mrkLst,ic=1){
  res = data.frame(
    Gene = mrkLst[[ic]]$Gene,
    Cluster = names(mrkLst)[ic],
    OverallDiff = mrkLst[[ic]][,2],
    Pval = mrkLst[[ic]][,"P.Value"],
    AdjPval = mrkLst[[ic]][,"adj.P.Val"],
    QuantExprs = CoeffQuantExprs[mrkLst[[ic]]$Gene,names(mrkLst)[ic]],
    PctExprs = CountsPctExprs[mrkLst[[ic]]$Gene,names(mrkLst)[ic]]
  )
  res$QuantXPct = res$QuantExprs*res$PctExprs
  
  res = res[order(res$QuantXPct, decreasing=T),]
  rownames(res) = NULL
  return(res)
}

# ==== Markers Up ====
ClusterMarkersUp = mclapply(
  seq_along(resubclustMrkUpLst),
  function(ic){
    HelperMarker(mrkLst=resubclustMrkUpLst,ic=ic)
  },
  mc.cores = 4
)
names(ClusterMarkersUp) = names(resubclustMrkUpLst)

saveRDS(
  ClusterMarkersUp,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/",
    "OverallMarkersUp.RDS"
  )
)

# ==== Markers Down ====
ClusterMarkersDown = mclapply(
  seq_along(resubclustMrkDownLst),
  function(ic){
    HelperMarker(mrkLst=resubclustMrkDownLst,ic=ic)
  },
  mc.cores = 4
)
names(ClusterMarkersDown) = names(resubclustMrkDownLst)

saveRDS(
  ClusterMarkersDown,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/",
    "OverallMarkersDown.RDS"
  )
)


# ===== Pathway Helper =====
HelperPathway = function(mrkLst,ic=1){
  res = data.frame(
    Pathway = mrkLst[[ic]]$Pathway,
    Cluster = names(mrkLst)[ic],
    OverallDiff = mrkLst[[ic]][,2],
    Pval = mrkLst[[ic]][,"P.Value"],
    AdjPval = mrkLst[[ic]][,"adj.P.Val"],
    QuantScore = PathwayQuantExprs[mrkLst[[ic]]$Pathway,names(mrkLst)[ic]],
    PctExprs = PathwayPctExprs[mrkLst[[ic]]$Pathway,names(mrkLst)[ic]]
  )
  res$QuantXPct = res$QuantScore*res$PctExprs
  
  res = res[order(res$QuantXPct, decreasing=T),]
  rownames(res) = NULL
  return(res)
}

# ==== Pathways Up ====
ClusterPathwaysUp = mclapply(
  seq_along(resubclustPathUpLst),
  function(ic){
    HelperPathway(mrkLst=resubclustPathUpLst,ic=ic)
  },
  mc.cores = 4
)
names(ClusterPathwaysUp) = names(resubclustPathUpLst)

saveRDS(
  ClusterPathwaysUp,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/",
    "OverallPathwaysUp.RDS"
  )
)


# ==== Pathways Down ====
ClusterPathwaysDown = mclapply(
  seq_along(resubclustPathDownLst),
  function(ic){
    HelperPathway(mrkLst=resubclustPathDownLst,ic=ic)
  },
  mc.cores = 4
)
names(ClusterPathwaysDown) = names(resubclustPathDownLst)

saveRDS(
  ClusterPathwaysDown,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/",
    "OverallPathwaysDown.RDS"
  )
)

# ==== Upload Results ====
# group clusters by compartment
clustSplit = split(
  x=AllSubClusters,
  f=sapply(strsplit(
    x=AllSubClusters,
    split="_"
  ),head,1)
)

dribblePaths = c(
  "0" = "https://drive.google.com/drive/u/0/folders/1uoRFeU7ifs-TJ2aJjEcJPbdOgAI90Roi",
  "1" = "https://drive.google.com/drive/u/0/folders/1bGynYyxcQSrD19PdBVWsTh8PLOCaZoNo",
  "2" = "https://drive.google.com/drive/u/0/folders/162RCBnT1JZ-94kZcy6y6WrWdXLvaHP5h",
  "3" = "https://drive.google.com/drive/u/0/folders/1rPvcSTsM6JlrjAPlHB8lR_Uuu1Qa5p2a",
  "4" = "https://drive.google.com/drive/u/0/folders/1S1kyFrp7ZnhPGfzSRQYZq4bNRJ5bLLR_"
)

# Markers Up
for(ck in as.character(0:4)){
  tmpLst = ClusterMarkersUp[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$PctExprs >= 15,]})
  tmpDF = do.call(rbind,tmpLst)
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallMarkersUp.tsv"
    )
  )
  
  drive_upload(
    media =paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallMarkersUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

# Markers Down
summary(sapply(ClusterMarkersUp,function(x){sum(x$PctExprs >= 15)}))
genes_atleast15 = rownames(CountsPctExprs)[apply(CountsPctExprs,1,function(x){any( x >= 15)}) ]
summary(sapply(ClusterMarkersDown,function(x){sum(x$Gene %in% genes_atleast15)}))

for(ck in as.character(0:4)){
  tmpLst = ClusterMarkersDown[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$Gene  %in% genes_atleast15,]})
  tmpLst = tmpLst[sapply(tmpLst,nrow) > 0]
  
  tmpLst = lapply(tmpLst,function(x){x[order(x$OverallDiff,decreasing=F),]})
  
  tmpDF = do.call(rbind,tmpLst)
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallMarkersDown.tsv"
    )
  )
  
  drive_upload(
    media =paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallMarkersDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}


# Pathways Up
for(ck in as.character(0:4)){
  tmpLst = ClusterPathwaysUp[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$PctExprs >= 10,]})
  tmpDF = do.call(rbind,tmpLst)
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallPathwaysUp.tsv"
    )
  )
  
  drive_upload(
    media =paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallPathwaysUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}


# Pathways Down
summary(sapply(ClusterPathwaysUp,function(x){sum(x$PctExprs >= 10)}))
pathways_atleast10 = rownames(PathwayPctExprs)[apply(PathwayPctExprs,1,function(x){any( x >= 10)}) ]
summary(sapply(ClusterPathwaysDown,function(x){sum(x$Pathway %in% pathways_atleast10)}))


for(ck in as.character(0:4)){
  tmpLst = ClusterPathwaysDown[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$Pathway %in% pathways_atleast10,]})
  tmpLst = tmpLst[sapply(tmpLst,nrow) > 0]
  tmpLst = lapply(tmpLst,function(x){x[order(x$OverallDiff,decreasing=F),]})
  
  tmpDF = do.call(rbind,tmpLst)
  rownames(tmpDF) = NULL
  
  write_tsv(
    x = tmpDF,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallPathwaysDown.tsv"
    )
  )
  
  drive_upload(
    media =paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/Compartment",ck,"/",
      "Compartment",ck,"OverallPathwaysDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

subclust_coeff = compartment_fit2$coefficients[,grep("SubCluster",colnames(compartment_fit2$coefficients))]
scaled_coeff = t(scale(t(subclust_coeff)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})

# ==== Interactive Table (Markers Up) ====
for(ck in as.character(0:4)){
  tmpLst = ClusterMarkersUp[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$PctExprs >= 15,]})
  ClusterMarkersUpDF = do.call(rbind,tmpLst)
  rownames(ClusterMarkersUpDF) = NULL
  
  MarkerUpSummaryDF = data.frame(
    ClusterMarkersUpDF,
    scaled_coeff[ClusterMarkersUpDF$Gene,]
  )
  rownames(MarkerUpSummaryDF) = NULL
  colnames(MarkerUpSummaryDF) = gsub("SubCluster",".",colnames(MarkerUpSummaryDF))
  colnames(MarkerUpSummaryDF) = gsub("Quant","Qnt",colnames(MarkerUpSummaryDF) )
  
  
  font.size <- "7.5pt"
  MarkerUpDT = datatable(
    MarkerUpSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
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
    formatStyle(which(grepl("^\\.",colnames(MarkerUpSummaryDF))), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerUpDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/DT/",
      "Compartment",ck,"OverallMarkersUp.html"
    )
    
  )
  
}

# ==== Interactive Table (Markers Down) ====
for(ck in as.character(0:4)){
  tmpLst = ClusterMarkersDown[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$Gene  %in% genes_atleast15,]})
  tmpLst = tmpLst[sapply(tmpLst,nrow) > 0]
  
  tmpLst = lapply(tmpLst,function(x){x[order(x$OverallDiff,decreasing=F),]})
  
  ClusterMarkersDownDF = do.call(rbind,tmpLst)
  rownames(ClusterMarkersDownDF) = NULL
  
  MarkerDownSummaryDF = data.frame(
    ClusterMarkersDownDF,
    scaled_coeff[ClusterMarkersDownDF$Gene,]
  )
  rownames(MarkerDownSummaryDF) = NULL
  colnames(MarkerDownSummaryDF) = gsub("SubCluster",".",colnames(MarkerDownSummaryDF))
  colnames(MarkerDownSummaryDF) = gsub("Quant","Qnt",colnames(MarkerDownSummaryDF) )
  
  
  font.size <- "7.5pt"
  MarkerDownDT = datatable(
    MarkerDownSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
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
    formatStyle(which(grepl("^\\.",colnames(MarkerDownSummaryDF))), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    MarkerDownDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/DT/",
      "Compartment",ck,"OverallMarkersDown.html"
    )
    
  )
  
}

pathway_coeff = comp_path_lfit2$coefficients[,grep("SubCluster",colnames(comp_path_lfit2$coefficients))]

scaled_coeff = t(scale(t(pathway_coeff)))
paletteLength <- 50
myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})


# ==== Interactive Table (Pathways Up) ====
for(ck in as.character(0:4)){
  tmpLst = ClusterPathwaysUp[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$PctExprs >= 10,]})
  ClusterPathwaysUpDF = do.call(rbind,tmpLst)
  rownames(ClusterPathwaysUpDF) = NULL
  
  PathwayUpSummaryDF = data.frame(
    ClusterPathwaysUpDF,
    scaled_coeff[ClusterPathwaysUpDF$Pathway,]
  )
  rownames(PathwayUpSummaryDF) = NULL
  colnames(PathwayUpSummaryDF) = gsub("SubCluster",".",colnames(PathwayUpSummaryDF))
  colnames(PathwayUpSummaryDF) = gsub("Quant","Qnt",colnames(PathwayUpSummaryDF) )
  
  
  font.size <- "7.5pt"
  PathwayUpDT = datatable(
    PathwayUpSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
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
    formatStyle(which(grepl("^\\.",colnames(PathwayUpSummaryDF))), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayUpDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/DT/",
      "Compartment",ck,"OverallPathwaysUp.html"
    )
    
  )
}

# ==== Interactive Table (Pathways Down) ====
for(ck in as.character(0:4)){
  tmpLst = ClusterPathwaysDown[clustSplit[[ck]]]
  tmpLst = lapply(tmpLst,function(x){x[x$Pathway %in% pathways_atleast10,]})
  tmpLst = tmpLst[sapply(tmpLst,nrow) > 0]
  tmpLst = lapply(tmpLst,function(x){x[order(x$OverallDiff,decreasing=F),]})
  
  ClusterPathwaysDownDF = do.call(rbind,tmpLst)
  rownames(ClusterPathwaysDownDF) = NULL
  
  PathwayDownSummaryDF = data.frame(
    ClusterPathwaysDownDF,
    scaled_coeff[ClusterPathwaysDownDF$Pathway,]
  )
  rownames(PathwayDownSummaryDF) = NULL
  colnames(PathwayDownSummaryDF) = gsub("SubCluster",".",colnames(PathwayDownSummaryDF))
  colnames(PathwayDownSummaryDF) = gsub("Quant","Qnt",colnames(PathwayDownSummaryDF) )
  
  
  font.size <- "7.5pt"
  PathwayDownDT = datatable(
    PathwayDownSummaryDF,
    extensions = 'FixedColumns',
    width=1800,
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
    formatStyle(which(grepl("^\\.",colnames(PathwayDownSummaryDF))), backgroundColor = styleInterval(cuts=myBreaks, values=my_rgb))
  
  DT::saveWidget(
    PathwayDownDT, 
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/DT/",
      "Compartment",ck,"OverallPathwaysDown.html"
    )
    
  )
}
#### Here ####


tmptab = liver_metadata[,c("Batch","SubCluster")] %>% 
  table 

tmpdf = as.matrix.data.frame(tmptab)
dimnames(tmpdf) = dimnames(tmptab)
tmpdf = tmpdf %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var="Sample")
tmpdf$Source = liver_metadata$Source[match(tmpdf$Sample,liver_metadata$Batch)]

tmpdf = tmpdf %>% dplyr::relocate("Source")

write_tsv(
  x = tmpdf,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/SampleNucleiCountsBySubCluster.tsv"
  )
)

drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/Clustering.v.1.1/SampleNucleiCountsBySubCluster.tsv"
  ),
  path=as_dribble("https://drive.google.com/drive/u/0/folders/1aIIw9k00gEU1icnowP-PBSAaMLWjpxRq"),
  type = "spreadsheet",
  overwrite = TRUE
)
