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


library(RColorBrewer)
library(DT)
Sys.setenv("RSTUDIO_PANDOC" = "/data/resources/tools/pandoc/")




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


# ===== Liver Data =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1.RDS"
)
liver_metadata$lognUMI = log(liver_metadata$nUMI)
liver_metadata$lognGene = log(liver_metadata$nGene)

liver_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDbl2a/data/SelectLiverNoDoubletsSeuratCounts.RDS"
)

# Pathway scores
comp_pathwayScoresAdj = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblPathways/data/AdjustedAllPathwayScoresAdjEmpZRanks.RDS"
)

base_obj = ls()

for(ck in as.character(0:4)){ 
  
  cat(
    paste0(
      "\n-----------------\n",
      ck,
      "\n-----------------\n"
    )
  )
  
  subset_metadata = liver_metadata[liver_metadata$seurat_clusters == ck,]
  AllSubClusters = sort(unique(subset_metadata$SubCluster))
  
  cat("\nFitting Pathway Linear Model ...\n")
  # ==== Pathway Linear Model =====
  pathway_design2 <- model.matrix(~ -1 + SubCluster, subset_metadata )
  comp_path_lfit2 = lmFit(
    object = comp_pathwayScoresAdj[,rownames(subset_metadata)],
    design = pathway_design2
  )
  comp_path_lfit2 = eBayes(comp_path_lfit2)
  rm(pathway_design2)
  
  saveRDS(
    comp_path_lfit2,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwayLm.RDS"
    )
  )
  
  
  cat("\nFitting Gene Linear Model ...\n")
  # ===== Gene Linear Model  =====
  # filter genes
  subset_counts = liver_counts[,rownames(subset_metadata)]
  subpct = mclapply(
    AllSubClusters,
    function(sck){
      Matrix::rowMeans(subset_counts[,subset_metadata$SubCluster == sck] > 0)
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
    counts=subset_counts[keep_ind,rownames(subset_metadata)], 
    group=subset_metadata$SubCluster
  )
  compartment_dge <- calcNormFactors(object=compartment_dge,method="TMM")
  
  compartment_design2 <- model.matrix(~ -1 + SubCluster + Batch, subset_metadata )
  colnames(compartment_design2) = make.names(colnames(compartment_design2) )
  
  y <- new("EList")
  y$E <- edgeR::cpm(compartment_dge, log = TRUE, prior.count = 3)
  compartment_fit2 <- lmFit(y, design = compartment_design2)
  compartment_fit2 <- eBayes(compartment_fit2, trend = TRUE, robust = TRUE)
  
  rm(compartment_design2,y,compartment_dge)
  
  saveRDS(
    compartment_fit2,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_GeneLm.RDS"
    )
  )
  
  cat("\nComputing expression summaries ...\n")
  # ==== Expression Summaries ====
  # - percent expressed
  CountsPctExprs = mclapply(
    AllSubClusters,
    function(sck){
      CurrentCells = rownames(subset_metadata)[subset_metadata$SubCluster == sck]
      rowMeans(subset_counts[,CurrentCells] > 0)*100
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
  
  cat("\nComputing pathway summaries ... \n")
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
      colMeans(CountsPctExprs[rownames(CountsPctExprs) %in% x,,drop=F])
    },
    mc.cores = 12
  )
  PathwayPctExprs = do.call(rbind,PathwayPctExprs)
  
  
  
  
  
  cat("\nGetting markers ...\n")
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
      
      
      if( length(subclust_up_markers) > 0 ){
        tmptab = cbind(
          subclust_contrast_toptab[subclust_up_markers,1:2],
          subclust_contrast2_toptab[subclust_up_markers,c("F","P.Value","adj.P.Val")]
        )
        
        
        tmptab = tmptab[order(tmptab[,1],decreasing=T),]
        
        tmptab = tibble::rownames_to_column(tmptab,var="Gene")
      }else{
        tmptab = NULL
      }
      
      
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
      
      if( length(subclust_down_markers) > 0 ){
        tmptab = cbind(
          subclust_contrast_toptab[subclust_down_markers,1:2],
          subclust_contrast2_toptab[subclust_down_markers,c("F","P.Value","adj.P.Val")]
        )
        
        
        tmptab = tmptab[order(tmptab[,1],decreasing=T),]
        tmptab = tibble::rownames_to_column(tmptab,var="Gene")
      }else{
        tmptab = NULL
      }
      
      
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
      
      if(length(subclust_up_markers) > 0){
        tmptab = cbind(
          subclust_contrast_toptab[subclust_up_markers,1:2],
          subclust_contrast2_toptab[subclust_up_markers,c("F","P.Value","adj.P.Val")]
        )
        
        tmptab = tmptab[order(tmptab[,1],decreasing=T),]
        tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
      }else{
        tmptab = NULL
      }
      
      
      
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
      
      if( length(subclust_down_markers) > 0 ){
        tmptab = cbind(
          subclust_contrast_toptab[subclust_down_markers,1:2],
          subclust_contrast2_toptab[subclust_down_markers,c("F","P.Value","adj.P.Val")]
        )
        
        tmptab = tmptab[order(tmptab[,1],decreasing=T),]
        tmptab = tibble::rownames_to_column(tmptab,var="Pathway")
      }else{
        tmptab = NULL
      }
      
      
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
  
  saveRDS(
    resubclustMrkUpLst,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersUpLst.RDS"
    )
  )
  saveRDS(
    resubclustMrkDownLst,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersDownLst.RDS"
    )
  )
  
  saveRDS(
    resubclustPathUpLst,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysUpLst.RDS"
    )
  )
  
  saveRDS(
    resubclustPathDownLst,
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysDownLst.RDS"
    )
  )
  
  resubclustMrkUpLst = resubclustMrkUpLst[!sapply(resubclustMrkUpLst,is.null)]
  resubclustMrkDownLst = resubclustMrkDownLst[!sapply(resubclustMrkDownLst,is.null)]
  resubclustPathUpLst = resubclustPathUpLst[!sapply(resubclustPathUpLst,is.null)]
  resubclustPathDownLst = resubclustPathDownLst[!sapply(resubclustPathDownLst,is.null)]
  

  
  cat("\nFormatting markers ...\n")
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersUp.RDS"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersDown.RDS"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysUp.RDS"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysDown.RDS"
    )
  )
  
  rm_obj = setdiff(ls(),base_obj)
  rm_obj = rm_obj[rm_obj != "base_obj"]
  rm(list=rm_obj)
  
  cat("\n\nDONE!\n\n")
}


# ===== Upload Results =====
library(googledrive)

dribblePaths = c(
  "0" = "https://drive.google.com/drive/u/0/folders/1kjxLmuNnGgR3qJ7CnMY3UsWtyMvx9qj1",
  "1" = "https://drive.google.com/drive/u/0/folders/1siF87RvpeM4-lZWIpbhYXFu7QUROI5Cr",
  "2" = "https://drive.google.com/drive/u/0/folders/1gF_xL_6JYRBLEh2mAlDx4sG662FTY4Nw",
  "3" = "https://drive.google.com/drive/u/0/folders/1li9eBkyTZSu1Op1xnAKXcZcRPJEbTLI8",
  "4" = "https://drive.google.com/drive/u/0/folders/1Ef70FPEvS-plo86Xfc99SheDLj3_jR9l"
)


for(ck in as.character(1:4)){
  ClusterMarkersUp = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersUp.RDS"
    )
  )
  
  ClusterMarkersUp = do.call(rbind,ClusterMarkersUp)
  rownames(ClusterMarkersUp) = NULL
  
  write_tsv(
    x = ClusterMarkersUp,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinMarkersUp.tsv"
    )
  )
  
  drive_upload(
    media =  paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinMarkersUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  
  ClusterMarkersDown = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersDown.RDS"
    )
    
  )
  
  ClusterMarkersDown = do.call(rbind,ClusterMarkersDown)
  rownames(ClusterMarkersDown) = NULL
  
  write_tsv(
    x = ClusterMarkersDown,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinMarkersDown.tsv"
    )
  )
  
  drive_upload(
    media =  paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinMarkersDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  ClusterPathwaysUp = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysUp.RDS"
    )
  )
  
  ClusterPathwaysUp = do.call(rbind,ClusterPathwaysUp)
  rownames(ClusterPathwaysUp) = NULL
  
  write_tsv(
    x = ClusterPathwaysUp,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinPathwaysUp.tsv"
    )
  )
  
  drive_upload(
    media =  paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinPathwaysUp.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
  
  
  ClusterPathwaysDown = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysDown.RDS"
    )
    
  )
  
  ClusterPathwaysDown = do.call(rbind,ClusterPathwaysDown)
  rownames(ClusterPathwaysDown) = NULL
  
  write_tsv(
    x = ClusterPathwaysDown,
    file = paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinPathwaysDown.tsv"
    )
  )
  
  drive_upload(
    media =  paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_WithinPathwaysDown.tsv"
    ),
    path=as_dribble(dribblePaths[ck]),
    type = "spreadsheet",
    overwrite = TRUE
  )
}

for(ck in as.character(0:4)){
  
  compartment_fit2 = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_GeneLm.RDS"
    )
  )
  
  subclust_coeff = compartment_fit2$coefficients[,grep("SubCluster",colnames(compartment_fit2$coefficients))]
  scaled_coeff = t(scale(t(subclust_coeff)))
  paletteLength <- 50
  myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
  my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})
  
  # ==== Interactive Table (Markers Up) ====
  ClusterMarkersUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersUp.RDS"
    )
  )
  
  ClusterMarkersUpDF = do.call(rbind,ClusterMarkersUpDF)
  rownames(ClusterMarkersUpDF) = NULL
  
  MarkerUpSummaryDF = data.frame(
    ClusterMarkersUpDF,
    scaled_coeff[ClusterMarkersUpDF$Gene,]
  )
  rownames(MarkerUpSummaryDF) = NULL
  colnames(MarkerUpSummaryDF) = gsub("SubCluster",".",colnames(MarkerUpSummaryDF))
  colnames(MarkerUpSummaryDF) = gsub("Quant","Qnt",colnames(MarkerUpSummaryDF) )
  
  font.size <- "7pt"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/DT/",
      "Compartment",ck,"WithinMarkersUp.html"
    )
    
  )
  
  
  
  
  
  # ==== Interactive Table (Markers Down) ====
  ClusterMarkersDownDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_MarkersDown.RDS"
    )
  )
  
  ClusterMarkersDownDF = do.call(rbind,ClusterMarkersDownDF)
  rownames(ClusterMarkersDownDF) = NULL
  
  MarkerDownSummaryDF = data.frame(
    ClusterMarkersDownDF,
    scaled_coeff[ClusterMarkersDownDF$Gene,]
  )
  rownames(MarkerDownSummaryDF) = NULL
  colnames(MarkerDownSummaryDF) = gsub("SubCluster",".",colnames(MarkerDownSummaryDF))
  colnames(MarkerDownSummaryDF) = gsub("Quant","Qnt",colnames(MarkerDownSummaryDF) )
  
  font.size <- "7pt"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/DT/",
      "Compartment",ck,"WithinMarkersDown.html"
    )
    
  )
  
  comp_path_lfit2 = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwayLm.RDS"
    )
  )
  
  pathway_coeff = comp_path_lfit2$coefficients[,grep("SubCluster",colnames(comp_path_lfit2$coefficients))]
  
  scaled_coeff = t(scale(t(pathway_coeff)))
  paletteLength <- 50
  myBreaks <- c(seq(min(scaled_coeff), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(scaled_coeff)/paletteLength, max(scaled_coeff), length.out=floor(paletteLength/2)))
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))( length(myBreaks) + 1  )
  my_rgb = apply(col2rgb(myColor) - 1,2,function(x){paste0("rgb(",x[1],",",x[2],",",x[3],")")})
  
  # ==== Interactive Table (Pathways Up) ====
  ClusterPathwaysUpDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysUp.RDS"
    )
  )
  
  ClusterPathwaysUpDF = do.call(rbind,ClusterPathwaysUpDF)
  rownames(ClusterPathwaysUpDF) = NULL
  
  PathwayUpSummaryDF = data.frame(
    ClusterPathwaysUpDF,
    scaled_coeff[ClusterPathwaysUpDF$Pathway,]
  )
  rownames(PathwayUpSummaryDF) = NULL
  colnames(PathwayUpSummaryDF) = gsub("SubCluster",".",colnames(PathwayUpSummaryDF))
  colnames(PathwayUpSummaryDF) = gsub("Quant","Qnt",colnames(PathwayUpSummaryDF) )
  
  font.size <- "7pt"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/DT/",
      "Compartment",ck,"WithinPathwaysUp.html"
    )
    
  )
  
  
  # ==== Interactive Table (Pathways Down) ====
  ClusterPathwaysDownDF = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/",
      "Compartment",ck,"_PathwaysDown.RDS"
    )
  )
  
  ClusterPathwaysDownDF = do.call(rbind,ClusterPathwaysDownDF)
  rownames(ClusterPathwaysDownDF) = NULL
  
  PathwayDownSummaryDF = data.frame(
    ClusterPathwaysDownDF,
    scaled_coeff[ClusterPathwaysDownDF$Pathway,]
  )
  rownames(PathwayDownSummaryDF) = NULL
  colnames(PathwayDownSummaryDF) = gsub("SubCluster",".",colnames(PathwayDownSummaryDF))
  colnames(PathwayDownSummaryDF) = gsub("Quant","Qnt",colnames(PathwayDownSummaryDF) )
  
  font.size <- "7pt"
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
      "SelectLiversNoDoublets/Clustering.v.1.1/WithinCompartment/DT/",
      "Compartment",ck,"WithinPathwaysDown.html"
    )
    
  )
  
}

