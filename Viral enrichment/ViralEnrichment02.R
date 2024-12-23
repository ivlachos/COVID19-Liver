enrichLst = list()
# =====> Enrichment II <====
for(ck in as.character(0:4)){
  cat(
    paste0(
      "\n\n------------\n",
      ck,
      "\n------------\n"
    )
  )
  
  compartment_metadata = liver_metadata[liver_metadata$seurat_clusters == ck,]
  compartment_metadata$COVID.plus = compartment_metadata$Boot.SARS.CoV.2 == "RNA.Plus" & compartment_metadata$SARS.CoV.2.UMI >= 2
  
  # remove doublets
  compartment_metadata = compartment_metadata[!grepl("DBL",compartment_metadata$ClusterName),]
  
  # keep only donors with at least 3 COVID+ nuclei
  donor_tab = tapply(X=compartment_metadata$COVID.plus,INDEX=compartment_metadata$Batch2,FUN = sum)
  keep_donors = names(which(donor_tab >= 3))
  compartment_metadata = compartment_metadata[ compartment_metadata$Batch2 %in%  keep_donors, ]
  
  # ===== Subclusters: Observed Enrichment Scores =====
  AllSubClusters = sort(unique(compartment_metadata$SubCluster))
  total_vcells = sum(compartment_metadata$COVID.plus)
  clust_vcells_vec = tapply(
    X=compartment_metadata$COVID.plus,
    INDEX=compartment_metadata$SubCluster,
    FUN=sum
  )
  clust_prop_vec = prop.table(table(compartment_metadata$SubCluster))
  
  clust_vcells_vec = clust_vcells_vec[AllSubClusters]
  clust_prop_vec = clust_prop_vec[AllSubClusters]
  
  clust_scores = sapply(
    AllSubClusters,
    function(x){
      EnrichScore(
        clust_vcells = clust_vcells_vec[x],
        total_vcells = total_vcells,
        clust_prop = clust_prop_vec[x]
      )
    }
  )
  names(clust_scores) = AllSubClusters
  
  mypar(1,1,mar=c(7,2,1,1))
  barplot(clust_scores,las=2)
  
  
  # ===== Subclusters: Permuted Enrichment Scores =====
  perm_enrich_scores = mclapply(
    1:1E5,
    function(ic){
      boot_status = rep("Other",nrow(compartment_metadata))
      boot_ind = sample(x=seq_along(boot_status),size = total_vcells)
      boot_status[boot_ind] = "RNA.Plus"
      
      clust_vcells_boot = tapply(
        X=boot_status == "RNA.Plus",
        INDEX=compartment_metadata$SubCluster,
        FUN=sum
      )
      
      clust_scores_boot = sapply(
        AllSubClusters,
        function(x){
          EnrichScore(
            clust_vcells = clust_vcells_boot[x],
            total_vcells = total_vcells,
            clust_prop = clust_prop_vec[x]
          )
        }
      )
      names(clust_scores_boot) = AllSubClusters
      return(clust_scores_boot)
    },
    mc.cores=8
  )
  
  perm_enrich_scores = do.call(rbind,perm_enrich_scores)
  
  
  p_boot = sapply(
    AllSubClusters,
    function(x){
      mean(perm_enrich_scores[,x,drop=TRUE] > clust_scores[x] )
    }
  )
  
  
  fdr_boot = p.adjust(p=p_boot,method="fdr")
  
  
  
  mypar(a=ceiling(length(AllSubClusters)/4),b=4)
  for(x in AllSubClusters){
    xrng = range(c(clust_scores[x],perm_enrich_scores[,x,drop=TRUE]))
    hist(perm_enrich_scores[,x,drop=TRUE],xlim=xrng,main=x,xlab="")
    abline(v=clust_scores[x],col="red",lwd=2)
  }
  
  
  perm_means = colMeans(perm_enrich_scores)
  perm_sd = colSds(perm_enrich_scores)
  names(perm_sd) = colnames(perm_enrich_scores)
  
  std_clust_scores = sapply(
    AllSubClusters,
    function(x){
      (clust_scores[x] - perm_means[x])/perm_sd[x]
    }
  )
  names(std_clust_scores) = AllSubClusters
  
  clusterDF = data.frame(
    Cluster = AllSubClusters,
    Compartment = ck,
    ClusterName = short_names[AllSubClusters],
    CompartmentName = compartment_names[ck],
    ES = clust_scores[AllSubClusters],
    NES = std_clust_scores[AllSubClusters],
    Pval = p_boot[AllSubClusters],
    FDR = fdr_boot[AllSubClusters]
  )
  
  
  mypar(1,1,mar=c(5,2,2,1))
  barplot(
    clusterDF$NES,
    names.arg = clusterDF$Cluster,
    las=2,
    col = ifelse(
      clusterDF$FDR < 0.01,
      "red","gray"
    ),
    main=paste0("SARS-CoV-2 RNA+ Enrichment Score (",ck,")")
  )
  
  enrichLst[[ck]] = clusterDF
  
  cat("\n\nDONE!\n\n")
}

enrichDF = do.call(rbind,enrichLst)
rownames(enrichDF) = NULL


saveRDS(
  enrichDF,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Clusteringv1.1_CovidEnrichment04.RDS"
  )
)


write_tsv(
  x = enrichDF,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Clusteringv1.1_CovidEnrichment04.tsv"
  )
)


drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Clusteringv1.1_CovidEnrichment04.tsv"
  ),
  path=as_dribble("https://drive.google.com/drive/u/0/folders/1eIQScIrui5UqucrBPuVcUbmMaChN96G9"),
  type = "spreadsheet",
  overwrite = TRUE
)
