
patientEnrichLst = list()

liver_metadata$COVID.plus= liver_metadata$Boot.SARS.CoV.2 == "RNA.Plus" & liver_metadata$SARS.CoV.2.UMI >= 2

# remove doublets
patient_metadata = liver_metadata[!grepl("DBL",liver_metadata$ClusterName),]


# keep only donors with at least 3 COVID+ nuclei
donor_tab = tapply(X=patient_metadata$COVID.plus,INDEX=patient_metadata$Batch2,FUN = sum)
keep_donors = names(which(donor_tab >= 3))
patient_metadata = patient_metadata[ patient_metadata$Batch2 %in%  keep_donors, ]

AllPatients = sort(unique(patient_metadata$Batch2))

# =====> Enrichment II <====
# ==== All Livers ====
for(my_patient in AllPatients){
  cat(
    paste0(
      "\n===================\n",
      my_patient,
      "\n===================\n"
    )
  )
  
  tmpLst = list()
  # subset by patient
  subset_metadata = patient_metadata[patient_metadata$Batch2 == my_patient,]
  for(ck in as.character(0:4)){
    cat(
      paste0(
        "\n\n------------\n",
        ck,
        "\n------------\n"
      )
    )
    
    compartment_metadata = subset_metadata[subset_metadata$seurat_clusters == ck,]
    
    # ===== Subclusters: Observed Enrichment Scores =====
    AllSubClusters = sort(unique(compartment_metadata$SubCluster))
    total_vcells = sum(compartment_metadata$COVID.plus)
    
    if(total_vcells == 0){
      clusterDF = data.frame(
        DonorID = donor_names[my_patient],
        Cluster = AllSubClusters,
        Compartment = ck,
        ClusterName = short_names[AllSubClusters],
        CompartmentName = compartment_names[ck],
        ES = NA,
        NES = NA,
        Pval = NA,
        FDR = NA
      )
    }
    
    if(total_vcells > 0){
      
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
      
      # mypar(1,1)
      # barplot(clust_scores,las=2)
      
      
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
      
      std_clust_scores[is.na(std_clust_scores)] = 0
      
      clusterDF = data.frame(
        DonorID = donor_names[my_patient],
        Cluster = AllSubClusters,
        Compartment = ck,
        ClusterName = short_names[AllSubClusters],
        CompartmentName = compartment_names[ck],
        ES = clust_scores[AllSubClusters],
        NES = std_clust_scores[AllSubClusters],
        Pval = p_boot[AllSubClusters],
        FDR = fdr_boot[AllSubClusters]
      )
      
      tmpLst[[ck]] = clusterDF
      
      
      cat("\n\nDONE!\n\n")
    }
    
    patientEnrichLst[[my_patient]] = tmpLst
  }
}

patientEnrichLst = lapply(patientEnrichLst,function(x){
  x = do.call(rbind,x)
  rownames(x) = NULL
  return(x)
})


saveRDS(
  patientEnrichLst,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_",
    "Clusteringv1.1_ByPatient_CovidEnrichment04.RDS"
  )
)


