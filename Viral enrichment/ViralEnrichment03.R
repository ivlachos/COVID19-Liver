# ==== All Livers ====

tmp00 = liver_metadata[,c("DonorID","Batch2")] %>% table()
donor_names = rownames(tmp00)
names(donor_names) = colnames(tmp00)[apply(tmp00,1,function(x){which(x > 0)})]
rm(tmp00)

donor_names = donor_names[order(as.numeric(gsub("L","",donor_names)))]
# ===== Observed Enrichment Scores =====
liver_metadata$COVID.plus= liver_metadata$Boot.SARS.CoV.2 == "RNA.Plus" & liver_metadata$SARS.CoV.2.UMI >= 2

# remove doublets
subset_metadata = liver_metadata[!grepl("DBL",liver_metadata$ClusterName),]

tmp00 = subset_metadata
tmp00$ClusterName = factor(tmp00$ClusterName, levels=short_names)
covid_umi_tab = lapply(
  split(
    x=tmp00,
    f=tmp00$DonorID
  ),
  function(x){
    tapply(X=x$COVID.plus,INDEX=x$ClusterName,FUN=sum,narm=T)
  }
)

covid_umi_tab = do.call(cbind,covid_umi_tab)
covid_umi_tab[is.na(covid_umi_tab)] = 0
covid_umi_tab = covid_umi_tab[,donor_names]

write_tsv(
  x = covid_umi_tab %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var="DonorID"),
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Donor_Clusteringv1.1_CovidUMICounts.tsv"
  )
)


drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Donor_Clusteringv1.1_CovidUMICounts.tsv"
  ),
  path=as_dribble("https://drive.google.com/drive/u/0/folders/1eIQScIrui5UqucrBPuVcUbmMaChN96G9"),
  type = "spreadsheet",
  overwrite = TRUE
)


# keep only donors with at least 3 COVID+ nuclei
donor_tab = tapply(X=subset_metadata$COVID.plus,INDEX=subset_metadata$Batch2,FUN = sum)
keep_donors = names(which(donor_tab >= 3))
subset_metadata = subset_metadata[ subset_metadata$Batch2 %in%  keep_donors, ]

AllPatients = sort(unique(subset_metadata$Batch2))

total_vcells = sum(subset_metadata$COVID.plus)
clust_vcells_vec = tapply(
  X=subset_metadata$COVID.plus,
  INDEX=subset_metadata$Batch2,
  FUN=sum
)
clust_prop_vec = prop.table(table(subset_metadata$Batch2))

clust_vcells_vec = clust_vcells_vec[AllPatients]
clust_prop_vec = clust_prop_vec[AllPatients]

clust_scores = sapply(
  AllPatients,
  function(x){
    EnrichScore(
      clust_vcells = clust_vcells_vec[x],
      total_vcells = total_vcells,
      clust_prop = clust_prop_vec[x]
    )
  }
)
names(clust_scores) = AllPatients

# ===== Permuted Enrichment Scores =====
perm_enrich_scores = mclapply(
  1:1E5,
  function(ic){
    boot_status = rep("Other",nrow(subset_metadata))
    boot_ind = sample(x=seq_along(boot_status),size = total_vcells)
    boot_status[boot_ind] = "RNA.Plus"
    
    clust_vcells_boot = tapply(
      X=boot_status == "RNA.Plus",
      INDEX=subset_metadata$Batch2,
      FUN=sum
    )
    
    clust_scores_boot = sapply(
      AllPatients,
      function(x){
        EnrichScore(
          clust_vcells = clust_vcells_boot[x],
          total_vcells = total_vcells,
          clust_prop = clust_prop_vec[x]
        )
      }
    )
    names(clust_scores_boot) = AllPatients
    return(clust_scores_boot)
  },
  mc.cores=8
)

perm_enrich_scores = do.call(rbind,perm_enrich_scores)

p_boot = sapply(
  AllPatients,
  function(x){
    mean(perm_enrich_scores[,x,drop=TRUE] > clust_scores[x] )
  }
)


fdr_boot = p.adjust(p=p_boot,method="fdr")


perm_means = colMeans(perm_enrich_scores)
perm_sd = colSds(perm_enrich_scores)
names(perm_sd) = colnames(perm_enrich_scores)

std_clust_scores = sapply(
  AllPatients,
  function(x){
    (clust_scores[x] - perm_means[x])/perm_sd[x]
  }
)
names(std_clust_scores) = AllPatients

BatchDF02a = data.frame(
  DonorID = donor_names[AllPatients],
  ES = clust_scores[AllPatients],
  NES = std_clust_scores[AllPatients],
  Pval = p_boot[AllPatients],
  FDR = fdr_boot[AllPatients]
)
BatchDF02a = BatchDF02a[order(as.numeric(gsub("L","",BatchDF02a$DonorID))),]
rownames(BatchDF02a) = NULL

clust_vcells_vec


saveRDS(
  BatchDF02a,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Donor_CovidEnrichment04.RDS"
  )
)


write_tsv(
  x = BatchDF02a,
  file = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Donor_CovidEnrichment04.tsv"
  )
)


drive_upload(
  media = paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "SelectLiversNoDoublets_Donor_CovidEnrichment04.tsv"
  ),
  path=as_dribble("https://drive.google.com/drive/u/0/folders/1eIQScIrui5UqucrBPuVcUbmMaChN96G9"),
  type = "spreadsheet",
  overwrite = TRUE
)


