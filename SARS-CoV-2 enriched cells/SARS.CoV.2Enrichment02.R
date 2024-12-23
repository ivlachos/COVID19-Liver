rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(Matrix)
library(parallel)


# ==== Liver Data ====
liver_metadata = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsMetadata6.RDS"
)

bootDF = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversCOVIDBootTest.RDS"
)

covid_counts = readRDS(
  "/data/shiny/shiny-server/apps/SelectLiverNoDblCOVID/data/SelectLiverNoDoubletsCOVIDcounts.RDS"
)

tmpLabs = unique(liver_metadata$SubCluster)
tmpLabs = tmpLabs[sapply(strsplit(tmpLabs,"_"),length) == 3]
tmpFac = sapply(strsplit(tmpLabs,"_"),function(x){paste(head(x,2),collapse="_")})

split_subclust = split(
  x=tmpLabs,
  f=tmpFac
)

missing_nuclei = rownames(liver_metadata)[which(!(rownames(liver_metadata) %in% bootDF$Nuclei))]

liver_metadata[missing_nuclei,]
covid_counts[,missing_nuclei]
# ==== Liver 1: Seurat Objects ====
# Pre-Cellbender
precellbender = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/pegasus_BI_liver_trimmed/Liver1/Liver1.GRCh38premrna-rna.Seurat.RDS"
)
# Post-Cellbender
postcellbender = readRDS(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver1/Liver1.GRCh38premrna-rna.Seurat.RDS"
)

# ==== COVID Genes ====
precovid_genes = grep("sars",rownames(precellbender@assays$RNA@counts),ignore.case = T,value=T)
postcovid_genes = grep("sars",rownames(postcellbender@assays$RNA@counts),ignore.case = T,value=T)

preclbdnr_covid_totalcounts = colSums(precellbender@assays$RNA@counts[precovid_genes,])
postclbdnr_covid_totalcounts = colSums(postcellbender@assays$RNA@counts[postcovid_genes,])

preclbdnr_totalcounts = colSums(precellbender@assays$RNA@counts)
postclbdnr_totalcounts = colSums(postcellbender@assays$RNA@counts)

ambient_nuclei = setdiff(names(preclbdnr_totalcounts),colnames(postcellbender))

rm(precellbender,postcellbender)

# estimate null p from ambient 
p_covid = preclbdnr_covid_totalcounts/preclbdnr_totalcounts
p_missing = postclbdnr_covid_totalcounts[missing_nuclei]/postclbdnr_totalcounts[missing_nuclei]

# ==== Null Distribution (From Ambient) ====
p_null_boot = mclapply(
  1:1E5,
  function(ic){
    mean(
      sample(
        x=p_covid[ambient_nuclei],size=length(ambient_nuclei),replace=TRUE
      )
    )
  },mc.cores = 8
)

p_null_boot = unlist(p_null_boot)

# ==== Bootstrap PValue ====
boot_pval = mclapply(
  missing_nuclei,
  function(x){
    mean(p_null_boot > p_missing[x] )
  },
  mc.cores =8
)

boot_pval = unlist(boot_pval)
boot_pval[boot_pval == 0] = 1/1E6

boot_missing = data.frame(
  Nuclei = missing_nuclei,
  COVID.UMI = postclbdnr_covid_totalcounts[missing_nuclei],
  Total.UMI = postclbdnr_totalcounts[missing_nuclei],
  COVID.Abundance = p_missing[missing_nuclei],
  Pval = boot_pval
)

boot_test = rbind(
  boot_missing,bootDF[bootDF$Batch == "ncLabLiver1",colnames(boot_missing)]
)

boot_test$FDR = p.adjust(p=boot_test$Pval,method="fdr")
boot_test$Batch = "ncLabLiver1"
rownames(boot_test) = NULL

# update results
resDF = rbind(
  boot_test,
  bootDF[bootDF$Batch != "ncLabLiver1",]
)

saveRDS(
  resDF,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversCOVIDBootTestFull.RDS"
)
