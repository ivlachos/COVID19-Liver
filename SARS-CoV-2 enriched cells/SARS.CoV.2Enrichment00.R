# conda activate r_enc
.libPaths("/home/yered/miniconda3/envs/r_enc/lib/R/library")


rm(list=ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(hdf5r)
library(Matrix)
library(parallel)

setwd("/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverGEO")
# ==== H5 Files ====
h5gz_files = system(
  command = 'find /data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverGEO -type f -name "*.h5.gz"',
  intern = TRUE
)

pb = txtProgressBar(
  min = 0, max = length(h5gz_files),style = 3
)

for(ic in seq_along(h5gz_files)){
  # decompress file
  system(
    command = paste0("gzip -d ",h5gz_files[ic]),
  )
  
  # ==== Pre-Cellbender Counts ====
  # get counts
  h5counts = Read10X_h5(
    filename = gsub(".gz","",h5gz_files[ic],fixed=T)
  )
  
  # compress file
  system(
    command = paste0("gzip ",gsub(".gz","",h5gz_files[ic],fixed=T)),
  )
  
  # get sample name
  sname = sapply(strsplit(h5gz_files[ic],"/"),tail,1)
  sname = unlist(strsplit(sname,"_"))[2]
  
  # rename nuclei
  tmpcnames = gsub("\\-[0-9]","",colnames(h5counts))
  colnames(h5counts) =  paste0(sname,"-",tmpcnames)
  
  # ==== Post-Cellbender Counts ====
  postcellbender = readRDS(
    paste0(
      "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/",
      sname,"/",sname,".GRCh38premrna-rna.Seurat.RDS"
    )
  )
  
  keep_nuclei = colnames(postcellbender)
  ambient_nuclei = setdiff(colnames(h5counts),colnames(postcellbender))
  
  # ==== COVID Genes ====
  covid_genes = grep("sars",rownames(h5counts),ignore.case = T,value=T)
  covid_totalcounts = colSums(h5counts[covid_genes,,drop=F])
  totalcounts = colSums(h5counts)
  
  keep_ind = totalcounts > 0
  totalcounts = totalcounts[keep_ind]
  covid_totalcounts = covid_totalcounts[keep_ind]
  
  p_covid = covid_totalcounts/totalcounts
  
  keep_nuclei = keep_nuclei[keep_nuclei %in% names(p_covid)]
  ambient_nuclei = ambient_nuclei[ambient_nuclei %in% names(p_covid)]
  
  # ==== Null Distribution (From Ambient) ====
  p_null_boot = mclapply(
    1:1E5,
    function(ic){
      mean(
        sample(
          x=p_covid[ambient_nuclei],size=length(ambient_nuclei),replace=TRUE
        )
      )
    }
  )
  
  p_null_boot = unlist(p_null_boot)
  
  
  # ==== Bootstrap PValue ====
  boot_pval = mclapply(
    keep_nuclei,
    function(x){
      mean(p_null_boot > p_covid[x] )
    },
    mc.cores =8
  )
  
  boot_pval = unlist(boot_pval)
  boot_pval[boot_pval == 0] = 1/1E6
  
  boot_test = data.frame(
    COVID.UMI = covid_totalcounts[keep_nuclei],
    Total.UMI = totalcounts[keep_nuclei],
    COVID.Abundance = p_covid[keep_nuclei],
    Pval = boot_pval
  )
  
  boot_test$FDR = p.adjust(p=boot_pval,method="fdr")
  
  # save results
  saveRDS(
    boot_test,
    paste0(
      "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/",
      sname,"/",sname,".COVID.BootTest.RDS"
    )
  )
  
  cat("\n\n")
  setTxtProgressBar(pb, ic)
  cat("\n\n")
  
}


