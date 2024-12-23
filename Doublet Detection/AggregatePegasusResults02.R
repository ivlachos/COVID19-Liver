rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(data.table)


# ====> ncLab Livers <=====
liver1_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver1/Liver1.Pegasus.csv"
)
liver1_pegasus$Batch = "ncLabLiver1"

liver2_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver2_and_2a/Liver2_and_2a.Pegasus.csv"
)
liver2_pegasus$Batch = "ncLabLiver2"

liver3_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver3/Liver3.Pegasus.csv"
)
liver3_pegasus$Batch = "ncLabLiver3"

# ====> Broad Livers <=====
csv_files = system(
  command = 'find /data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus -type f -name "*.Pegasus.csv"',
  intern = TRUE
)
batch_labels = gsub(".Pegasus.csv","",sapply(strsplit(csv_files,"/"),tail,1))

pegasusLst = lapply(
  csv_files,
  fread
)

for(ic in seq_along(pegasusLst)){
  pegasusLst[[ic]]$Batch = batch_labels[ic]
}


colsLst = lapply(pegasusLst,colnames)
common_cols = Reduce(intersect,colsLst)

pegasusLst = lapply(
  pegasusLst,
  function(x){x[,..common_cols]}
)


liverBroadPegasus = do.call(
  rbind,
  pegasusLst
)


# ====> Ben Izar <====
liver_02_cov_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_02_cov_Ben/liver_02_cov_Ben.Pegasus.csv"
)
liver_02_cov_pegasus$Batch = "Izarliver_02_cov"

liver_12_cov_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_12_cov_Ben/liver_12_cov_Ben.Pegasus.csv"
)
liver_12_cov_pegasus$Batch = "Izarliver_12_cov"

liver_01_cov_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_01_cov_Ben/liver_01_cov_Ben.Pegasus.csv"
)
liver_01_cov_pegasus$Batch = "Izarliver_01_cov"

liver_08_cov_pegasus = fread(
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_08_cov_Ben/liver_08_cov_Ben.Pegasus.csv"
)
liver_08_cov_pegasus$Batch = "Izarliver_08_cov"

# ==== Merge Pegasus Results =====
select_cols = c("barcodekey","doublet_score", "pred_dbl","Batch" )

PegasusRes = rbind(
  liver1_pegasus[,..select_cols],
  liver2_pegasus[,..select_cols],
  liver3_pegasus[,..select_cols],
  liverBroadPegasus[,..select_cols],
  liver_02_cov_pegasus[,..select_cols],
  liver_12_cov_pegasus[,..select_cols],
  liver_01_cov_pegasus[,..select_cols],
  liver_08_cov_pegasus[,..select_cols]
)



saveRDS(
  PegasusRes,
  "/data/work/Projects/BrScRNAseq/data/BroadTerra/CumulusCellbenderLiverAllPegasusScoresBatch.RDS"
)
