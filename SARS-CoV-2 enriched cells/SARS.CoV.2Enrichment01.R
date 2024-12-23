rm(list=ls())
options(stringsAsFactors = FALSE)

# ==== Bootstrap COVID Abundance Test Results ====
boot_files = c(
  ## ncLabBI
  "ncLabLiver1" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver1/Liver1.COVID.BootTest.RDS",
  "ncLabLiver2" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver2_and_2a/Liver2_and_2a.COVID.BootTest.RDS",
  "ncLabLiver3" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_BI_liver_trimmed/Liver3/Liver3.COVID.BootTest.RDS",
  ## BroadBI
  "BroadLiver1" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/12-P230638-S003-R01/12-P230638-S003-R01.COVID.BootTest.RDS",
  "BroadLiver2" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/12-P485759-S020-R01/12-P485759-S020-R01.COVID.BootTest.RDS",
  "BroadLiver3" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/12-P617758-S003-R01/12-P617758-S003-R01.COVID.BootTest.RDS",
  "BroadLiver4" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/12-P852049-S003-R01/12-P852049-S003-R01.COVID.BootTest.RDS",
  "BroadLiver5" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/12-P890292-S003-R01/12-P890292-S003-R01.COVID.BootTest.RDS",
  ## BWH
  "04-P054921-S080-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/04-P054921-S080-R01/04-P054921-S080-R01.COVID.BootTest.RDS",
  "04-P079042-S079-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/04-P079042-S079-R01/04-P079042-S079-R01.COVID.BootTest.RDS",
  ## MGH
  "02-P005175-S022-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P005175-S022-R01/02-P005175-S022-R01.COVID.BootTest.RDS",
  "02-P118946-S070-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P118946-S070-R01/02-P118946-S070-R01.COVID.BootTest.RDS",
  "02-P166169-S022-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P166169-S022-R01/02-P166169-S022-R01.COVID.BootTest.RDS",
  "02-P240970-S016-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P240970-S016-R01/02-P240970-S016-R01.COVID.BootTest.RDS",
  "02-P248880-S041-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P248880-S041-R01/02-P248880-S041-R01.COVID.BootTest.RDS",
  "02-P348762-S041-R01" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/BroadLiverCumulus/02-P348762-S041-R01/02-P348762-S041-R01.COVID.BootTest.RDS",
  ## Izar
  "Izarliver_12_cov" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_12_cov_Ben/liver_12_cov_Ben.COVID.BootTest.RDS",
  "Izarliver_01_cov" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_01_cov_Ben/liver_01_cov_Ben.COVID.BootTest.RDS",
  "Izarliver_08_cov" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_08_cov_Ben/liver_08_cov_Ben.COVID.BootTest.RDS",
  "Izarliver_02_cov" = "/data/work/Projects/BrScRNAseq/data/BroadTerra/cellbenderv2_pegasus_ben_and_BI/liver_02_cov_Ben/liver_02_cov_Ben.COVID.BootTest.RDS"
)


bootLst = lapply(
  boot_files,readRDS
)
# add batch information
for(ic in seq_along(bootLst)){
  bootLst[[ic]]$Batch = names(bootLst)[ic]
  bootLst[[ic]] = tibble::rownames_to_column( bootLst[[ic]],var="Nuclei")
}

bootDF = do.call(rbind,bootLst)
rownames(bootDF) = NULL

saveRDS(
  bootDF,
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversCOVIDBootTest.RDS"
)
