###################################################################################################
## This Rscript defines Cell-cell communication with CellChat

###################################################################################################
## Installing and loading required packages
###################################################################################################
###################################################################################################

library("Seurat")
library("data.table")
library("writexl")
library("future")
library("NMF")
library("patchwork")
library("ggalluvial")
library("Matrix")
library("CellChat")

run_cellChat <- function(cellchat, meta, trim_var=trim_var[1], result_dir){
  
  tmp <- paste0(result_dir, "/", trim_var)
  dir.create(tmp, recursive = T)
  
  ################### Communication probability and infer cellular communication network ############
  
  options(future.globals.maxSize= 700000000)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = trim_var)  
 
  
  ############ Cellular communication network as a data frame 
  
  df.ligand.receptor <- subsetCommunication(cellchat)
  df.signaling.pathway <- subsetCommunication(cellchat, slot.name = "netP")
  
  
  write_xlsx(df.ligand.receptor,paste0(tmp, "/ligand_receptor_", trim_var, ".xlsx"))
  write_xlsx(df.signaling.pathway,paste0(tmp, "/signaling_pathway_", trim_var, ".xlsx"))
  
  ### Infer the cell-cell communication at a signaling pathway level
  #NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
  cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
  
  
  ############################ Aggregated cell-cell communication network ############################
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  saveRDS(cellchat, file = paste0(tmp,"/cellchat_COVID19_Liver_", trim_var, ".rds"))
  
  
  return(cellchat)
}



###################################################################################################
#Initialize Seurat object
custom_counts <- readRDS("SelectLiverNoDoubletsCustomSeqCounts.RDS")
metadata <- readRDS("/Data/COVID19LiverMetadata.RDS")
metadata$Cell.id <- rownames(metadata)

annotation <- setnames(fread("/Data/Cell_populations.txt"),
                       c("SubCluster","Label","Alias","Major_compartments"))

metadata <- merge(metadata, annotation, by = "SubCluster")
rownames(metadata) <- metadata$Cell.id

sarsCov2 <- readRDS("SelectLiversNoDoublets_Clusteringv1.1_SARS_Cov2_Enrichment_v1.0.RDS")
sarsCov2$Cell.id <- rownames(sarsCov2)

metadata <- merge(setDT(metadata), setDT(sarsCov2), by="Cell.id")

metadata[SARS.CoV.2.Plus == "TRUE" & Major_compartments == "Hepatocytes"]$Alias <- "Hepatocytes_infected"


metadata[SARS.CoV.2.Plus == "TRUE" & Major_compartments == "Hepatocytes"]$Label <- "HEP Inf"

metadata[SARS.CoV.2.Plus == "TRUE" & Major_compartments == "Hepatocytes"]$ClusterName <- "HEP Inf"

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Cell.id


LiverObject <- CreateSeuratObject(counts = custom_counts, project = "COVIDLiver", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 1,
                                  names.delim = "-", meta.data = metadata)

LiverObject[["CustomSeq"]] = CreateAssayObject(data = custom_counts)

DefaultAssay(LiverObject) = "CustomSeq"


cluster2_to_remove <- c("1_3", "1_8_1", "2_9", "3_5")
LiverObject <- LiverObject[,!LiverObject@meta.data$SubCluster %in% cluster2_to_remove]
LiverObject <- LiverObject[,!(LiverObject@meta.data$SARS.CoV.2.Plus=="TRUE" & LiverObject@meta.data$Major_compartments != "Hepatocytes")]

#Initialization
triMean_var <- 1 #Default trimMean value to produce fewer, still stronger interactions. Set this value to 0 for truncatedMean
min.cells <- 10 #Minimum number of cells per group


LiverObject.assay <- GetAssayData(LiverObject, assay = "CustomSeq", slot = "data")
labels <- LiverObject$Label
meta <- data.frame(group = labels, row.names = names(labels)) 

###################################################################################################
########################## MAIN code ##############################################################
###################################################################################################

############ Create CellChat object from Seurat ############

#future::plan("multiprocess", workers =4)
# do parallel

cellchat <- createCellChat(object = LiverObject.assay, meta = meta, group.by = "group", do.sparse = F)


cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")

#Set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "labels") 

#Show factor levels of the cell labels
levels(cellchat@idents)

#Number of cells in each cell group
groupSize <- as.numeric(table(cellchat@idents))

############ Set CellCharDB for ligand-receptor interactions ############

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Set CellChatDB in the object
CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

trim_var <- 0.2

chellchat.run <- run_cellChat(cellchat, meta, trim_var, result_dir)
