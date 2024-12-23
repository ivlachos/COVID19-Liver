library("Seurat")
library("data.table")
library("future")
library("slingshot")
library("scales")
library("harmony")

#### Initialization
annotation_label <- fread("/Data/Cell_populations.txt")
annotation_label_dbl <- annotation_label[regexpr("DBL.*",annotation_label$Label) == TRUE]
cluster2remove <- annotation_label_dbl$Cluster
Label2remove <- annotation_label_dbl$Label


#### Healthy
setnames(annotation_label,"Label","PredictedCluster")

name <- "Trajectory_Hep_Chol"
healthy <- readRDS("Toronto_Metadata.RDS")
healthy_bec <- readRDS("TorontoBECs_HarmonyEmbeddings.RDS")
healthy_hep <- readRDS("TorontoHepatocytes_HarmonyEmbeddings.RDS")
healthy_embeddings <- rbind(healthy_hep,healthy_bec)

healthy <- healthy[healthy$cell_ID %in% rownames(healthy_embeddings),]
healthy_embeddings <- healthy_embeddings[rownames(healthy_embeddings) %in% healthy$cell_ID,]

healthy_ids <- rownames(healthy)
healthy$id <- 1:nrow(healthy)
healthy <- merge(annotation_label, healthy, by = "PredictedCluster")
healthy <- healthy[order(id)]
healthy_hep_bec_labels <- healthy$Cluster

#### Covid
covid_hep_bec_metada <- readRDS("HepatocyteCholangiocytes_metadata.RDS")
covid_embeddings <- readRDS("HepatocyteCholangiocytes_HarmonyEmbeddings.RDS")


covid_hep_bec_metada$cell.id <- rownames(covid_hep_bec_metada)
covid_hep_bec_metada <- as.data.table(covid_hep_bec_metada)

covid_embeddings_new <- covid_embeddings
covid_embeddings_new <- as.data.frame(covid_embeddings_new)
covid_embeddings_new$cell.id <- rownames(covid_embeddings_new)
covid_embeddings_new <- as.data.table(covid_embeddings_new)
covid_embeddings_new <- merge(covid_embeddings_new,covid_hep_bec_metada, by = "cell.id")
covid_embeddings_new <- as.data.frame(covid_embeddings_new)
rownames(covid_embeddings_new) <- covid_embeddings_new$cell.id
covid_embeddings_new$cell.id <- NULL

embeddings <- rbind(healthy_embeddings,covid_embeddings_new[,1:50])
labels <- c(healthy_hep_bec_labels,covid_embeddings_new$SubCluster)

sds <- slingshot(embeddings[,1:20], clusterLabels = labels, omega=TRUE, start.clus = "0_2")

saveRDS(sds, file = paste0(result_dir,"/sds_",name,"_hepatocytes_Becs.rds"))

