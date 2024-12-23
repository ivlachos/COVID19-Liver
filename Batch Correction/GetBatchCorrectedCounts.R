#!/usr/bin/env Rscript

# conda activate r_enc
.libPaths("/home/yered/miniconda3/envs/r_enc/lib/R/library")


rm(list=ls())
options(stringsAsFactors = FALSE)

library(optparse)
library(Seurat)
library(glmGamPoi)
library(sva)
library(edgeR)
# library(stackoverflow)



source("/data/work/Projects/BrScRNAseq/src/CustomCombatSeq.R")
option_list <- list(
  make_option(
    opt_str = c("--CHUNK"),
    action="store",
    type="integer",
    default=27,
    help="Chunk index",
    metavar = "integer"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# ==== Function ====
chunk2 = function (x, n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}


cat("\nReading data ...\n")
# ==== Liver Data =====
liver_select <- readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLivers/SelectLiversHarmonySeuratObjNoDoublets.RDS"
)

cat("\nSplitting data ...\n")
ici = as.numeric( opt$CHUNK )
# ==== Split Genes ====
AllGenes = rownames(liver_select@assays$RNA@counts)
# shuffle genes 
set.seed(664)
AllGenes = sample(AllGenes,length(AllGenes))

GeneChunk = chunk2(AllGenes,10)



# ==== GLM ====
cat("\nFitting GLM ...\n")
CurrentGenes = GeneChunk[[ ici ]]

batch_factor = liver_select@meta.data$Batch


liver_counts <- as.matrix(liver_select@assays$RNA@counts[CurrentGenes,])

rm(liver_select)
gc()


liver_adjusted <- ComBatCustomSeq(
  as.matrix(liver_counts), 
  batch=batch_factor, group=NULL,
  shrink=TRUE, shrink.disp=TRUE, gene.subset.n=ceiling(0.1*nrow(liver_counts))
)

cat("\nSaving Results ...\n")
saveRDS(
  liver_adjusted,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/",
    "SelectLiversNoDoublets/CustomSeq/Shrink/AllLivers_CustomSeq_Shrink",ici,".RDS"
  )
)

cat("\nDONE! >:( \n")

