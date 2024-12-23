#!/usr/bin/env Rscript

rm(list=ls())
options(stringsAsFactors = FALSE)

library(optparse)
library(Seurat)
library(stackoverflow)
library(glmGamPoi)


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


cat("\nReading data ...\n")
# ===== Liver Data =====
liver_all = readRDS(
  "/n/scratch3/users/y/yhp2/LiverAll/CumulusCellbenderLiverAll.RDS"
)


cat("\nSplitting data ...\n")
# ===== Split Cells ====
Idents(liver_all) = "Batch"
sampleLst = CellsByIdentities(liver_all)
# shuffle cells
sampleLst = lapply(
  sampleLst,
  function(x){
    set.seed(667)
    sample(x,length(x))
  }
)


sampleChunkLst = lapply(
  sampleLst,
  function(x){
    chunk2(x,10)
  }
)

SampleChunk = lapply(
  1:10,
  function(ic){
    unlist(
      lapply(
        sampleChunkLst,
        function(x){
          x[[ic]]
        }
      )
      ,use.names = FALSE
    )
  }
)


ici = as.numeric( opt$CHUNK )
# ==== Split Genes ====
AllGenes = rownames(liver_all@assays$RNA@counts)
# shuffle genes 
set.seed(664)
AllGenes = sample(AllGenes,length(AllGenes))

GeneChunk = chunk2(AllGenes,10)


ChunkIndex = expand.grid(
  seq_along(SampleChunk),
  seq_along(GeneChunk)
)
colnames(ChunkIndex) = c("Cells","Genes")


cat("\nFitting GLM ...\n")
# ==== GLM ====
CurrentGenes = GeneChunk[[ChunkIndex$Genes[ici]]]
CurrentCells = SampleChunk[[ChunkIndex$Cells[ici]]]

glm_model00a = model.matrix(~  Batch + nCount_RNA,liver_all@meta.data[CurrentCells,])

liver_counts <- as.matrix(liver_all@assays$RNA@counts[,CurrentCells])

rm(liver_all)
gc()

tmp = rowSums(liver_counts)
liver_counts <- liver_counts[tmp > 5,]

cat(paste0("\n",sum(tmp <= 5)," genes removed\n"))

rm(tmp)

glmfit =  glm_gp(liver_counts, design = glm_model00a)

saveRDS(
  glmfit,
  paste0(
    "/n/scratch3/users/y/yhp2/LiverAll/GLM2/glm_gpFitv2_",
    ici,".RDS"
  )
)

