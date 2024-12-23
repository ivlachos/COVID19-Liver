rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(rafalib)
library(ggrepel)
library(RColorBrewer)
library(DT)
library(rmarkdown)
library(kableExtra)
library(DirichletReg)
library(readxl)
library(emmeans)
library(lme4)



# ==== Cluster Manual Annotation ====
subcluster_annot = read_xlsx(
  path="/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/subcluster_annotation.xlsx",
  sheet=1
) %>% as.data.frame()

# cluster names used in manuscript
short_names = subcluster_annot$Label
names(short_names) = subcluster_annot$`Clustering v.1.1`

compartment_names = c(
  "0" = "Hepatocytes",
  "1" = "Immune",
  "2" = "Endothelial",
  "3" = "Mesenchymal",
  "4" = "BECs"
)



# ==== COVID Liver =====
liver_metadata = readRDS(
  "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/SelectLiversNoDoublets_Clusteringv1.1_DonorIDs.RDS"
)

liver_metadata$Batch2 = gsub("ncLab","",liver_metadata$Batch)
liver_metadata$Batch2 = gsub("Broad","",liver_metadata$Batch2)


# ===== Projection Results ====
NN20PCsNoDblProj = readRDS(
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "TorontoNN20PCsNoDblProj.RDS"
  )
)


# ===== Binomial GLMM ====
HelperBinomGLMM = function(x){
  BinomGLMByCompartment = lapply(
    compartment_names,
    function(ck){
      # samples x cell types counts
      tmpDF01 = liver_metadata[liver_metadata$Compartment == ck,c("Batch","SubCluster")]
      colnames(tmpDF01) = c("Sample","Cluster")
      
      eval(
        parse(
          text = paste0(
            "tmpDF02 = ",x,"[",x,'$PredictedCompartment == ck,c("sample","SubCluster")]'
          )
        )
      )
      colnames(tmpDF02) = c("Sample","Cluster")
      
      sample_annot = rbind(tmpDF01,tmpDF02)
      sample_annot$Condition = rep("Control",nrow(sample_annot))
      sample_annot$Condition[sample_annot$Sample %in% unique(liver_metadata$Batch)] = "COVID"
      
      # remove doublets
      clust_rm = names(short_names)[grepl("DBL",short_names)]
      sample_annot = sample_annot[!(sample_annot$Cluster %in% clust_rm),]
      
      
      # observed cluster frequencies
      clustertab_obs = table(sample_annot[,c("Cluster","Condition")]) %>% prop.table(margin=1)
      covid_specific = rownames(clustertab_obs)[clustertab_obs[,"Control"] == 0]
      
      # remove COVID specific clusters (model won't work well with these)
      df = sample_annot %>% 
        filter(!(Cluster %in% covid_specific)) %>%
        select(Condition,Sample,Cluster) %>%
        group_by(.dots = c("Condition","Sample","Cluster")) %>%
        summarise(Count=n())
      
      total_vec = table(sample_annot$Sample)
      
      df$Total = total_vec[df$Sample]
      df = df %>% mutate(Other = Total - Count)
      
      # Fit binomial GLM
      formula01 = cbind(Count,Other) ~ Cluster * Condition  +  Sample * Condition
      
      # make subject a random effect
      glmm01 = glmer(
        formula = cbind(Count,Other) ~ Cluster * Condition  +  (1|Sample), 
        family="binomial", 
        data=df,
        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
      )
      
      emm2 <- emmeans(glmm01, specs = revpairwise ~ Condition | Cluster)
      emm2$contrasts %>%
        summary(infer = TRUE, type = 'response') %>%
        rbind() %>%
        as.data.frame() -> c2_results
      
      resDF = c2_results %>% 
        # mutate(FDR = p.adjust(p=p.value,method="fdr")) %>%
        arrange(desc(odds.ratio)) %>%
        # filter(FDR < 0.001) %>%
        select(-c(df,asymp.LCL,asymp.UCL,null))
      
      resDF$Cluster = short_names[as.character(resDF$Cluster)]
      resDF
    }
  )
  
  BinomGLMByCompartment = do.call(rbind,BinomGLMByCompartment)
  BinomGLMByCompartment$FDR = p.adjust(p=BinomGLMByCompartment$p.value,method="fdr")
  rownames(BinomGLMByCompartment) = NULL
  BinomGLMByCompartment
}

NN20PCsNoDblBinomGLMM = HelperBinomGLMM(x = "NN20PCsNoDblProj")




