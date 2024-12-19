
rm(list=ls())
options(stringsAsFactors = FALSE)

library(rafalib)
library(tidyverse)
library(pheatmap)
library(parallel)
library(ggplot2)
library(cowplot)


# ===== Viral UMI Results ====
viral_umi_enrich = readRDS(
  "COVID19LiverSARS_Cov2_Enrichment.RDS"
)

# ==== Viral Enrichment Score (Donor) =====
donor_viralenrich = readRDS(
  "COVID19LiverDonorViralEnrichment.RDS"
)


# ==== Figure: Viral Enrichment Barplot =====
tmpdf = viral_umi_enrich %>% filter(!duplicated(DonorID)) %>% select(DonorID,MedicalCenter)
tmpdf$DonorID = factor(tmpdf$DonorID,levels=paste0("L",1:17))
tmpdf = tmpdf[order(tmpdf$DonorID),]

tmptab = tapply(
  X=viral_umi_enrich$SARS.CoV.2.Plus,INDEX=viral_umi_enrich$DonorID,FUN=function(x){100*mean(x)}
)

tmptab = tmptab[paste0("L",1:17)]


missing_donors = setdiff(paste0("L",1:17),donor_viralenrich$DonorID)

plotdf = rbind(
  donor_viralenrich,
  data.frame(
    DonorID = missing_donors,
    ES = NA,
    NES = NA,
    Pval = NA,
    FDR = NA,
    scNES = NA,
    scES=NA,
    qnES = NA
  )
)
plotdf$COVID.Pct = tmptab[plotdf$DonorID]
plotdf$DonorID = factor(plotdf$DonorID,levels=paste0("L",1:17))

plotdf$FDR[is.na(plotdf$FDR)]=1

ves_donorbp = ggplot(
  plotdf,
  aes(
    x=DonorID,y=COVID.Pct,fill=scNES,
    label = ifelse(FDR < 0.05,"*","")
  )
) +
  geom_bar(stat="identity") +
  scale_fill_gradient2(
    low="blue",high="red",mid="grey86",
    name="Enrichment Score"
  ) +
  geom_text(
    vjust=0.5,
    hjust = "center",
    angle=0,
    size=11
    
  ) +
  xlab("") +
  ylab("(%) SARS-CoV2 +") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      size=17,
      angle=90,
      vjust=0.5,
      hjust=1,
      color=ifelse(
        plotdf$FDR < 0.05,
        "red","black"
      ),
      face = ifelse(
        plotdf$FDR < 0.05,
        "bold.italic","plain"
      )
    ),
    legend.position = c(0.77,0.75),
    legend.margin = margin(10,10,10,10),
    legend.background = element_rect(color = "black")
  )

ves_donorbp



