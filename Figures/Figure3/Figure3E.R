
rm(list=ls())
options(stringsAsFactors = F)

library(Matrix)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)




# ==== Liver Data ====
# meta data
liver_metadata = readRDS(
  "COVID19LiverMetadata.RDS"
)


# Viral UMI Bootstrap test
viral_bootres = readRDS(
  "COVID19LiverSARS_Cov2_Enrichment.RDS"
)

ves_cluster = readRDS(
  "COVID19LiverViralEnrichment.RDS"
)


# ==== COVID+ Bar Plots ====
cluster_covid_plus_tab = with(
  viral_bootres,
  100*tapply(X=SARS.CoV.2.Plus,INDEX=ClusterName,FUN=mean)
)
ves_cluster$PCT = cluster_covid_plus_tab[ves_cluster$ClusterName]

ves_cluster$Label = ifelse(ves_cluster$FDR < 0.05,"*","")
ves_cluster$Label[ves_cluster$FDR < 0.01] = "**"
pbp_vescluster = ggplot(
  ves_cluster,
  aes(
    x=ClusterName,y=PCT,fill=qnES,
    label = Label
  )
) +
  geom_bar(stat="identity") +
  scale_fill_gradient2(
    low="blue",high="red",mid="grey86",
    name="Enrichment Score"
  ) +
  geom_text(
    hjust=0,
    vjust=0.85,
    angle=90,
    size=7.5
  ) +
  xlab("") +
  ylab("(%) SARS-CoV2 +") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle=90,
      vjust=0.5,
      hjust=1,
      size=15,
      color=ifelse(
        ves_cluster$FDR < 0.05,
        "red","black"
      ),
      face = ifelse(
        ves_cluster$FDR < 0.05,
        "bold.italic","plain"
      )
    ),
    legend.position = c(0.77,0.85),
    legend.margin = margin(10,10,10,10),
    legend.background = element_rect(color = "black")
  )


pbp_vescluster


