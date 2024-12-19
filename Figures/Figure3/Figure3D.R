rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(ggpubr)

# ==== Differential Abundance Results ====
ResultsDA = read_tsv(
  "COVIDDifferentialAbundance.tsv"
)


# ==== Barplot ====
tmpDF = reshape2::melt(
  ResultsDA[,c("Control","COVID","Cluster","Compartment")]
)
# scale relative abundances to (%)
tmpDF$value = 100*tmpDF$value
# set cluster names as facto to preserve order
tmpDF$Cluster = factor(tmpDF$Cluster,levels=ResultsDA$Cluster)
# mark significant clusters
tmpDF$Significant = tmpDF$Cluster %in% ResultsDA$Cluster[ResultsDA$GLMM.FDR < 0.05]
tmpDF$lab = ifelse(tmpDF$Significant,"*","")
tmpDF$lab[tmpDF$Cluster %in% ResultsDA$Cluster[ResultsDA$GLMM.FDR < 0.01]] = "**"

# barplot
ggplot(
  tmpDF,
  aes(x=Cluster,y=value,fill=variable)
) +
  geom_bar(stat="identity", position=position_dodge(1)) +
  labs(fill="") +
  ylab("Relative Abundance (%)") +
  xlab("") +
  geom_text(
    data=tmpDF,
    aes(y=value,x=Cluster,label=lab,group=variable),
    hjust=0,
    vjust=0.8,
    angle=90,
    position = position_dodge(1),
    inherit.aes = FALSE,
    size=6.5
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle = 90, 
      vjust = 0.5, 
      hjust = 1,
      color=ifelse(
        tmpDF$Significant,
        "red","black"
      ),
      face = ifelse(
        tmpDF$Significant,
        "bold.italic","plain"
      )
    ),
    legend.position = "top",
    legend.justification = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,-30,0)
  ) 
