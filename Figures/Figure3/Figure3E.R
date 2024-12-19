rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(ggsignif)

# ==== Differential Abundance Results ====
DifferentialAbundanceResults = readRDS(
  "DifferentialAbundanceResults.RDS"
)


# ==== Plots ====
tmpDF = reshape2::melt(
  DifferentialAbundanceResults[,c("Control","COVID","Cluster","Compartment")]
)

tmpDF$value = 100*tmpDF$value
tmpDF$Significant = tmpDF$Cluster %in% DifferentialAbundanceResults$Cluster[DifferentialAbundanceResults$GLMM.FDR < 0.05]
tmpDF$Cluster = factor(tmpDF$Cluster,levels=DifferentialAbundanceResults$Cluster)
tmpDF$lab = ifelse(tmpDF$Significant,"*","")
tmpDF$lab[tmpDF$Cluster %in% DifferentialAbundanceResults$Cluster[DifferentialAbundanceResults$GLMM.FDR < 0.01]] = "**"

ggplot(
  tmpDF,
  aes(x=Cluster,y=value,fill=variable)
) +
  geom_bar(stat="identity", position=position_dodge2(1)) +
  scale_fill_manual(
    values = c(
      "Control" = "#00BFC4",
      "COVID" = "#F8766D"
    )
  ) + 
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

clusttemp = tmpDF %>%
  filter(variable == "Control",Significant == TRUE) %>%
  select(Cluster) %>%
  unlist()

xtemp = tmpDF %>%
  filter(variable == "Control") %>%
  select(Significant) %>%
  unlist() %>%
  which()

ytemp = apply(
  cbind(
    tmpDF %>%
      filter(variable == "Control",Significant == TRUE) %>%
      select(value) %>%
      unlist() ,
    tmpDF %>%
      filter(variable == "COVID",Significant == TRUE) %>%
      select(value) %>%
      unlist() 
  ),1,function(x){max(x) + 0.5}
)

labtemp = tmpDF %>%
  filter(variable == "Control",Significant == TRUE) %>%
  select(lab) %>%
  unlist()

signifdat = data.frame(
  Cluster = clusttemp,
  x = xtemp - 0.3,
  xend = xtemp + 0.3,
  y = ytemp,
  annotation = labtemp
)

tmp_ind = tmpDF %>%
  filter(variable == "COVID",Significant == TRUE) %>%
  select(value) %>%
  unlist() >
  tmpDF %>%
  filter(variable == "Control",Significant == TRUE) %>%
  select(value) %>%
  unlist() 


covid_up = (tmpDF %>%
              filter(variable == "COVID",Significant == TRUE) %>%
              select(Cluster) %>%
              unlist() %>% 
              as.character()
)[tmp_ind]

covid_down = (tmpDF %>%
                filter(variable == "COVID",Significant == TRUE) %>%
                select(Cluster) %>%
                unlist() %>% 
                as.character()
)[!tmp_ind]


# Covid specific clusters
covid_specific = DifferentialAbundanceResults %>%
  filter(
    Control == 0, COVID > 0
  ) %>%
  select(Cluster) %>%
  unlist()

# %>%
#   select(Cluster) %>%
#   unlist() %>% 
#   as.character()


tmpDF$xlabColor = "black"
tmpDF$xlabColor[tmpDF$Cluster %in% covid_up] = "#F8766D"
tmpDF$xlabColor[tmpDF$Cluster %in% covid_down] = "#00BFC4"
tmpDF$xlabColor[tmpDF$Cluster %in% covid_specific] = "red4"


pbar_base = ggplot(
  tmpDF,
  aes(x=Cluster,y=value,fill=variable)
) +
  geom_bar(stat="identity", position=position_dodge2(1)) +
  scale_fill_manual(
    values = c(
      "Control" = "#00BFC4",
      "COVID" = "#F8766D"
    )
  ) + 
  labs(fill="") +
  ylab("Relative Abundance (%)") +
  xlab("") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color=tmpDF$xlabColor,
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


for(ic in 1:nrow(signifdat)){
  pbar_base = pbar_base +
    geom_signif(
      stat = "identity",
      data = signifdat[ic,,drop=F],
      aes(
        x = x,
        xend = xend,
        y = y,
        yend = y,
        annotation = annotation
      ),
      inherit.aes = FALSE
    ) 
}

pbar_base 

