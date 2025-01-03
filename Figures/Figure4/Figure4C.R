library("RColorBrewer")
library("ggplot2")
library("viridis")
library("scales")
library("ggrepel")
library("googledrive")
library("googlesheets4")
library("pheatmap")
library("pryr")
library("data.table")
library("dplyr")
library("ggsci")
library("cowplot")
library("grid")
library("gridExtra")
library("scater")
library("ggplotify")
library("ggtree")
library("ComplexHeatmap")
library("circlize")


############ Initialization

source("pheatmap_fixed.R")

hexbinlst <- readRDS("/Data/COVID19LiverHexbinLst.RDS")
metadata <- readRDS("/Data/COVID19LiverMetadata.RDS")
umap_embeddings <- readRDS("/Data/COVID19LiverUMAPembeddingsLst.RDS")
clustername.table <- setnames(fread("/Data/Cell_populations.txt"),c("SubCluster","Label","Alias","Major_compartments"))
cholangocytes <- fread("/Data/markers_cholangiocytes.txt")
hepatocytes <- fread("/Data/markers_hepatocytes.txt")
summaries <- readRDS("/Data/COVID19Liver_SubClusterLimmaSummaries.RDS")


############ Add Major cluster Labels

compartments <- data.table("SubCluster"= metadata$SubCluster)

setDF(compartments)
rownames(compartments) <- rownames(metadata)
compartments$sample <- rownames(compartments)
compartments <- merge(compartments, setDF(clustername.table), by = "SubCluster")
rownames(compartments) <- compartments$sample
compartments$sample <- NULL

compartments$Cell.id <- rownames(compartments)

##################### UMAP plot

umap_1 <- as.data.frame(umap_embeddings$`4`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$Label <- as.factor(plot_df$Label )


# Make color scale
getPalette = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))
cols <- getPalette(length(unique(plot_df$Label)))

plot_df2 <- plot_df %>% group_by(Label) %>% summarize_all(mean)

figureB <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = plot_df, pch = 21, color = 'black', size = 0.05,alpha = 0.4) +
  geom_point(aes(color = Label), alpha = 0.4, size = 0.05, show.legend = T) +
  ggrepel::geom_text_repel(data = plot_df2,aes(label = Label), size=3, max.iter = 60) +
  scale_color_manual(values = cols, 
                     labels = paste0(levels(factor(plot_df$Label)))) +
  
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_blank(), 
        axis.title=element_text(size=8,color="black",face="bold"),
        legend.title = element_blank(),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        legend.text = element_text(colour="black", size = 7),
        legend.key = element_rect(fill = "transparent", size = 2,colour = "transparent"),
        legend.key.size = unit(1, "cm"), legend.position = "none", 
        axis.line = element_line(size = 0.5, colour = "black"), plot.margin = unit(c(0,0.5,0,0), "cm"))


figureB <- figureB +  theme(axis.line = element_blank(), axis.title=element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(), axis.ticks.y = element_blank())

min_y <- min(umap_1[,2])
min_x <- min(umap_1[,1])
max_y <- max(umap_1[,2])
max_x <- max(umap_1[,1])
# 
figureB <- figureB + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y))

figureB <-figureB + geom_segment(aes(y = min_y,
                                     yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                                     x = min_x - 1,
                                     xend = min_x - 1),
                                 arrow = arrow(length = unit(0.1, "cm")))

figureB <- figureB + geom_segment(aes(y = min_y,
                                      yend = min_y,
                                      x = min_x - 1,
                                      xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1),
                                  arrow = arrow(length = unit(0.1, "cm")))


statement1 <- 'atop(bold("UMAP1"))'
statement2 <- 'atop(bold("UMAP2"))'
statement3 <- 'atop(bold("'
statement4 <- '"))'
repel_face <- "bold"


figureB <- figureB + annotate(geom="text", x=min_x + 1.4, y=min_y - 0.35,
                              color="black", label = statement1,size= 3, parse = TRUE)

figureB <- figureB + annotate(geom="text", x=min_x - 0.6, y=min_y + 2.5,
                              color="black", label = statement2,size= 3, parse = TRUE,
                              angle = 90)

############## Save Figures

png(filename = "Figure_4C.png",height=4, width=5,res=300, units = "in")

print(figureB)

dev.off()

