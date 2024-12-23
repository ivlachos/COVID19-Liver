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
##################### Figure 3a UMAP plot

umap_1 <- as.data.frame(umap_embeddings$`0`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$Label <- as.factor(plot_df$Label )

# Make color scale
getPalette = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))
cols <- getPalette(length(unique(plot_df$SubCluster)))

plot_df2 <- plot_df %>% group_by(Label) %>% summarize_all(mean)

figureA <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
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


figureA <- figureA +  theme(axis.line = element_blank(), axis.title=element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(), axis.ticks.y = element_blank())

min_y <- min(umap_1[,2])
min_x <- min(umap_1[,1])
max_y <- max(umap_1[,2])
max_x <- max(umap_1[,1])
# 
figureA <- figureA + coord_cartesian(xlim=c(min_x-0.9,max_x),ylim=c(min_y-0.9,max_y))

figureA <-figureA + geom_segment(aes(y = min_y,
                                     yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                                     x = min_x - 1,
                                     xend = min_x - 1),
                                 arrow = arrow(length = unit(0.1, "cm")))

figureA <- figureA + geom_segment(aes(y = min_y,
                                      yend = min_y,
                                      x = min_x - 1,
                                      xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1),
                                  arrow = arrow(length = unit(0.1, "cm")))


statement1 <- 'atop(bold("UMAP1"))'
statement2 <- 'atop(bold("UMAP2"))'
statement3 <- 'atop(bold("'
statement4 <- '"))'
repel_face <- "bold"

figureA <- figureA + annotate(geom="text", x=min_x + 1.4, y=min_y - 0.65,
                              color="black", label = statement1,size= 3, parse = TRUE)

figureA <- figureA + annotate(geom="text", x=min_x - 0.9, y=min_y + 2.1,
                              color="black", label = statement2,size= 3, parse = TRUE,
                              angle = 90)


########################## Heatmap plots

cholangocytes <- as.data.frame(t(cholangocytes))
cholangocytes$Lineage <- rownames(cholangocytes)
cholangocytes$Lineage <- gsub("\\..*","",cholangocytes$Lineage)
cholangocytes <- setnames(as.data.table(cholangocytes), c("Gene_rnaSeq","Gene","Cell Type"))
cholangocytes$`Cell Type` <- gsub("X","", cholangocytes$`Cell Type`)
setnames(clustername.table,"SubCluster","Cell Type")
cholangocytes <- merge(cholangocytes, clustername.table, by = "Cell Type")
cholangocytes <- unique(cholangocytes)

hepatocytes <- as.data.frame(t(hepatocytes))
hepatocytes$Lineage <- rownames(hepatocytes)
hepatocytes$Lineage <- gsub("\\..*","",hepatocytes$Lineage)
hepatocytes <- setnames(as.data.table(hepatocytes), c("Gene_rnaSeq","Gene","Cell Type"))
hepatocytes$`Cell Type` <- gsub("X","", hepatocytes$`Cell Type`)
hepatocytes <- merge(hepatocytes, clustername.table, by = "Cell Type")
hepatocytes <- unique(hepatocytes)

hepatocytes <- rbind(hepatocytes,cholangocytes)
hepatocytes <- hepatocytes[hepatocytes[, .I[sample(.N, 1)], by=Gene_rnaSeq]$V1]
names <- colnames(summaries)

colnames.hepatocytes <- names[regexpr("0_.*.Avg|4_.*.Avg", names)==TRUE]
summaries.hepatocytes <- summaries[, colnames(summaries) %in% colnames.hepatocytes]
summaries.hepatocytes <- summaries.hepatocytes[rownames(summaries.hepatocytes) %in% hepatocytes$Gene_rnaSeq,]
colnames(summaries.hepatocytes) <- gsub(".Avg","",colnames(summaries.hepatocytes))
names.hep <- setnames(data.table(colnames(summaries.hepatocytes)),"Cell Type")
names.hep <- merge(names.hep, clustername.table, by="Cell Type")
colnames(summaries.hepatocytes) <- names.hep$Label



colors <-  c("#170C3AFF","#420A68FF","#6B186EFF","#932667FF",
             "#BB3754FF","#000004FF" ,"#DD513AFF","#F3771AFF","#FCA50AFF","#F6D645FF","#FCFFA4FF")

annDF = data.frame(Cluster = sapply(strsplit(colnames(summaries.hepatocytes),"_"),head,1))
rownames(annDF) = colnames(summaries.hepatocytes)

dt <- as.data.frame(summaries.hepatocytes[rownames(summaries.hepatocytes) %in% as.character(hepatocytes$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, hepatocytes, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Label` <- NULL

annRow <- setDT(hepatocytes[,c("Major_compartments"), with=F])
annRow <- annRow[order(-annRow$Major_compartments)]
annRow <- setDF(annRow)
rownames(annRow) <- hepatocytes$Gene
hepatocytes$`Label` <- as.factor(hepatocytes$`Label`)
rownames(annRow) <- factor(rownames(annRow), levels = rownames(annRow))

dt <- dt[rownames(annRow),]
dt <- na.omit(dt)
dt$`Cell Type` <- NULL
dt$SubCluster <- NULL
dt$Major_compartments <- NULL
dt$Alias <- NULL

setnames(annRow, "Major_compartments", "Cell Type")

getPalette = colorRampPalette(brewer.pal(n = 8, name = "Spectral"))
mycolors <- getPalette(length(unique(annRow$`Cell Type`)))
names(mycolors) <- unique(annRow$`Cell Type`)
mycolors <- list(`Cell Type` = mycolors)


scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}
dt <- scale_mat(dt, "row")

colors <-  c("#000004FF","#280B54FF","#65156EFF","#9F2A63FF","#D44842FF","#F57D15FF","#FAC127FF","#FCFFA4FF")


figureC <- pheatmap_fixed(mat = as.matrix(dt),cluster_rows=F,cluster_cols =T,
                          annotation_row = annRow,fontsize_row = 6, fontsize_col =7,annotation_colors=mycolors,
                          cellwidth=12, cellheight=5,color =colors,name="Log(Expr)",fontsize = 10,
                          heatmap_legend_param = list(title_gp=gpar(fontsize=10,fontface="bold"),
                                                      labels_gp = gpar(fontsize = 10)))


############## Save Figures


png(filename = "Figure_3_A.png",height=4, width=5,res=300, units = "in")

print(figureA)

dev.off()


png(filename = "Figure_3_B.png",height=12, width=4,res=300, units = "in")

print(figureC)

dev.off()


