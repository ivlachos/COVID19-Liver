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
library("cowplot")
library("scater")
library("ggplotify")

############# Initialization

source("pheatmap_fixed.R")


hexbinlst <- readRDS("/Data/COVID19LiverHexbinLst.RDS")
metadata <- readRDS("/Data/COVID19LiverMetadata.RDS")
umap_embeddings <- readRDS("/Data/COVID19LiverUMAPembeddingsLst.RDS")
summaries <- readRDS("/Data/COVID19Liver_SubClusterLimmaSummaries.RDS")

immune <- fread("/Data/markers_immune.txt")
endothelial <- fread("/Data/markers_endothelial.txt")
mesenchymal <- fread("/Data/markers_mesenchymal.txt")
clustername.table <- setnames(fread("/Data/Cell_populations.txt"),c("SubCluster","Label","Alias","Major_compartments"))


############# Add Major cluster Labels

names <- colnames(summaries)

compartments <- data.table("SubCluster"= metadata$SubCluster)

setDF(compartments)
rownames(compartments) <- rownames(metadata)
compartments$sample <- rownames(compartments)
compartments <- merge(compartments, setDF(clustername.table), by = "SubCluster")
rownames(compartments) <- compartments$sample
compartments$sample <- NULL

compartments$Cell.id <- rownames(compartments)


############# UMAP plots

umap_1 <- as.data.frame(umap_embeddings$`1`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$Label <- as.factor(plot_df$Label )

# Make color scale
cols <- c("deeppink3","brown2","brown","darkcyan","darkgoldenrod",
          "darkgoldenrod1","darkgoldenrod4","burlywood3","blueviolet","cornflowerblue","chocolate1",
          "chocolate4","chartreuse4","aquamarine1","chartreuse","black")

plot_df2 <- plot_df %>% group_by(Label) %>% summarize_all(mean)

figureA <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = plot_df, pch = 21, color = 'black', size = 0.01,alpha = 0.2) +
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
        axis.line = element_line(size = 0.5, colour = "black"))

figureA <- figureA +  theme(axis.line = element_blank(), axis.title=element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(), axis.ticks.y = element_blank())

min_y <- min(umap_1[,2])
min_x <- min(umap_1[,1])
max_y <- max(umap_1[,2])
max_x <- max(umap_1[,1])
# 
figureA <- figureA + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y))

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

figureA <- figureA + annotate(geom="text", x=min_x + 3.7, y=min_y - 0.5,
                              color="black", label = statement1,size= 3, parse = TRUE)

figureA <- figureA + annotate(geom="text", x=min_x - 0.4, y=min_y + 4.5,
                              color="black", label = statement2,size= 3, parse = TRUE,
                              angle = 90)



############# Heatmap

immune <- as.data.frame(t(immune))
immune$Lineage <- rownames(immune)
immune$Lineage <- gsub("\\..*","",immune$Lineage)
immune <- setnames(as.data.table(immune), c("Gene_rnaSeq","Gene","SubCluster"))
immune$`SubCluster` <- gsub("X","", immune$`SubCluster`)
immune <- merge(immune, clustername.table, by = "SubCluster")
immune <- unique(immune)


colnames.immune <- names[regexpr("1_.*.Avg", names)==TRUE]
summaries.immune <- summaries[, colnames(summaries) %in% colnames.immune]
summaries.immune <- summaries.immune[rownames(summaries.immune) %in% immune$Gene_rnaSeq,]
colnames(summaries.immune) <- gsub(".Avg","",colnames(summaries.immune))
names.hep <- setnames(data.table(colnames(summaries.immune)),"SubCluster")
names.hep <- merge(names.hep, clustername.table, by="SubCluster")
colnames(summaries.immune) <- names.hep$Label

annDF = data.frame(Cluster = sapply(strsplit(colnames(summaries.immune),"_"),head,1))
rownames(annDF) = colnames(summaries.immune)

dt <- as.data.frame(summaries.immune[rownames(summaries.immune) %in% as.character(immune$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, immune, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Label` <- NULL

annRow <- setDF(immune[,c("Major_compartments"), with=F])
rownames(annRow) <- immune$Gene
immune$`Label` <- as.factor(immune$`Label`)
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


figureB <- pheatmap_fixed(mat = t(dt),cluster_rows=F,cluster_cols =T,legend = T,
         fontsize_row = 5, fontsize_col =5,scale="column",
         cellwidth=4.5, cellheight=7,color =inferno(10), name="Log(Expr)",fontsize = 7,
         heatmap_legend_param = list(title_gp=gpar(fontsize=7,fontface="bold"),
                                     labels_gp = gpar(fontsize = 7)))


############# Figure 5b UMAP plot
umap_1 <- as.data.frame(umap_embeddings$`2`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$Label <- as.factor(plot_df$Label )

# Make color scale
getPalette = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))
cols <- getPalette(length(unique(plot_df$Label)))

plot_df2 <- plot_df %>% group_by(Label) %>% summarize_all(mean)

figureC <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
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
        axis.line = element_line(size = 0.5, colour = "black"))



figureC <- figureC +  theme(axis.line = element_blank(), axis.title=element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(), axis.ticks.y = element_blank())

min_y <- min(umap_1[,2])
min_x <- min(umap_1[,1])
max_y <- max(umap_1[,2])
max_x <- max(umap_1[,1])
# 
figureC <- figureC + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y))

figureC <-figureC + geom_segment(aes(y = min_y,
                                     yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                                     x = min_x - 1,
                                     xend = min_x - 1),
                                 arrow = arrow(length = unit(0.1, "cm")))

figureC <- figureC + geom_segment(aes(y = min_y,
                                      yend = min_y,
                                      x = min_x - 1,
                                      xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1),
                                  arrow = arrow(length = unit(0.1, "cm")))


statement1 <- 'atop(bold("UMAP1"))'
statement2 <- 'atop(bold("UMAP2"))'
statement3 <- 'atop(bold("'
statement4 <- '"))'
repel_face <- "bold"

figureC <- figureC +annotate(geom="text", x=min_x + 3, y=min_y-0.3,
                             color="black", label = statement1,size= 3, parse = TRUE)

figureC <- figureC + annotate(geom="text", x=min_x + 0.3, y=min_y + 3,
                              color="black", label = statement2,size= 3, parse = TRUE,
                              angle = 90)

############# Heatmap

endothelial <- as.data.frame(t(endothelial))
endothelial$Lineage <- rownames(endothelial)
endothelial$Lineage <- gsub("\\..*","",endothelial$Lineage)
endothelial <- setnames(as.data.table(endothelial), c("Gene_rnaSeq","Gene","SubCluster"))
endothelial$`SubCluster` <- gsub("X","", endothelial$`SubCluster`)
endothelial <- merge(endothelial, clustername.table, by = "SubCluster")
endothelial <- unique(endothelial)



colnames.endothelial <- names[regexpr("2_.*.Avg", names)==TRUE]
summaries.endothelial <- summaries[, colnames(summaries) %in% colnames.endothelial]
summaries.endothelial <- summaries.endothelial[rownames(summaries.endothelial) %in% endothelial$Gene_rnaSeq,]
colnames(summaries.endothelial) <- gsub(".Avg","",colnames(summaries.endothelial))
names.hep <- setnames(data.table(colnames(summaries.endothelial)),"SubCluster")
names.hep <- merge(names.hep, clustername.table, by="SubCluster")
colnames(summaries.endothelial) <- names.hep$Label



annDF = data.frame(Cluster = sapply(strsplit(colnames(summaries.endothelial),"_"),head,1))
rownames(annDF) = colnames(summaries.endothelial)

dt <- as.data.frame(summaries.endothelial[rownames(summaries.endothelial) %in% as.character(endothelial$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, endothelial, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Label` <- NULL

annRow <- setDF(endothelial[,c("Major_compartments"), with=F])
rownames(annRow) <- endothelial$Gene
endothelial$`Label` <- as.factor(endothelial$`Label`)
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

figureD <- pheatmap_fixed(mat = as.matrix(dt),cluster_rows=F,cluster_cols =T,scale="row",
                    fontsize_row = 4, fontsize_col =4,name = "Log(Expr)",
                    cellwidth=7, cellheight=3.5,color = inferno(10),fontsize = 5,
                    heatmap_legend_param = list(title_gp=gpar(fontsize=5,fontface="bold"),
                                                labels_gp = gpar(fontsize = 5)))



############# Figure 5c
umap_1 <- as.data.frame(umap_embeddings$`3`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$Label <- as.factor(plot_df$Label )

# Make color scale
getPalette = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))
cols <- getPalette(length(unique(plot_df$Label)))

plot_df2 <- plot_df %>% group_by(Label) %>% summarize_all(mean)

figureE <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(data = plot_df, pch = 21, color = 'black', size = 0.05,alpha = 0.4) +
  geom_point(aes(color = Label), alpha = 0.4, size = 0.05, show.legend = T) +
  ggrepel::geom_text_repel(data = plot_df2,aes(label = Label), size=2.7, max.iter = 60) +
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
        axis.line = element_line(size = 0.5, colour = "black"),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


figureE <- figureE +  theme(axis.line = element_blank(), axis.title=element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            axis.text.y = element_blank(), axis.ticks.y = element_blank())

min_y <- min(umap_1[,2])
min_x <- min(umap_1[,1])
max_y <- max(umap_1[,2])
max_x <- max(umap_1[,1])
# 
figureE<- figureE + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y))

figureE <-figureE + geom_segment(aes(y = min_y,
                                     yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                                     x = min_x - 1,
                                     xend = min_x - 1),
                                 arrow = arrow(length = unit(0.1, "cm")))

figureE <- figureE + geom_segment(aes(y = min_y,
                                      yend = min_y,
                                      x = min_x - 1,
                                      xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1),
                                  arrow = arrow(length = unit(0.1, "cm")))


statement1 <- 'atop(bold("UMAP1"))'
statement2 <- 'atop(bold("UMAP2"))'
statement3 <- 'atop(bold("'
statement4 <- '"))'
repel_face <- "bold"

figureE <- figureE + annotate(geom="text", x=min_x + 5, y=min_y-0.4,
                              color="black", label = statement1,size= 3, parse = TRUE)

figureE <- figureE + annotate(geom="text", x=min_x+0.5, y=min_y + 4.3,
                              color="black", label = statement2,size= 3, parse = TRUE,
                              angle = 90)


############# Heatmap

mesenchymal <- as.data.frame(t(mesenchymal))
mesenchymal$Lineage <- rownames(mesenchymal)
mesenchymal$Lineage <- gsub("\\..*","",mesenchymal$Lineage)
mesenchymal <- setnames(as.data.table(mesenchymal), c("Gene_rnaSeq","Gene","SubCluster"))
mesenchymal$`SubCluster` <- gsub("X","", mesenchymal$`SubCluster`)
mesenchymal <- merge(mesenchymal, clustername.table, by = "SubCluster")
mesenchymal <- unique(mesenchymal)



colnames.Stellate <- names[regexpr("3_.*.Avg", names)==TRUE]
summaries.Stellate <- summaries[, colnames(summaries) %in% colnames.Stellate]
summaries.Stellate <- summaries.Stellate[rownames(summaries.Stellate) %in% mesenchymal$Gene_rnaSeq,]
colnames(summaries.Stellate) <- gsub(".Avg","",colnames(summaries.Stellate))
names.hep <- setnames(data.table(colnames(summaries.Stellate)),"SubCluster")
names.hep <- merge(names.hep, clustername.table, by="SubCluster")
colnames(summaries.Stellate) <- names.hep$Label


annDF = data.frame(Cluster = sapply(strsplit(colnames(summaries.Stellate),"_"),head,1))
rownames(annDF) = colnames(summaries.Stellate)

dt <- as.data.frame(summaries.Stellate[rownames(summaries.Stellate) %in% as.character(mesenchymal$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, mesenchymal, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Label` <- NULL

annRow <- setDF(mesenchymal[,c("Major_compartments"), with=F])
rownames(annRow) <- mesenchymal$Gene
mesenchymal$`Label` <- as.factor(mesenchymal$`Label`)
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


figureF <-pheatmap_fixed(mat = as.matrix(dt),cluster_rows=F,cluster_cols =T,scale="row",
                    fontsize_row = 4, fontsize_col =4,name = "Log(Expr)",
                    cellwidth=9, cellheight=3.5,color = inferno(10),fontsize = 6,
                    heatmap_legend_param = list(title_gp=gpar(fontsize=6,fontface="bold"),
                                                labels_gp = gpar(fontsize = 6)))




png("Figure_5A_1.png", res=300,  height=4, width=4, units="in")

print(figureA)

dev.off()


png(filename = "Figure_5A_2.png", res=300,  height=4, width=7, units="in")

print(figureB)

dev.off()


png(filename = "Figure_5B_1.png", res=300,  height=4, width=2.1, units="in")

print(figureC)

dev.off()

png(filename = "Figure_5B_2.png", res=300,  height=5, width=3, units="in")

print(figureD)

dev.off()


png(filename = "Figure_5C_1.png", res=300,  height=5, width=2.1, units="in")

print(figureE)

dev.off()

png(filename = "Figure_5C_2.png", res=300,  height=5, width=5, units="in")

print(figureF)

dev.off()
