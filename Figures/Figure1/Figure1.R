library("RColorBrewer")
library("ggplot2")
library("viridis")
library("scales")
library("ggrepel")
library("googledrive")
library("googlesheets4")
library("ComplexHeatmap")
library("pryr")
library("data.table")
library("dplyr")
library("ggsci")
library("cowplot")
library("grid")
library("gridExtra")
library("scater")
library("ggplotify")
library("gtable")

########## Initialization ##############

hexbinlst <- readRDS("/Data/COVID19LiverHexbinLst.RDS")
metadata <- readRDS("/Data/COVID19LiverMetadata.RDS")
umap_embeddings <- readRDS("/Data/COVID19LiverUMAPembeddingsLst.RDS")
logcounts <- readRDS("COVID19LiverCorrectedLogCounts_Fig1.RDS")
markers_majorClusters <- fread("/Data/major_markers.txt", header=T)
clustername.table <- setnames(fread("/Data/Cell_populations.txt"),c("SubCluster","Label","Alias","Major_compartments"))

seur_obj <- CreateSeuratObject(counts = logcounts, project = "COVIDLiver", assay = "RNA",
                               min.cells = 0, min.features = 0, names.field = 1,
                               names.delim = "-", meta.data = metadata)

seur_obj[["logcounts"]] = CreateAssayObject(data = logcounts)

seur_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_embeddings$All,key="UMAP_",assay = "logcounts")


#################### Add Major cluster Labels

compartments <- data.table("SubCluster"= metadata$SubCluster)

setDF(compartments)
rownames(compartments) <- rownames(metadata)
compartments$sample <- rownames(compartments)
compartments <- merge(compartments, setDF(clustername.table), by = "SubCluster")
rownames(compartments) <- compartments$sample
compartments$sample <- NULL

seur_obj <- AddMetaData(object = seur_obj, metadata = compartments)

compartments$Cell.id <- rownames(compartments)
##################### Figure 3a UMAP plot

umap_1 <- as.data.frame(umap_embeddings$`All`)
umap_1$Cell.id <- rownames(umap_1)
plot_df <- merge(umap_1, compartments, by = "Cell.id")
plot_df$Cell.id <- NULL

colnames(plot_df) <- c(paste0('UMAP', 1:2), 'SubCluster',"Label","Alias","Major_compartments")

plot_df$SubCluster <- as.factor(plot_df$SubCluster )
plot_df$Major_compartments <- as.factor(plot_df$Major_compartments )

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

plot_df2 <- plot_df %>% group_by(Major_compartments) %>% summarize_all(mean)


figureA <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) + ggtitle("n=80,808") + 
  geom_point(data = plot_df, pch = 21, color = 'black', size = 0.005,alpha = 0.4) +
  geom_point(aes(color = SubCluster), alpha = 0.4, size = 0.005, show.legend = T) +
   ggrepel::geom_text_repel(data = plot_df2,aes(label = Major_compartments), size=2.5, max.iter = 60) +
  scale_color_manual(values = cols, 
                     labels = paste0(levels(factor(plot_df$SubCluster)))) +
  
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=8, face='bold', hjust = 0.08),
        axis.title=element_text(size=8,color="black",face="bold"),
        legend.title = element_blank(),axis.text.x = element_text(size=8),axis.text.y = element_text(size=8),
        legend.text = element_text(colour="black", size = 6),
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


figureA <- figureA + annotate(geom="text", x=min_x + 5.15, y=min_y - 0.65,
                    color="black", label = statement1,size= 3, parse = TRUE)

figureA <- figureA + annotate(geom="text", x=min_x - 0.1, y=min_y + 5.7,
                    color="black", label = statement2,size= 3, parse = TRUE,
                    angle = 90)


############## Figure 1c MARKERs plots in UMAP

markers_majorClusters <- as.data.frame(t(markers_majorClusters))
markers_majorClusters$Lineage <- rownames(markers_majorClusters)
markers_majorClusters$Lineage <- gsub("\\..*","",markers_majorClusters$Lineage)
markers_majorClusters <- setnames(as.data.table(markers_majorClusters), c("Gene_rnaSeq","Gene","Cell Type"))

seur_obj@active.ident <- factor(seur_obj$Major_compartments)

cell_type_1 <- markers_majorClusters[`Cell Type` == "Hepatocytes"][Gene_rnaSeq %in% c("HNF4A","APOA1","TF")]$Gene_rnaSeq
cell_type_2 <- markers_majorClusters[`Cell Type` == "Immune"][Gene_rnaSeq %in% c("PTPRC","MARCO","THEMIS")]$Gene_rnaSeq
cell_type_3 <- markers_majorClusters[`Cell Type` == "Endothelial"][Gene_rnaSeq %in% c("FCGR2B","VWF","FLT1")]$Gene_rnaSeq
cell_type_4 <- markers_majorClusters[`Cell Type` == "Mesenchymal"]$Gene_rnaSeq
cell_type_5 <- markers_majorClusters[`Cell Type` == "BECs"][Gene_rnaSeq %in% c("CFTR","KRT7","KRT19")]$Gene_rnaSeq

type_list = c("cell_type_1","cell_type_2","cell_type_3","cell_type_4","cell_type_5")
cellEmbed = seur_obj@reductions$umap@cell.embeddings

plot <- list();
k <- 1
for(j in type_list){
  cellType=j
  for(i in get(cellType)){
   # cols <- c("grey" ,brewer.pal(9,"YlOrRd"))
    tryCatch(
      assay_data <- GetAssayData(object = seur_obj, assay= "logcounts")[i,],
      error = function(e) print(paste0("no ",i)),
      warning = function(w) print(paste0("no ",i)))
    df = data.frame(
      x=cellEmbed[, 1], 
      y=cellEmbed[, 2])
    
    
    df$`Log2(Expr)`=assay_data
   
    df$alpha <- 1
    df[df$`Log2(Expr)`!=0,"alpha"] <- 0.5
    df$size <- 0.005
    df[df$`Log2(Expr)`!=0,"size"] <- 2.5
    data<-df[order(df$`Log2(Expr)`, decreasing=FALSE),]
    
    
    name <- markers_majorClusters[Gene_rnaSeq == i]$Gene     
     plot[[k]] <- ggplot(data,aes(x=x, y=y, colour=`Log2(Expr)`),mar=c(0,0,3,0)) + 
      ggtitle(name) +
      geom_point(data = data, pch = 21, size = 0.0005,alpha = 0.3) + 
      ylab("UMAP2") + xlab("UMAP1") + 
       scale_color_gradientn(colours = colorRampPalette(c("darkblue","yellow"))(3), name = "Log(Expr)") + #scale_size(range=c(0,0.7))+
    
       theme(panel.border = element_blank(), panel.background = element_blank(),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             plot.title = element_text(size=6,color="black",face="bold"), 
             axis.title=element_blank(),
             legend.title = element_text(size=6,color="black",face="bold"),
             legend.text = element_text(colour="black", size = 5),
             legend.key = element_rect(fill = "transparent", size = 5,colour = "transparent"),
             legend.key.size = unit(0.2, "cm"), legend.direction="horizontal", 
             legend.key.width = unit(0.2,"cm"), legend.position = c(0.6,0.9999),
             axis.line = element_blank(),
             legend.background = element_blank(),
             legend.box.background = element_blank(),
             axis.text.x = element_blank(), axis.ticks.x = element_blank(),
             axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             plot.margin = unit(c(0,0,0,0), "cm")) 
     
     min_y <- min(cellEmbed[,2])
     min_x <- min(cellEmbed[,1])
     max_y <- max(cellEmbed[,2])
     max_x <- max(cellEmbed[,1])
     
     statement1 <- 'atop(bold("UMAP1"))'
     statement2 <- 'atop(bold("UMAP2"))'

     
     plot[[k]] <- plot[[k]] + coord_cartesian(xlim=c(min_x-0.5,max_x),ylim=c(min_y-0.5,max_y)) + 
       geom_segment(aes(y = min_y,
                        yend = min_y+1/6*(abs(max_y) + abs(min_y)),
                        x = min_x - 1,
                        xend = min_x - 1),
                    arrow = arrow(length = unit(0.1, "cm")), color="black")+
       geom_segment(aes(y = min_y,
                        yend = min_y,
                        x = min_x - 1,
                        xend = min_x+1/6*(abs(max_y) + abs(min_y)) - 1), color="black",
                    arrow = arrow(length = unit(0.1, "cm"))) +
       annotate(geom="text", x=min_x + 5.5, y=min_y - 2,
                color="black", label = statement1,size= 2, parse = TRUE) + 
       annotate(geom="text", x=min_x + 0.5, y=min_y + 7,color="black", label = statement2,size= 2, parse = TRUE,angle = 90)
     
     if(k == 1){
       
       plot[[k]] <- plot[[k]] + theme(plot.margin = unit(c(0,0,0,0), "cm")) 
       
     }
     
    if(!k %in% c(1,4,7,10,13)){
      
      plot[[k]] <- plot[[k]] + ylab("") 
      
    }
    if(!k %in% c(3,6,9,12,15)){
      
      plot[[k]] <- plot[[k]] + theme( #legend.position = "none",
           plot.margin = unit(c(0,0,0,0), "cm")) 
    } 
     
     
     k <- k+1

}}

title <- ggdraw() + 
  draw_label("Hepatocytes",fontface = 'bold',x = 0,hjust = 0, size=7) +
  theme(plot.margin =  unit(c(0,0,0,0.5), "cm"))

hepatocytes <- plot_grid(plotlist = list(plot[[1]],plot[[2]],plot[[3]]), 
                               ncol = 3, nrow=1)


figureB.hepatocytes <- plot_grid(title, hepatocytes, ncol = 1, rel_heights = c(0.1, 1))


title <- ggdraw() + 
  draw_label("Immune",fontface = 'bold',x = 0,hjust = 0, size=7) +
  theme(plot.margin = unit(c(0,0,0,0.5), "cm"))

Immune <- plot_grid(plotlist = list(plot[[4]],plot[[5]],plot[[6]]), 
                         ncol = 3, nrow=1)


figureB.Immune <- plot_grid(title, Immune, ncol = 1, rel_heights = c(0.1, 1))


title <- ggdraw() + 
  draw_label("Endothelial",fontface = 'bold',x = 0,hjust = 0, size=7) +
  theme(plot.margin = unit(c(0,0,0,0.5), "cm"))

Endothelial <- plot_grid(plotlist = list(plot[[7]],plot[[8]],plot[[9]]), 
                         ncol = 3, nrow=1)


figureB.Endothelial <- plot_grid(title, Endothelial, ncol = 1, rel_heights = c(0.1, 1))


title <- ggdraw() + 
  draw_label("Mesenchymal",fontface = 'bold',x = 0,hjust = 0, size=7) +
  theme(plot.margin =  unit(c(0,0,0,0.5), "cm"))


Mesenchymal <- plot_grid(plotlist = list(plot[[10]],plot[[11]],plot[[12]]), 
                         ncol = 3, nrow=1)

figureB.Mesenchymal <- plot_grid(title, Mesenchymal, ncol = 1, rel_heights = c(0.1, 1))


title <- ggdraw() + 
  draw_label("BECs",fontface = 'bold',x = 0,hjust = 0, size=7) +
  theme(plot.margin = margin(0, 0, 0, 7))


Cholangiocytes <- plot_grid(plotlist = list(plot[[13]],plot[[14]],plot[[15]]), 
                            ncol = 3, nrow=1)

figureB.Cholangiocytes <- plot_grid(title, Cholangiocytes, ncol = 1, rel_heights = c(0.1, 1))


figureB <- plot_grid(figureB.hepatocytes, figureB.Immune, figureB.Endothelial, figureB.Mesenchymal, figureB.Cholangiocytes,
                         ncol = 1, nrow=5, label_size = 10)


####################### Figure 1d, Heatmap

markers_majorClusters <- markers_majorClusters[!Gene_rnaSeq %in% c("KRT19","KRT7")]

SubClusterAvg = scuttle::summarizeAssayByGroup(
  x = seur_obj@assays$logcounts@data,
  ids = seur_obj$Major_compartments,
  statistics = "mean")

SubClusterAvg <- SubClusterAvg@assays@data$mean
SubClusterAvg <- SubClusterAvg[,colnames(SubClusterAvg) %in% seur_obj$Major_compartments]

annDF = data.frame(Cluster = c("BECs","Endothelial","Hepatocytes","Immune","Mesenchymal"))
rownames(annDF) = colnames(SubClusterAvg)

dt <- as.data.frame(SubClusterAvg[rownames(SubClusterAvg) %in% as.character(markers_majorClusters$Gene_rnaSeq),])
dt$Gene_rnaSeq <- rownames(dt)
dt <- merge(dt, markers_majorClusters, by = "Gene_rnaSeq")
dt <- setDF(dt)
rownames(dt) <- dt$Gene 
dt$Gene_rnaSeq <- NULL
dt$Gene <- NULL
dt$`Cell Type` <- NULL

annRow <- setDF(markers_majorClusters[,c("Cell Type"), with=F])
rownames(annRow) <- markers_majorClusters$Gene
markers_majorClusters$`Cell Type` <- factor(markers_majorClusters$`Cell Type`, levels = c("Hepatocytes", "Immune","Endothelial","Mesenchymal","BECs"))
rownames(annRow) <- factor(rownames(annRow), levels = rownames(annRow))

dt <- dt[rownames(annRow),]
dt <- na.omit(dt)

getPalette = colorRampPalette(brewer.pal(n = 8, name = "Spectral"))
mycolors <- getPalette(length(unique(annDF$Cluster)))
names(mycolors) <- unique(annDF$Cluster)
mycolors <- list(`Cell Type` = mycolors)

annRow$`Cell Type` <- factor(annRow$`Cell Type`, levels = c("Hepatocytes", "Immune","Endothelial","Mesenchymal","BECs"))

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



figureC <- ComplexHeatmap::pheatmap(mat = as.matrix(dt),cluster_rows=FALSE,cluster_cols =FALSE,fontsize = 7, 
                                    annotation_row = annRow,fontsize_row = 7, fontsize_col = 7,annotation_colors=mycolors,
                                    cellwidth=9, cellheight=9,name="Log(Expr)", color = colors)


png(filename = "Figure1B.png", res=300, height=4, width=4, units="in")

print(figureA)

dev.off()


png(filename = "Figure1C.png", res=300, height=8, width=5, units="in")

print(figureB)

dev.off()

png(filename = "Figure1D.png", res=300, height=4, width=4, units="in")

print(figureC)

dev.off()
